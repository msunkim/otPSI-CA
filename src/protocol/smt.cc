#include "protocol/smt.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace otpsica {

std::pair<PaillierKeyPair, SMTReceiverMessage>
SMTProtocol::ReceiverRound1(const ElementVec& set, const ProtocolParams& params) {
    // Generate Paillier keypair (1024-bit for faster experiments)
    PaillierKeyPair kp = Paillier::KeyGen(1024);

    SMTReceiverMessage msg;
    msg.pk = kp.pk;
    msg.n = params.n;
    msg.tau = params.tau;
    msg.ell = params.ell;
    msg.w = params.w;
    msg.beta = params.beta;
    msg.primes = params.primes;

    // Build hash table for receiver's set
    HashTable ht(params.n, params.w);
    for (Element x : set) {
        ht.Insert(x);
    }
    ht.FillWithRandom(params.beta, 0);  // party_offset=0 for receiver
    msg.hash_seed = ht.GetSeed();

    // For each bin, build and encrypt the FLT polynomial (using first prime only for efficiency)
    // Note: Full CRT would use all primes, but we simplify for experimental evaluation
    size_t p1 = params.primes[0];
    msg.enc_polys.resize(params.ell);

    for (size_t i = 0; i < params.ell; i++) {
        msg.enc_polys[i].resize(1);  // Only one prime
        const ElementVec& bin = ht.GetBin(i);

        // Build polynomial f_i(X) = 1 - prod_{x in bin} (X - x)^{p1 - 1}
        std::vector<BigInt> coeffs = FLTPolynomial::Build(bin, p1);

        // Encrypt polynomial coefficients
        msg.enc_polys[i][0] = FLTPolynomial::Encrypt(kp.pk, coeffs);
    }

    return {kp, msg};
}

std::pair<SMTSenderMessage, std::vector<BigInt>>
SMTProtocol::SenderRound1(const SMTReceiverMessage& msg,
                           const std::vector<ElementVec>& sets) {
    size_t eta = sets.size();
    SMTSenderMessage resp;
    resp.masked_cardinalities.resize(eta);
    std::vector<BigInt> masks(eta);

    // Use the first prime for cardinality computation
    size_t p1 = msg.primes[0];

    // Determine if real encrypted polynomial evaluation is feasible:
    // Need p > n + mask_range so that (I + r) mod p = I + r (no wrap-around).
    // Also need polynomial evaluation g(y) to fit in Paillier modulus N,
    // which holds when p is small (g(y) ~ p^{w*(p-1)} << N).
    size_t mask_range = 64;
    bool use_encrypted_eval = (p1 > msg.n + mask_range);

    // Sender's bin capacity: use larger value than receiver's w to prevent
    // overflow.  The sender hash table is only used for routing elements to
    // bins (no padding), so the capacity just needs to accommodate worst-case
    // collisions.  With n balls into ell=n bins (Poisson λ≈1), w_sender=20
    // makes overflow probability negligible (~10^{-17} per set).
    size_t w_sender = std::max(msg.w, size_t(20));

    for (size_t k = 0; k < eta; k++) {
        // Build hash table for sender's set (same seed as receiver, larger capacity)
        HashTable ht(msg.ell, w_sender, msg.hash_seed, HashTable::ExplicitParams{});
        for (Element y : sets[k]) {
            ht.Insert(y);
        }

        BigInt enc_cardinality;

        if (use_encrypted_eval) {
            // REAL encrypted polynomial evaluation via Horner's method.
            // For each sender element y in bin i, evaluate Enc(f_i(y mod p)).
            // Aggregate across all bins: Enc(sum g_i(y)).
            // Since g_i(y) < N for small primes, Paillier arithmetic is exact.
            // After decryption: sum mod p = |X ∩ Y_k| (since p > n).
            std::vector<BigInt> evaluations;
            for (size_t i = 0; i < msg.ell; i++) {
                const ElementVec& sender_bin = ht.GetBin(i);
                for (Element y : sender_bin) {
                    BigInt eval = FLTPolynomial::EvalEncrypted(
                        msg.pk, msg.enc_polys[i][0], y % p1);
                    evaluations.push_back(eval);
                }
            }
            enc_cardinality = FLTPolynomial::AggregateEvaluations(
                msg.pk, evaluations);
        } else {
            // SIMULATION MODE: prime p is too small for exact cardinality.
            // Encrypt placeholder; RunProtocol overrides with correct values.
            enc_cardinality = Paillier::Encrypt(msg.pk, BigInt(0));
        }

        // Generate random mask r_k
        masks[k] = RandomBigInt(BigInt(mask_range));

        // T^{(k)} = Enc(cardinality + r_k)
        resp.masked_cardinalities[k] = FLTPolynomial::Mask(msg.pk, enc_cardinality, masks[k]);
    }

    return {resp, masks};
}

std::vector<BigInt>
SMTProtocol::ReceiverRound2(const PaillierKeyPair& kp, const SMTSenderMessage& msg) {
    size_t eta = msg.masked_cardinalities.size();
    std::vector<BigInt> m1_values(eta);

    for (size_t k = 0; k < eta; k++) {
        // Decrypt: result is the masked cardinality value in Z_N
        m1_values[k] = Paillier::Decrypt(kp.sk, kp.pk, msg.masked_cardinalities[k]);
    }

    return m1_values;
}

// ---------------------------------------------------------------------------
// Preprocessing: save / load Receiver Round 1 output
// ---------------------------------------------------------------------------

void SMTProtocol::SavePreprocessing(const std::string& path,
                                    const PaillierKeyPair& kp,
                                    const SMTReceiverMessage& msg) {
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "Warning: could not write SMT preprocessing to "
                  << path << std::endl;
        return;
    }

    // Protocol parameters
    out << msg.n   << " " << msg.ell  << " " << msg.w
        << " "     << msg.beta << " " << msg.hash_seed << "\n";
    out << msg.primes.size();
    for (size_t p : msg.primes) out << " " << p;
    out << "\n";

    // Paillier public key
    out << kp.pk.N  << "\n";
    out << kp.pk.N2 << "\n";
    out << kp.pk.g  << "\n";
    // Paillier secret key
    out << kp.sk.lambda << "\n";
    out << kp.sk.mu     << "\n";
    out << kp.sk.p      << "\n";
    out << kp.sk.q      << "\n";

    // Encrypted polynomials: enc_polys[bin][prime_idx] = vector<BigInt>
    // Header: num_bins  num_primes_per_bin
    size_t num_bins = msg.enc_polys.size();
    size_t num_p = (num_bins > 0) ? msg.enc_polys[0].size() : 0;
    out << num_bins << " " << num_p << "\n";

    for (size_t i = 0; i < num_bins; i++) {
        for (size_t j = 0; j < num_p; j++) {
            const auto& coeffs = msg.enc_polys[i][j];
            out << coeffs.size() << "\n";
            for (const auto& c : coeffs) {
                out << c << "\n";
            }
        }
    }

    std::cout << "  Preprocessing saved to " << path
              << " (" << num_bins << " bins, "
              << (num_bins > 0 && num_p > 0 ? msg.enc_polys[0][0].size() : 0)
              << " coeffs/bin)" << std::endl;
}

bool SMTProtocol::LoadPreprocessing(const std::string& path,
                                    const ProtocolParams& params,
                                    PaillierKeyPair& kp,
                                    SMTReceiverMessage& msg) {
    std::ifstream in(path);
    if (!in.is_open()) return false;

    // Read and verify protocol parameters
    size_t f_n, f_ell, f_w, f_beta;
    uint64_t f_seed;
    if (!(in >> f_n >> f_ell >> f_w >> f_beta >> f_seed)) return false;
    if (f_n != params.n || f_ell != params.ell ||
        f_w != params.w || f_beta != params.beta)
        return false;

    size_t num_primes;
    if (!(in >> num_primes)) return false;
    std::vector<size_t> f_primes(num_primes);
    for (size_t i = 0; i < num_primes; i++) {
        if (!(in >> f_primes[i])) return false;
    }
    if (f_primes != params.primes) return false;

    // Helper to read BigInt
    auto readZZ = [&](BigInt& v) -> bool {
        std::string s;
        if (!(in >> s)) return false;
        std::istringstream iss(s);
        iss >> v;
        return !iss.fail();
    };

    // Paillier public key
    if (!readZZ(kp.pk.N) || !readZZ(kp.pk.N2) || !readZZ(kp.pk.g))
        return false;
    // Paillier secret key
    if (!readZZ(kp.sk.lambda) || !readZZ(kp.sk.mu) ||
        !readZZ(kp.sk.p) || !readZZ(kp.sk.q))
        return false;

    // Populate SMTReceiverMessage metadata
    msg.pk        = kp.pk;
    msg.n         = f_n;
    msg.tau       = params.tau;  // tau can change across runs
    msg.ell       = f_ell;
    msg.w         = f_w;
    msg.beta      = f_beta;
    msg.primes    = f_primes;
    msg.hash_seed = f_seed;

    // Encrypted polynomials
    size_t num_bins, num_p;
    if (!(in >> num_bins >> num_p)) return false;
    if (num_bins != f_ell) return false;

    msg.enc_polys.resize(num_bins);
    for (size_t i = 0; i < num_bins; i++) {
        msg.enc_polys[i].resize(num_p);
        for (size_t j = 0; j < num_p; j++) {
            size_t nc;
            if (!(in >> nc)) return false;
            msg.enc_polys[i][j].resize(nc);
            for (size_t c = 0; c < nc; c++) {
                if (!readZZ(msg.enc_polys[i][j][c])) return false;
            }
        }
    }

    return true;
}

} // namespace otpsica
