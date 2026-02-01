#include "protocol/gs19.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <chrono>
#include <iostream>

namespace otpsica {

using namespace NTL;

// ============================================================
// Parameter initialization
// ============================================================

BigInt GS19PICT::FindPrime(const BigInt& bound) {
    BigInt p = bound + 1;
    // Make odd if even
    if (IsZero(p % BigInt(2))) p += 1;
    while (!ProbPrime(p)) {
        p += 2;
    }
    return p;
}

GS19Params GS19PICT::InitParams(size_t n, size_t tau, size_t kappa) {
    GS19Params params;
    params.n = n;
    params.tau = tau;
    params.kappa = kappa;
    params.t = n - tau;
    params.m = 2 * params.t + 1;

    // q > (4t^2 + 2t) * 2^kappa
    BigInt t_zz = conv<ZZ>(params.t);
    BigInt bound = (4 * t_zz * t_zz + 2 * t_zz) * power2_ZZ(kappa);
    params.q = FindPrime(bound);

    return params;
}

// ============================================================
// Power-sum polynomial encoding
// ============================================================

std::vector<BigInt> GS19PICT::ComputePowerSums(
    const ElementVec& S, const BigInt& u,
    size_t max_power, const BigInt& q) {

    // Set field modulus
    ZZ_p::init(q);

    ZZ_p u_p = conv<ZZ_p>(u);

    // result[i] = p_S(u^i) = Σ_{a∈S} u^{i*a} = Σ_{a∈S} (u^a)^i
    std::vector<BigInt> result(max_power + 1, BigInt(0));

    // For each element a in S, compute base_a = u^a mod q,
    // then accumulate base_a^i into result[i]
    for (Element a : S) {
        ZZ_p base_a = power(u_p, a);  // u^a mod q
        ZZ_p running = conv<ZZ_p>(1); // (u^a)^0 = 1

        for (size_t i = 0; i <= max_power; i++) {
            result[i] += rep(running);
            running *= base_a;
        }
    }

    // Reduce mod q
    for (size_t i = 0; i <= max_power; i++) {
        result[i] %= q;
    }

    return result;
}

// ============================================================
// Plaintext verification (no encryption)
// ============================================================

bool GS19PICT::PlaintextTest(const ElementVec& X, const ElementVec& Y,
                              const GS19Params& params) {
    ZZ_p::init(params.q);

    // Pick random evaluation point
    ZZ_p u = random_ZZ_p();
    while (IsZero(u)) u = random_ZZ_p();

    BigInt u_zz = rep(u);
    size_t diag_len = 4 * params.t;

    // Compute power sums for both sets
    auto c_A = ComputePowerSums(X, u_zz, diag_len, params.q);
    auto c_B = ComputePowerSums(Y, u_zz, diag_len, params.q);

    // Build H_C = H_A - H_B as mat_ZZ_p
    ZZ_p::init(params.q);
    mat_ZZ_p H_C;
    H_C.SetDims(params.m, params.m);

    for (size_t i = 0; i < params.m; i++) {
        for (size_t j = 0; j < params.m; j++) {
            ZZ_p ca = conv<ZZ_p>(c_A[i + j]);
            ZZ_p cb = conv<ZZ_p>(c_B[i + j]);
            H_C[i][j] = ca - cb;
        }
    }

    // Singular iff determinant = 0
    ZZ_p det = determinant(H_C);
    return IsZero(det);  // singular => |sym diff| <= 2t => |intersection| >= tau
}

// ============================================================
// Alice offline phase
// ============================================================

GS19PICT::AliceState GS19PICT::AliceOffline(const ElementVec& X,
                                              const GS19Params& params) {
    AliceState state;
    state.params = params;

    // Generate Paillier keypair (1024-bit for experiments)
    state.kp = Paillier::KeyGen(1024);

    // Pick random evaluation point u in F_q
    ZZ_p::init(params.q);
    ZZ_p u = random_ZZ_p();
    while (IsZero(u)) u = random_ZZ_p();
    state.u_val = rep(u);

    // Compute power-sum polynomial evaluations
    // c[i] = p_A(u^i) for i = 0, ..., 4t
    size_t diag_len = 4 * params.t;
    auto c = ComputePowerSums(X, state.u_val, diag_len, params.q);

    // Encrypt the Hankel diagonal entries
    // Only 4t+1 values needed (Hankel reuses entries: H[i][j] = c[i+j])
    state.enc_diag.resize(diag_len + 1);
    for (size_t i = 0; i <= diag_len; i++) {
        state.enc_diag[i] = Paillier::Encrypt(state.kp.pk, c[i]);
    }

    return state;
}

// ============================================================
// Encrypted Hankel matrix-vector multiplication
// ============================================================

std::vector<BigInt> GS19PICT::EncHankelVecMul(
    const PaillierPublicKey& pk,
    const std::vector<BigInt>& enc_diag,
    const std::vector<BigInt>& w_plain,
    size_t m) {

    std::vector<BigInt> result(m);

    for (size_t r = 0; r < m; r++) {
        // result[r] = Σ_{c=0}^{m-1} Enc(H[r][c]) * w[c]
        //           = Σ_{c} Enc(enc_diag[r+c])^{w[c]}
        BigInt acc = Paillier::Encrypt(pk, BigInt(0));  // Enc(0) as identity

        for (size_t c = 0; c < m; c++) {
            // Enc(H[r][c])^{w[c]} = ScalarMul(enc_diag[r+c], w[c])
            BigInt term = Paillier::ScalarMul(pk, enc_diag[r + c], w_plain[c]);
            acc = Paillier::Add(pk, acc, term);
        }

        result[r] = acc;
    }

    return result;
}

// ============================================================
// FINV: Run for one sender set (full protocol)
// ============================================================

bool GS19PICT::RunForSender(const AliceState& state, const ElementVec& Y) {
    const auto& params = state.params;
    const auto& pk = state.kp.pk;
    size_t m = params.m;
    size_t diag_len = 4 * params.t;

    // --- Bob: compute power sums for Y and homomorphic subtraction ---
    auto c_B = ComputePowerSums(Y, state.u_val, diag_len, params.q);

    // enc_diag_C[i] = Enc(c_A[i] - c_B[i]) = Enc(c_A[i]) * Enc(-c_B[i])
    std::vector<BigInt> enc_diag_C(diag_len + 1);
    for (size_t i = 0; i <= diag_len; i++) {
        // -c_B[i] mod N (in Paillier plaintext space)
        BigInt neg_cB = (-c_B[i]) % pk.N;
        if (neg_cB < 0) neg_cB += pk.N;
        BigInt enc_neg = Paillier::Encrypt(pk, neg_cB);
        enc_diag_C[i] = Paillier::Add(pk, state.enc_diag[i], enc_neg);
    }

    // --- FINV protocol (sequential iterations) ---
    ZZ_p::init(params.q);

    // Bob picks random vectors u_vec, v_vec in F_q^m
    vec_ZZ_p u_vec, v_vec;
    u_vec.SetLength(m);
    v_vec.SetLength(m);
    for (size_t i = 0; i < m; i++) {
        u_vec[i] = random_ZZ_p();
        v_vec[i] = random_ZZ_p();
    }

    // w starts as v_vec (plaintext, as BigInt for Paillier interface)
    std::vector<BigInt> w(m);
    for (size_t i = 0; i < m; i++) {
        w[i] = rep(v_vec[i]);
    }

    // Compute sequence a_0, ..., a_{2m-1}
    vec_ZZ_p sequence;
    sequence.SetLength(2 * m);
    BigInt half_N = pk.N / 2;

    for (size_t iter = 0; iter < 2 * m; iter++) {
        // a[iter] = u_vec · w (plaintext dot product in F_q)
        // Note: ZZ_p modulus is already set to params.q above; do NOT re-init
        // as that would invalidate u_vec, v_vec, and sequence.
        ZZ_p dot = conv<ZZ_p>(0);
        for (size_t j = 0; j < m; j++) {
            dot += u_vec[j] * conv<ZZ_p>(w[j]);
        }
        sequence[iter] = dot;

        // Last iteration: don't need to compute next w
        if (iter == 2 * m - 1) break;

        // Encrypted mat-vec: enc_w_new = Enc(H_C · w)
        auto enc_w_new = EncHankelVecMul(pk, enc_diag_C, w, m);

        // Alice decrypts and reduces mod q
        // Paillier plaintext space is Z_N. Negative values are represented
        // as N + val. Since |true value| << N/2, we center around 0.
        for (size_t j = 0; j < m; j++) {
            BigInt plain = Paillier::Decrypt(state.kp.sk, pk, enc_w_new[j]);
            if (plain > half_N) {
                plain -= pk.N;  // Negative value in centered representation
            }
            w[j] = plain % params.q;
            if (w[j] < 0) w[j] += params.q;
        }
    }

    // --- Berlekamp-Massey: compute minimal polynomial ---
    // ZZ_p modulus is still params.q from above
    ZZ_pX min_poly;
    MinPolySeq(min_poly, sequence, m);

    // Singular iff constant coefficient = 0
    bool singular = IsZero(coeff(min_poly, 0));
    return singular;  // singular => |intersection| >= tau
}

// ============================================================
// FINV: Partial run for benchmarking (large parameters)
// ============================================================

bool GS19PICT::RunForSenderPartial(const AliceState& state, const ElementVec& Y,
                                    size_t max_iters, double& per_iter_ms) {
    const auto& params = state.params;
    const auto& pk = state.kp.pk;
    size_t m = params.m;
    size_t diag_len = 4 * params.t;

    // Bob: compute power sums and homomorphic subtraction
    auto c_B = ComputePowerSums(Y, state.u_val, diag_len, params.q);

    std::vector<BigInt> enc_diag_C(diag_len + 1);
    for (size_t i = 0; i <= diag_len; i++) {
        BigInt neg_cB = (-c_B[i]) % pk.N;
        if (neg_cB < 0) neg_cB += pk.N;
        BigInt enc_neg = Paillier::Encrypt(pk, neg_cB);
        enc_diag_C[i] = Paillier::Add(pk, state.enc_diag[i], enc_neg);
    }

    // Initialize FINV
    ZZ_p::init(params.q);
    std::vector<BigInt> w(m);
    for (size_t i = 0; i < m; i++) {
        ZZ_p r = random_ZZ_p();
        w[i] = rep(r);
    }

    // Run limited iterations and measure time
    size_t iters = std::min(max_iters, 2 * m);
    auto t0 = std::chrono::high_resolution_clock::now();

    BigInt half_N = pk.N / 2;
    for (size_t iter = 0; iter < iters; iter++) {
        auto enc_w_new = EncHankelVecMul(pk, enc_diag_C, w, m);

        for (size_t j = 0; j < m; j++) {
            BigInt plain = Paillier::Decrypt(state.kp.sk, pk, enc_w_new[j]);
            if (plain > half_N) plain -= pk.N;
            w[j] = plain % params.q;
            if (w[j] < 0) w[j] += params.q;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    per_iter_ms = total_ms / iters;

    return true;  // Not a real result; just benchmarking
}

// ============================================================
// Convenience wrappers
// ============================================================

bool GS19PICT::RunProtocol(const ElementVec& X, const ElementVec& Y,
                            const GS19Params& params) {
    auto state = AliceOffline(X, params);
    return RunForSender(state, Y);
}

std::vector<bool> GS19PICT::RunProtocolMultiSender(
    const ElementVec& X, const std::vector<ElementVec>& Y_sets,
    const GS19Params& params) {

    auto state = AliceOffline(X, params);

    std::vector<bool> results(Y_sets.size());
    for (size_t k = 0; k < Y_sets.size(); k++) {
        results[k] = RunForSender(state, Y_sets[k]);
    }

    return results;
}

} // namespace otpsica
