#include "crypto/cek.h"
#include <NTL/ZZ.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

namespace otpsica {

using namespace NTL;

// ---------------------------------------------------------------------------
// File-cache helpers for CEK key generation.
//
// Cache file format (all values decimal, one per line):
//   line 0 : kappa_bits  N_bits  rho  d          (parameters for matching)
//   line 1 : p_s
//   line 2 : q_s
//   line 3 : p_t
//   line 4 : q_t
//   line 5 : N
//   line 6 : g1
//   line 7 : g2
//   line 8 : x
// ---------------------------------------------------------------------------

static bool TryLoadCEKCache(const std::string& path,
                            size_t kappa_bits, size_t N_bits,
                            size_t rho, size_t d,
                            CEKKeyPair& kp) {
    std::ifstream in(path);
    if (!in.is_open()) return false;

    // Verify parameters on the first line
    size_t f_kappa, f_Nbits, f_rho, f_d;
    if (!(in >> f_kappa >> f_Nbits >> f_rho >> f_d)) return false;
    if (f_kappa != kappa_bits || f_Nbits != N_bits ||
        f_rho != rho || f_d != d)
        return false;

    // Helper: read one NTL BigInt from the next non-empty token
    auto readZZ = [&](BigInt& v) -> bool {
        std::string s;
        if (!(in >> s)) return false;
        std::istringstream iss(s);
        iss >> v;
        return !iss.fail();
    };

    BigInt p_s, q_s, p_t, q_t, N, g1, g2, x;
    if (!readZZ(p_s) || !readZZ(q_s) ||
        !readZZ(p_t) || !readZZ(q_t) ||
        !readZZ(N)   || !readZZ(g1)  ||
        !readZZ(g2)  || !readZZ(x))
        return false;

    kp.pk.rho   = conv<BigInt>(rho);
    kp.pk.d     = d;
    kp.pk.kappa = kappa_bits;
    kp.pk.N     = N;
    kp.pk.g1    = g1;
    kp.pk.g2    = g2;
    kp.sk.x     = x;
    kp.sk.p_s   = p_s;
    kp.sk.q_s   = q_s;
    return true;
}

static void SaveCEKCache(const std::string& path,
                         size_t kappa_bits, size_t N_bits,
                         size_t rho, size_t d,
                         const BigInt& p_s, const BigInt& q_s,
                         const BigInt& p_t, const BigInt& q_t,
                         const CEKKeyPair& kp) {
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "Warning: could not write CEK cache to "
                  << path << std::endl;
        return;
    }
    out << kappa_bits << " " << N_bits << " " << rho << " " << d << "\n";
    out << p_s       << "\n";
    out << q_s       << "\n";
    out << p_t       << "\n";
    out << q_t       << "\n";
    out << kp.pk.N   << "\n";
    out << kp.pk.g1  << "\n";
    out << kp.pk.g2  << "\n";
    out << kp.sk.x   << "\n";
}

// ---------------------------------------------------------------------------

CEKKeyPair CEK::KeyGen(size_t kappa_bits, size_t N_bits, size_t rho, size_t d) {
    // --- Try loading from cache file first ---
    const std::string cache_path = "cek_kg.txt";
    CEKKeyPair kp;
    if (TryLoadCEKCache(cache_path, kappa_bits, N_bits, rho, d, kp)) {
        std::cout << "  [CEK] Loaded cached key (d=" << d
                  << ") from " << cache_path << std::endl;
        return kp;
    }

    std::cout << "  [CEK] Generating key (d=" << d
              << ") — this may be slow..." << std::endl;

    kp.pk.rho = conv<BigInt>(rho);
    kp.pk.d = d;
    kp.pk.kappa = kappa_bits;

    // Compute rho^d
    BigInt rho_d = power(kp.pk.rho, d);

    // Generate special primes p, q of the form:
    // p = 2 * rho^d * p_s * p_t + 1
    // q = 2 * rho^d * q_s * q_t + 1
    // where p_s, q_s are kappa-bit primes and p_t, q_t are chosen
    // to make p, q the appropriate size (N_bits / 2 each).

    size_t remaining_bits = N_bits / 2 - NumBits(rho_d) - 1;
    size_t ps_bits = kappa_bits;
    size_t pt_bits = remaining_bits - ps_bits;

    BigInt p_s, q_s, p_t, q_t, p, q;

    // Generate p
    size_t attempts = 0;
    do {
        GenPrime(p_s, ps_bits);
        GenPrime(p_t, pt_bits);
        p = 2 * rho_d * p_s * p_t + 1;
        ++attempts;
        if (attempts % 100 == 0)
            std::cout << "  [CEK]   p: " << attempts
                      << " attempts..." << std::endl;
    } while (!ProbPrime(p));
    std::cout << "  [CEK]   Found p after " << attempts
              << " attempts" << std::endl;

    // Generate q (different from p)
    attempts = 0;
    do {
        GenPrime(q_s, ps_bits);
        GenPrime(q_t, pt_bits);
        q = 2 * rho_d * q_s * q_t + 1;
        ++attempts;
        if (attempts % 100 == 0)
            std::cout << "  [CEK]   q: " << attempts
                      << " attempts..." << std::endl;
    } while (!ProbPrime(q) || p == q);
    std::cout << "  [CEK]   Found q after " << attempts
              << " attempts" << std::endl;

    kp.pk.N = p * q;

    // Compute secret key: x = p_s * q_s * ((p_s * q_s)^{-1} mod rho^d)
    BigInt ps_qs = p_s * q_s;
    BigInt inv_ps_qs = InvMod(ps_qs % rho_d, rho_d);
    kp.sk.x = ps_qs * inv_ps_qs;
    kp.sk.p_s = p_s;
    kp.sk.q_s = q_s;

    // Find generator g1 of order rho^d in Z_N^*
    // g1 = h^{2 * p_s * p_t * q_s * q_t} mod N for random h
    BigInt order_cofactor = 2 * p_s * p_t * q_s * q_t;
    BigInt h;
    do {
        h = RandomBigIntNonZero(kp.pk.N);
        kp.pk.g1 = PowerMod(h, order_cofactor, kp.pk.N);
    } while (IsOne(kp.pk.g1) ||
             IsOne(PowerMod(kp.pk.g1, rho_d / kp.pk.rho, kp.pk.N)));

    // Find generator g2 of order p_s * q_s in Z_N^*
    // g2 = h^{2 * rho^d * p_t * q_t} mod N for random h
    BigInt order_cofactor2 = 2 * rho_d * p_t * q_t;
    do {
        h = RandomBigIntNonZero(kp.pk.N);
        kp.pk.g2 = PowerMod(h, order_cofactor2, kp.pk.N);
    } while (IsOne(kp.pk.g2));

    // Save to cache for future runs
    SaveCEKCache(cache_path, kappa_bits, N_bits, rho, d,
                 p_s, q_s, p_t, q_t, kp);
    std::cout << "  [CEK] Key saved to " << cache_path << std::endl;

    return kp;
}

BigInt CEK::Encrypt(const CEKPublicKey& pk, size_t m) {
    BigInt r = RandomBigInt(power2_ZZ(pk.kappa));
    return Encrypt(pk, m, r);
}

BigInt CEK::Encrypt(const CEKPublicKey& pk, size_t m, const BigInt& r) {
    // c = g1^{rho^m} * g2^r mod N
    BigInt rho_m = power(pk.rho, m);
    BigInt g1_part = PowerMod(pk.g1, rho_m, pk.N);
    BigInt g2_part = PowerMod(pk.g2, r, pk.N);
    return MulMod(g1_part, g2_part, pk.N);
}

long CEK::Decrypt(const CEKSecretKey& sk, const CEKPublicKey& pk, const BigInt& c) {
    // a1 = c^x mod N
    // This projects onto the subgroup <g1> and yields g1^{rho^m}
    BigInt a1 = PowerMod(c, sk.x, pk.N);

    if (IsOne(a1)) {
        return 0;
    }

    // Compute discrete log: find m such that a1 = g1^{rho^m}
    // First find a2 = log_{g1}(a1), then m = log_rho(a2)
    long log_val = DiscreteLogSmall(pk, a1);
    if (log_val < 0) {
        return -1;  // Decryption failed
    }

    return LogRho(conv<BigInt>(log_val), conv<long>(pk.rho), pk.d);
}

long CEK::DecryptRaw(const CEKSecretKey& sk, const CEKPublicKey& pk, const BigInt& c) {
    // Returns the raw discrete log without rho-adic conversion
    // Used when ciphertext has been blinded with arbitrary exponent
    BigInt a1 = PowerMod(c, sk.x, pk.N);

    if (IsOne(a1)) {
        return 0;
    }

    return DiscreteLogSmall(pk, a1);
}

BigInt CEK::Multiply(const CEKPublicKey& pk, const BigInt& c1, const BigInt& c2) {
    return MulMod(c1, c2, pk.N);
}

BigInt CEK::Power(const CEKPublicKey& pk, const BigInt& c, const BigInt& exp) {
    return PowerMod(c, exp, pk.N);
}

BigInt CEK::Blind(const CEKPublicKey& pk, const BigInt& c,
                  const BigInt& sigma, const BigInt& r) {
    // c * g1^sigma * g2^r mod N
    BigInt g1_sigma = PowerMod(pk.g1, sigma, pk.N);
    BigInt g2_r = PowerMod(pk.g2, r, pk.N);
    BigInt result = MulMod(c, g1_sigma, pk.N);
    return MulMod(result, g2_r, pk.N);
}

BigInt CEK::GetRhoPowD(const CEKPublicKey& pk) {
    return power(pk.rho, pk.d);
}

long CEK::DiscreteLogSmall(const CEKPublicKey& pk, const BigInt& h) {
    // Baby-step giant-step in the small subgroup of order rho^d
    // Since rho is small (typically 2) and d is small (typically <= 17),
    // we can use exhaustive search or Pollard's rho

    BigInt rho_d = GetRhoPowD(pk);
    BigInt current = IsOne(h) ? BigInt(1) : h;

    // For very small subgroups, exhaustive search is efficient
    // Check powers: g1^0, g1^1, g1^2, ..., g1^{rho^d - 1}
    BigInt target = h;
    BigInt power_val = BigInt(1);

    for (long i = 0; i < conv<long>(rho_d); i++) {
        if (IsOne(power_val) && i > 0) {
            break;  // Wrapped around
        }
        if (power_val == target) {
            return i;
        }
        power_val = MulMod(power_val, pk.g1, pk.N);
    }

    return -1;  // Not found
}

long CEK::LogRho(const BigInt& x, size_t rho, size_t d) {
    // Find m such that x = rho^m (for m in [0, d])
    BigInt rho_power = BigInt(1);
    for (size_t m = 0; m <= d; m++) {
        if (rho_power == x) {
            return static_cast<long>(m);
        }
        rho_power *= rho;
    }
    return -1;  // Not a power of rho
}

} // namespace otpsica
