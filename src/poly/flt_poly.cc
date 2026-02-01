#include "poly/flt_poly.h"
#include <NTL/ZZ_pX.h>

namespace otpsica {

using namespace NTL;

ZZ_pX FLTPolynomial::PowerPoly(const ZZ_p& a, size_t exp) {
    // Compute (X - a)^exp using fast polynomial exponentiation
    ZZ_pX base;
    SetCoeff(base, 1, 1);      // X
    SetCoeff(base, 0, -a);     // X - a

    ZZ_pX result;
    SetCoeff(result, 0, 1);    // 1

    while (exp > 0) {
        if (exp & 1) {
            result = result * base;
        }
        base = base * base;
        exp >>= 1;
    }

    return result;
}

std::vector<BigInt> FLTPolynomial::Build(const ElementVec& elements, size_t p) {
    // Set the modulus for polynomial arithmetic
    ZZ_p::init(conv<ZZ>(p));

    // Build f(X) = 1 - prod_{x in elements} (X - x)^{p-1}
    ZZ_pX product;
    SetCoeff(product, 0, 1);  // Start with 1

    for (Element x : elements) {
        ZZ_p x_mod = conv<ZZ_p>(x % p);
        ZZ_pX term = PowerPoly(x_mod, p - 1);
        product = product * term;
    }

    // f(X) = 1 - product
    ZZ_pX f;
    SetCoeff(f, 0, 1);
    f = f - product;

    // Extract coefficients
    size_t degree = deg(f);
    if (degree < 0) {
        // f is zero polynomial (shouldn't happen in normal use)
        return {BigInt(0)};
    }

    std::vector<BigInt> coeffs(degree + 1);
    for (size_t i = 0; i <= degree; i++) {
        coeffs[degree - i] = rep(coeff(f, i));  // Store high to low degree
    }

    return coeffs;
}

std::vector<BigInt> FLTPolynomial::Encrypt(const PaillierPublicKey& pk,
                                            const std::vector<BigInt>& coeffs) {
    std::vector<BigInt> enc_coeffs(coeffs.size());

    for (size_t i = 0; i < coeffs.size(); i++) {
        enc_coeffs[i] = Paillier::Encrypt(pk, coeffs[i]);
    }

    return enc_coeffs;
}

BigInt FLTPolynomial::EvalEncrypted(const PaillierPublicKey& pk,
                                     const std::vector<BigInt>& enc_coeffs,
                                     Element y) {
    // Evaluate encrypted polynomial at y using Horner's method
    // f(y) = f_d*y^d + ... + f_1*y + f_0
    //      = ((... (f_d * y + f_{d-1}) * y + ...) * y + f_0)
    // Horner keeps scalar multiplication bounded to y (not y^d)

    if (enc_coeffs.empty()) {
        return Paillier::Encrypt(pk, BigInt(0));
    }

    // coeffs are stored [f_d, f_{d-1}, ..., f_1, f_0]
    BigInt y_zz;
    conv(y_zz, static_cast<unsigned long>(y));

    // Horner's method: start with highest degree coefficient
    BigInt result = enc_coeffs[0];  // Enc(f_d)

    // Iterate from f_{d-1} down to f_0
    for (size_t i = 1; i < enc_coeffs.size(); i++) {
        // result = result * y + f_{d-i}
        // Enc(result * y) = Enc(result)^y
        result = Paillier::ScalarMul(pk, result, y_zz);
        // Enc(result * y + f_{d-i}) = Enc(result * y) * Enc(f_{d-i})
        result = Paillier::Add(pk, result, enc_coeffs[i]);
    }

    return result;
}

BigInt FLTPolynomial::AggregateEvaluations(const PaillierPublicKey& pk,
                                            const std::vector<BigInt>& evaluations) {
    // Aggregate: Enc(sum f(y_i)) = prod Enc(f(y_i))
    if (evaluations.empty()) {
        return Paillier::Encrypt(pk, BigInt(0));
    }

    BigInt result = evaluations[0];
    for (size_t i = 1; i < evaluations.size(); i++) {
        result = Paillier::Add(pk, result, evaluations[i]);
    }

    return result;
}

BigInt FLTPolynomial::Mask(const PaillierPublicKey& pk,
                            const BigInt& enc_val,
                            const BigInt& r) {
    // Enc(v + r) = Enc(v) * Enc(r)
    BigInt enc_r = Paillier::Encrypt(pk, r);
    return Paillier::Add(pk, enc_val, enc_r);
}

// CRT Decomposition

std::vector<size_t> CRTDecomposition::SelectPrimes(size_t beta, size_t nu) {
    // Find primes such that their product > 2^beta
    // Use smaller primes for smaller beta to reduce polynomial degree
    // nu is a hint for minimum primes, but we'll add more if needed

    // All candidate primes - smaller ones first for efficiency
    std::vector<size_t> candidates;

    if (beta <= 8) {
        // For small element ranges (up to 256), use small primes
        // Smaller primes = smaller polynomial degree = faster computation
        candidates = {11, 13, 17, 19, 23, 29, 31, 37, 41, 43};
    } else if (beta <= 16) {
        // Medium range (up to 65536)
        // Keep small primes to bound polynomial degree w*(p-1).
        // Product of enough primes easily exceeds 2^16 = 65536.
        candidates = {31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    } else {
        // Large range - original primes
        candidates = {
            257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
            313, 317, 331, 337, 347, 349, 353, 359, 367, 373
        };
    }

    std::vector<size_t> selected;
    BigInt product = BigInt(1);
    BigInt target = power2_ZZ(beta);

    // Add primes from candidates until product > target
    for (size_t p : candidates) {
        selected.push_back(p);
        product *= p;
        if (product > target) break;
    }

    // If still not enough (shouldn't happen with good candidates), extend
    if (product <= target) {
        size_t next_prime = candidates.back() + 2;
        while (product <= target) {
            while (!ProbPrime(conv<ZZ>(next_prime))) {
                next_prime++;
            }
            selected.push_back(next_prime);
            product *= next_prime;
            next_prime += 2;
        }
    }

    return selected;
}

bool CRTDecomposition::Validate(const std::vector<size_t>& primes, size_t beta) {
    BigInt product = BigInt(1);
    for (size_t p : primes) {
        product *= p;
    }
    return product > power2_ZZ(beta);
}

size_t CRTDecomposition::MaxDegree(const std::vector<size_t>& primes, size_t w) {
    size_t max_p = 0;
    for (size_t p : primes) {
        if (p > max_p) max_p = p;
    }
    return w * (max_p - 1);
}

} // namespace otpsica
