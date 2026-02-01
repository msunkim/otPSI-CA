#ifndef OTPSICA_FLT_POLY_H
#define OTPSICA_FLT_POLY_H

#include "common.h"
#include "crypto/paillier.h"

namespace otpsica {

/**
 * FLT (Fermat's Little Theorem) Polynomial for set membership testing.
 *
 * For a set V = {x_1, ..., x_w} and prime p, the polynomial:
 *   f(X) = 1 - prod_{i=1}^{w} (X - x_i)^{p-1}
 *
 * satisfies (by Fermat's Little Theorem):
 *   f(y) = 1 if y in V
 *   f(y) = 0 if y not in V (and p does not divide y - x_i for all i)
 *
 * Degree: deg(f) = w(p-1)
 *
 * For encrypted polynomial evaluation:
 *   F = Enc(f) = (Enc(f_{w(p-1)}), ..., Enc(f_1), Enc(f_0))
 *   F(y) = prod_{i=0}^{w(p-1)} F[i]^{y^i} = Enc(f(y))
 */

class FLTPolynomial {
public:
    /**
     * Build the FLT polynomial for a set of elements.
     * @param elements Set elements {x_1, ..., x_w}
     * @param p Prime for FLT (p-1 is the exponent)
     * @return Polynomial coefficients [f_{d}, f_{d-1}, ..., f_1, f_0]
     *         where d = w(p-1)
     */
    static std::vector<BigInt> Build(const ElementVec& elements, size_t p);

    /**
     * Encrypt a polynomial under Paillier.
     * @param pk Paillier public key
     * @param coeffs Polynomial coefficients [f_d, ..., f_0]
     * @return Encrypted coefficients [Enc(f_d), ..., Enc(f_0)]
     */
    static std::vector<BigInt> Encrypt(const PaillierPublicKey& pk,
                                       const std::vector<BigInt>& coeffs);

    /**
     * Evaluate an encrypted polynomial at a point.
     * @param pk Paillier public key
     * @param enc_coeffs Encrypted coefficients [Enc(f_d), ..., Enc(f_0)]
     * @param y Evaluation point
     * @return Enc(f(y)) = prod_{i=0}^{d} Enc(f_i)^{y^i}
     */
    static BigInt EvalEncrypted(const PaillierPublicKey& pk,
                                const std::vector<BigInt>& enc_coeffs,
                                Element y);

    /**
     * Aggregate encrypted evaluations (for cardinality counting).
     * @param pk Paillier public key
     * @param evaluations Vector of Enc(f(y_i)) values
     * @return Enc(sum of f(y_i)) = prod of Enc(f(y_i))
     */
    static BigInt AggregateEvaluations(const PaillierPublicKey& pk,
                                       const std::vector<BigInt>& evaluations);

    /**
     * Mask an encrypted value with randomness.
     * @param pk Paillier public key
     * @param enc_val Encrypted value Enc(v)
     * @param r Random mask
     * @return Enc(v + r)
     */
    static BigInt Mask(const PaillierPublicKey& pk,
                       const BigInt& enc_val,
                       const BigInt& r);

private:
    // Compute (X - a)^{p-1} mod p using fast exponentiation
    static NTL::ZZ_pX PowerPoly(const NTL::ZZ_p& a, size_t exp);
};

/**
 * CRT (Chinese Remainder Theorem) decomposition for degree reduction.
 *
 * For beta-bit elements, we need primes with product > 2^beta.
 * Using nu small primes {p_1, ..., p_nu} instead of one large prime
 * reduces polynomial degree from w(p-1) to max_j{w(p_j-1)}.
 */

class CRTDecomposition {
public:
    /**
     * Select CRT primes for given element bit-size.
     * @param beta Element bit-size (e.g., 16)
     * @param nu Number of primes (e.g., 2)
     * @return Vector of primes with product > 2^beta
     */
    static std::vector<size_t> SelectPrimes(size_t beta, size_t nu = 2);

    /**
     * Check if the prime selection is valid.
     * @param primes Selected primes
     * @param beta Element bit-size
     * @return true if product of primes > 2^beta
     */
    static bool Validate(const std::vector<size_t>& primes, size_t beta);

    /**
     * Get maximum polynomial degree for given primes and bin size.
     */
    static size_t MaxDegree(const std::vector<size_t>& primes, size_t w);
};

} // namespace otpsica

#endif // OTPSICA_FLT_POLY_H
