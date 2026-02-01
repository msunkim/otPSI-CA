#ifndef OTPSICA_CEK_H
#define OTPSICA_CEK_H

#include "common.h"

namespace otpsica {

/**
 * CEK (Carlton-Essex-Kapulkin) encryption scheme for integer comparison.
 *
 * Based on: Carlton, R., Essex, A., Kapulkin, K. "Threshold Properties of
 * Prime Power Subgroups with Application to Secure Integer Comparisons."
 * CT-RSA 2018.
 *
 * Key insight: Uses subgroups of Z_N^* with special structure:
 *   - G_1 = <g_1> has order rho^d (small subgroup for comparison)
 *   - G_2 = <g_2> has order p_s * q_s (large subgroup for hiding)
 *
 * Encryption: c = g_1^{rho^m} * g_2^r
 * Decryption: Compute discrete log in small subgroup G_1
 */

struct CEKPublicKey {
    BigInt N;       // RSA modulus N = p*q
    BigInt rho;     // Small prime base
    size_t d;       // Exponent such that rho^d bounds message space
    BigInt g1;      // Generator of order rho^d
    BigInt g2;      // Generator of order p_s * q_s
    size_t kappa;   // Security parameter for DLog hardness
};

struct CEKSecretKey {
    BigInt x;       // Secret key: x = p_s * q_s * (1/(p_s*q_s) mod rho^d)
    BigInt p_s;     // Prime factor component
    BigInt q_s;     // Prime factor component
};

struct CEKKeyPair {
    CEKPublicKey pk;
    CEKSecretKey sk;
};

class CEK {
public:
    /**
     * Generate a CEK key pair.
     * @param kappa_bits Security parameter for subgroup (typically 256)
     * @param N_bits Size of RSA modulus (typically 2048)
     * @param rho Small prime (typically 2)
     * @param d Exponent (message space is [0, rho^d - 1])
     */
    static CEKKeyPair KeyGen(size_t kappa_bits = 256, size_t N_bits = 2048,
                             size_t rho = 2, size_t d = 17);

    /**
     * Encrypt a message m in [0, d-1].
     * @param pk Public key
     * @param m Message (must be < d)
     * @return Ciphertext in Z_N^*
     */
    static BigInt Encrypt(const CEKPublicKey& pk, size_t m);

    /**
     * Encrypt with specified randomness.
     */
    static BigInt Encrypt(const CEKPublicKey& pk, size_t m, const BigInt& r);

    /**
     * Decrypt a ciphertext.
     * @param sk Secret key
     * @param pk Public key
     * @param c Ciphertext
     * @return Message in [0, d-1], or -1 if decryption fails
     */
    static long Decrypt(const CEKSecretKey& sk, const CEKPublicKey& pk, const BigInt& c);

    /**
     * Decrypt to get raw discrete log (without rho-adic conversion).
     * Used when the ciphertext has been blinded with arbitrary exponent.
     * @return Raw discrete log in [0, rho^d - 1], or -1 if decryption fails
     */
    static long DecryptRaw(const CEKSecretKey& sk, const CEKPublicKey& pk, const BigInt& c);

    /**
     * Homomorphic addition in the exponent (limited).
     * Note: Enc(m1) * Enc(m2) = Enc(?) where ? depends on rho^{m1} + rho^{m2}
     */
    static BigInt Multiply(const CEKPublicKey& pk, const BigInt& c1, const BigInt& c2);

    /**
     * Raise ciphertext to power (for comparison protocol).
     * c^{rho^k} shifts the message: Enc(m)^{rho^k} relates to m + k
     */
    static BigInt Power(const CEKPublicKey& pk, const BigInt& c, const BigInt& exp);

    /**
     * Blind a ciphertext with random elements from both subgroups.
     * Result: c * g1^sigma * g2^r where sigma in [1, rho^d-1], r random
     */
    static BigInt Blind(const CEKPublicKey& pk, const BigInt& c,
                        const BigInt& sigma, const BigInt& r);

    /**
     * Get rho^d (the order of g1).
     */
    static BigInt GetRhoPowD(const CEKPublicKey& pk);

private:
    // Compute discrete log base g1 in the small subgroup
    static long DiscreteLogSmall(const CEKPublicKey& pk, const BigInt& h);

    // Compute log_rho(x) for small x
    static long LogRho(const BigInt& x, size_t rho, size_t d);
};

} // namespace otpsica

#endif // OTPSICA_CEK_H
