#ifndef OTPSICA_PAILLIER_H
#define OTPSICA_PAILLIER_H

#include "common.h"

namespace otpsica {

/**
 * Paillier cryptosystem for additively homomorphic encryption.
 *
 * Based on: Paillier, P. "Public-Key Cryptosystems Based on Composite
 * Degree Residuosity Classes." EUROCRYPT 1999.
 *
 * Message space: Z_N
 * Ciphertext space: Z_{N^2}
 *
 * Homomorphic properties:
 *   Enc(m1) * Enc(m2) = Enc(m1 + m2)
 *   Enc(m)^k = Enc(k * m)
 */

struct PaillierPublicKey {
    BigInt N;       // RSA modulus N = p*q
    BigInt N2;      // N^2 (ciphertext space modulus)
    BigInt g;       // Generator (typically 1 + N)
};

struct PaillierSecretKey {
    BigInt lambda;  // lcm(p-1, q-1)
    BigInt mu;      // L(g^lambda mod N^2)^{-1} mod N
    BigInt p;       // Prime factor
    BigInt q;       // Prime factor
};

struct PaillierKeyPair {
    PaillierPublicKey pk;
    PaillierSecretKey sk;
};

class Paillier {
public:
    /**
     * Generate a Paillier key pair.
     * @param bits Number of bits for N (should be >= 2048 for security)
     */
    static PaillierKeyPair KeyGen(size_t bits = 2048);

    /**
     * Encrypt a plaintext message.
     * @param pk Public key
     * @param m Plaintext in Z_N
     * @return Ciphertext in Z_{N^2}
     */
    static BigInt Encrypt(const PaillierPublicKey& pk, const BigInt& m);

    /**
     * Encrypt with specified randomness (for testing/verification).
     */
    static BigInt Encrypt(const PaillierPublicKey& pk, const BigInt& m, const BigInt& r);

    /**
     * Decrypt a ciphertext.
     * @param sk Secret key
     * @param pk Public key
     * @param c Ciphertext in Z_{N^2}
     * @return Plaintext in Z_N
     */
    static BigInt Decrypt(const PaillierSecretKey& sk, const PaillierPublicKey& pk, const BigInt& c);

    /**
     * Homomorphic addition: Enc(m1 + m2) = Enc(m1) * Enc(m2)
     */
    static BigInt Add(const PaillierPublicKey& pk, const BigInt& c1, const BigInt& c2);

    /**
     * Homomorphic scalar multiplication: Enc(k * m) = Enc(m)^k
     */
    static BigInt ScalarMul(const PaillierPublicKey& pk, const BigInt& c, const BigInt& k);

    /**
     * Homomorphic negation: Enc(-m) = Enc(m)^{-1}
     */
    static BigInt Negate(const PaillierPublicKey& pk, const BigInt& c);

    /**
     * Re-randomize a ciphertext (produces fresh randomness).
     */
    static BigInt Rerandomize(const PaillierPublicKey& pk, const BigInt& c);

private:
    // L function: L(x) = (x - 1) / N
    static BigInt L(const BigInt& x, const BigInt& N);
};

} // namespace otpsica

#endif // OTPSICA_PAILLIER_H
