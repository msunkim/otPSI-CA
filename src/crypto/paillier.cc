#include "crypto/paillier.h"
#include <NTL/ZZ.h>

namespace otpsica {

using namespace NTL;

BigInt Paillier::L(const BigInt& x, const BigInt& N) {
    return (x - 1) / N;
}

PaillierKeyPair Paillier::KeyGen(size_t bits) {
    PaillierKeyPair kp;

    // Generate two large primes p, q of equal size
    size_t prime_bits = bits / 2;
    BigInt p, q;

    GenPrime(p, prime_bits);
    do {
        GenPrime(q, prime_bits);
    } while (p == q);

    // Compute N = p * q
    kp.pk.N = p * q;
    kp.pk.N2 = kp.pk.N * kp.pk.N;

    // g = 1 + N (simplified generator)
    kp.pk.g = 1 + kp.pk.N;

    // lambda = lcm(p-1, q-1)
    BigInt p1 = p - 1;
    BigInt q1 = q - 1;
    BigInt gcd_pq;
    GCD(gcd_pq, p1, q1);
    kp.sk.lambda = (p1 * q1) / gcd_pq;

    // mu = L(g^lambda mod N^2)^{-1} mod N
    BigInt g_lambda = PowerMod(kp.pk.g, kp.sk.lambda, kp.pk.N2);
    BigInt L_val = L(g_lambda, kp.pk.N);
    kp.sk.mu = InvMod(L_val, kp.pk.N);

    // Store factors for CRT optimization
    kp.sk.p = p;
    kp.sk.q = q;

    return kp;
}

BigInt Paillier::Encrypt(const PaillierPublicKey& pk, const BigInt& m) {
    // Random r in Z_N^*
    BigInt r = RandomBigIntNonZero(pk.N);
    return Encrypt(pk, m, r);
}

BigInt Paillier::Encrypt(const PaillierPublicKey& pk, const BigInt& m, const BigInt& r) {
    // c = (1 + mN) * r^N mod N^2
    // Using g = 1 + N: g^m = (1 + N)^m = 1 + mN mod N^2 (by binomial theorem)
    BigInt gm = (1 + (m % pk.N) * pk.N) % pk.N2;
    BigInt rN = PowerMod(r, pk.N, pk.N2);
    return MulMod(gm, rN, pk.N2);
}

BigInt Paillier::Decrypt(const PaillierSecretKey& sk, const PaillierPublicKey& pk,
                         const BigInt& c) {
    // m = L(c^lambda mod N^2) * mu mod N
    BigInt c_lambda = PowerMod(c, sk.lambda, pk.N2);
    BigInt L_val = L(c_lambda, pk.N);
    BigInt m = MulMod(L_val, sk.mu, pk.N);

    // Handle negative numbers (result should be in [0, N-1])
    if (m < 0) {
        m += pk.N;
    }

    return m;
}

BigInt Paillier::Add(const PaillierPublicKey& pk, const BigInt& c1, const BigInt& c2) {
    // Enc(m1 + m2) = Enc(m1) * Enc(m2) mod N^2
    return MulMod(c1, c2, pk.N2);
}

BigInt Paillier::ScalarMul(const PaillierPublicKey& pk, const BigInt& c, const BigInt& k) {
    // Enc(k * m) = Enc(m)^k mod N^2
    BigInt k_mod = k % pk.N;
    if (k_mod < 0) {
        k_mod += pk.N;
    }
    return PowerMod(c, k_mod, pk.N2);
}

BigInt Paillier::Negate(const PaillierPublicKey& pk, const BigInt& c) {
    // Enc(-m) = Enc(m)^{-1} mod N^2
    return InvMod(c, pk.N2);
}

BigInt Paillier::Rerandomize(const PaillierPublicKey& pk, const BigInt& c) {
    // c' = c * Enc(0) = c * r^N mod N^2
    BigInt r = RandomBigIntNonZero(pk.N);
    BigInt rN = PowerMod(r, pk.N, pk.N2);
    return MulMod(c, rN, pk.N2);
}

} // namespace otpsica
