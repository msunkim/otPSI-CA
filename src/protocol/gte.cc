#include "protocol/gte.h"
#include <NTL/ZZ_pX.h>

namespace otpsica {

using namespace NTL;

std::pair<CEKKeyPair, GTESenderMessage1>
GTEProtocol::SenderRound1(const std::vector<BigInt>& m2_values,
                          const SecurityParams& sec_params) {
    // Compute d based on maximum m2 value
    // m2 = r + tau, and m1 = I + r where I <= n
    // We need d > max(m1, m2) for all values
    // Since mask is small, d = 2 * max(m2) + 128 should suffice
    BigInt max_m2 = BigInt(0);
    for (const auto& m : m2_values) {
        if (m > max_m2) max_m2 = m;
    }
    size_t d = conv<long>(max_m2) * 2 + 128;
    if (d < 20) d = 20;  // Minimum d for security

    CEKKeyPair kp = CEK::KeyGen(sec_params.kappa, sec_params.N_bits, 2, d);

    GTESenderMessage1 msg;
    msg.pk_cek = kp.pk;

    size_t eta = m2_values.size();
    msg.enc_values.resize(eta);

    for (size_t k = 0; k < eta; k++) {
        // M_{2,k} = Enc^CEK(d - m_{2,k})
        long m2_mod = conv<long>(m2_values[k] % power(kp.pk.rho, kp.pk.d));
        long val = kp.pk.d - m2_mod;
        if (val < 0) val = 0;
        msg.enc_values[k] = CEK::Encrypt(kp.pk, static_cast<size_t>(val));
    }

    return {kp, msg};
}

std::pair<GTEReceiverMessage, std::vector<BigInt>>
GTEProtocol::ReceiverRound1(const GTESenderMessage1& msg,
                            const std::vector<BigInt>& m1_values,
                            const PaillierPublicKey& pai_pk) {
    size_t eta = m1_values.size();
    GTEReceiverMessage resp;
    resp.M_values.resize(eta);

    std::vector<BigInt> sigma_values(eta);
    BigInt rho_d = CEK::GetRhoPowD(msg.pk_cek);

    for (size_t k = 0; k < eta; k++) {
        // Generate random sigma_k in [1, rho^d - 1], not divisible by rho
        do {
            sigma_values[k] = RandomBigIntNonZero(rho_d);
        } while ((sigma_values[k] % msg.pk_cek.rho) == 0);

        // Generate random s_k for blinding in G_2
        BigInt s_k = RandomBigInt(power2_ZZ(msg.pk_cek.kappa));

        // M_k = M_{2,k}^{rho^{m_{1,k}}} * g1^{sigma_k} * g2^{s_k}
        long m1_mod = conv<long>(m1_values[k] % power(msg.pk_cek.rho, msg.pk_cek.d));
        BigInt rho_m1 = power(msg.pk_cek.rho, m1_mod);

        BigInt powered = CEK::Power(msg.pk_cek, msg.enc_values[k], rho_m1);
        resp.M_values[k] = CEK::Blind(msg.pk_cek, powered, sigma_values[k], s_k);
    }

    // Build GTE indicator polynomial: f(X) = prod_{k=1}^{eta} (X - sigma_k)
    ZZ_p::init(pai_pk.N);  // Work modulo N

    ZZ_pX f;
    SetCoeff(f, 0, 1);  // Start with 1

    for (size_t k = 0; k < eta; k++) {
        ZZ_pX term;
        SetCoeff(term, 1, 1);                           // X
        SetCoeff(term, 0, -conv<ZZ_p>(sigma_values[k])); // X - sigma_k
        f = f * term;
    }

    // Encrypt polynomial coefficients
    size_t degree = deg(f);
    resp.enc_poly.resize(degree + 1);
    for (size_t i = 0; i <= degree; i++) {
        resp.enc_poly[degree - i] = Paillier::Encrypt(pai_pk, rep(coeff(f, i)));
    }

    return {resp, sigma_values};
}

GTESenderMessage2
GTEProtocol::SenderRound2(const GTEReceiverMessage& msg,
                          const CEKKeyPair& kp,
                          const PaillierPublicKey& pai_pk) {
    size_t eta = msg.M_values.size();
    GTESenderMessage2 resp;
    resp.enc_evals.resize(eta);

    for (size_t k = 0; k < eta; k++) {
        // Decrypt M_k to get v_k (raw discrete log, not rho-adic)
        long v_k = CEK::DecryptRaw(kp.sk, kp.pk, msg.M_values[k]);

        // Evaluate encrypted polynomial F at v_k
        BigInt v_k_elem = (v_k >= 0) ? conv<ZZ>(v_k) : ZZ(0);

        // F(v_k) = prod_{i=0}^{d} Enc(f_i)^{v_k^i}
        size_t poly_degree = msg.enc_poly.size() - 1;
        BigInt result = Paillier::Encrypt(pai_pk, BigInt(0));
        BigInt v_power = BigInt(1);

        for (size_t i = 0; i <= poly_degree; i++) {
            BigInt term = Paillier::ScalarMul(pai_pk, msg.enc_poly[poly_degree - i], v_power);
            result = Paillier::Add(pai_pk, result, term);
            v_power *= v_k_elem;
        }

        // Blind with random gamma_k
        BigInt gamma_k = RandomBigIntNonZero(pai_pk.N);
        resp.enc_evals[k] = Paillier::ScalarMul(pai_pk, result, gamma_k);
    }

    return resp;
}

std::vector<bool>
GTEProtocol::ReceiverRound2(const GTESenderMessage2& msg,
                            const PaillierKeyPair& pai_kp) {
    size_t eta = msg.enc_evals.size();
    std::vector<bool> results(eta);

    for (size_t k = 0; k < eta; k++) {
        // Decrypt: if result is 0, then m_{1,k} >= m_{2,k}
        BigInt decrypted = Paillier::Decrypt(pai_kp.sk, pai_kp.pk, msg.enc_evals[k]);
        results[k] = IsZero(decrypted);
    }

    return results;
}

} // namespace otpsica
