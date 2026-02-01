#ifndef OTPSICA_GTE_H
#define OTPSICA_GTE_H

#include "common.h"
#include "crypto/paillier.h"
#include "crypto/cek.h"

namespace otpsica {

/**
 * Protocol 2: Integer Multi-Comparison Protocol Pi_GTE
 *
 * Functionality F_GTE:
 *   - Receiver input: (m_{1,1}, ..., m_{1,eta})
 *   - Sender input: (m_{2,1}, ..., m_{2,eta})
 *   - Output: Receiver learns b_k = 1 if m_{1,k} >= m_{2,k}, else b_k = 0
 *             for all k in [1..eta]; Sender learns nothing
 *
 * Key idea:
 *   - Use CEK encryption for the comparison mechanism
 *   - Build GTE indicator polynomial f(X) = prod_k (X - sigma_k)
 *   - If m_{1,k} >= m_{2,k}, then v_k = sigma_k, so f(v_k) = 0
 *
 * Protocol flow:
 *   1. R <- S: CEK encryptions of (d - m_{2,k}) for each k
 *   2. R -> S: Transformed CEK values M_k and encrypted indicator polynomial F
 *   3. R <- S: Encrypted evaluations F(v_k) for each k
 *   4. R: Decrypt and output b_k = 1 if decryption is 0, else 0
 */

// Message from Sender to Receiver (Step 1)
struct GTESenderMessage1 {
    CEKPublicKey pk_cek;
    std::vector<BigInt> enc_values;  // M_{2,k} = Enc^CEK(d - m_{2,k})
};

// Message from Receiver to Sender (Step 2)
struct GTEReceiverMessage {
    std::vector<BigInt> M_values;    // M_k = M_{2,k}^{rho^{m_{1,k}}} * g1^{sigma_k} * g2^{s_k}
    std::vector<BigInt> enc_poly;    // F = Enc^Pai(f) where f(X) = prod(X - sigma_k)
};

// Message from Sender to Receiver (Step 3)
struct GTESenderMessage2 {
    std::vector<BigInt> enc_evals;   // F(v_k)^{gamma_k} for each k
};

class GTEProtocol {
public:
    /**
     * Sender: Generate CEK keypair and encrypt (d - m_{2,k}) values.
     * @param m2_values Sender's comparison values (m_{2,1}, ..., m_{2,eta})
     * @param sec_params Security parameters
     * @return CEK keypair and message
     */
    static std::pair<CEKKeyPair, GTESenderMessage1>
    SenderRound1(const std::vector<BigInt>& m2_values,
                 const SecurityParams& sec_params);

    /**
     * Receiver: Transform CEK values and build indicator polynomial.
     * @param msg Sender's message with CEK encryptions
     * @param m1_values Receiver's values (m_{1,1}, ..., m_{1,eta})
     * @param pai_pk Paillier public key (for encrypting indicator polynomial)
     * @return Message with M_k values and encrypted indicator polynomial
     */
    static std::pair<GTEReceiverMessage, std::vector<BigInt>>
    ReceiverRound1(const GTESenderMessage1& msg,
                   const std::vector<BigInt>& m1_values,
                   const PaillierPublicKey& pai_pk);
    // Returns (message, sigma_values) where sigma_values are the random blinding factors

    /**
     * Sender: Decrypt CEK values, evaluate polynomial, and blind results.
     * @param msg Receiver's message
     * @param kp CEK keypair
     * @param pai_pk Paillier public key
     * @return Blinded polynomial evaluations
     */
    static GTESenderMessage2
    SenderRound2(const GTEReceiverMessage& msg,
                 const CEKKeyPair& kp,
                 const PaillierPublicKey& pai_pk);

    /**
     * Receiver: Decrypt and determine comparison results.
     * @param msg Sender's final message
     * @param pai_kp Paillier keypair
     * @return Vector of bits: b_k = 1 if m_{1,k} >= m_{2,k}
     */
    static std::vector<bool>
    ReceiverRound2(const GTESenderMessage2& msg,
                   const PaillierKeyPair& pai_kp);
};

} // namespace otpsica

#endif // OTPSICA_GTE_H
