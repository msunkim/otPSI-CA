#ifndef OTPSICA_PROTOCOL_H
#define OTPSICA_PROTOCOL_H

#include "common.h"
#include "protocol/smt.h"
#include "protocol/gte.h"

namespace otpsica {

/**
 * Protocol 3: Over-threshold Private Set Intersection Cardinality (otPSI-CA*)
 *
 * Functionality F_otPSI-CA*:
 *   - Receiver input: Set X
 *   - Sender input: Sets Y^(1), ..., Y^(eta)
 *   - Public input: Threshold tau
 *   - Output: Receiver learns b_k = 1 if |X ∩ Y^(k)| >= tau, else b_k = 0
 *             for all k in [1..eta]; Sender learns nothing
 *
 * Protocol structure:
 *   Phase 1 (Oblivious Count): Use Pi_SMT to compute masked cardinalities
 *   Phase 2 (Oblivious Comparison): Use Pi_GTE to compare with threshold
 *
 * Round complexity: 2 rounds (constant)
 * Computation: O(n) for receiver, O(eta * n) for sender
 * Communication: O(n + eta) ciphertexts
 */

// Complete protocol state for Receiver
struct ReceiverState {
    PaillierKeyPair pai_kp;
    ProtocolParams params;
    std::vector<BigInt> m1_values;     // Masked cardinalities
    std::vector<BigInt> sigma_values;  // GTE blinding factors
};

// Complete protocol state for Sender
struct SenderState {
    CEKKeyPair cek_kp;
    std::vector<BigInt> r_values;      // Random masks from Phase 1
    std::vector<BigInt> m2_values;     // r_k + tau values
};

// Round 1 message: Receiver -> Sender
struct Round1Message {
    SMTReceiverMessage smt_msg;
};

// Round 1 response: Sender -> Receiver
struct Round1Response {
    SMTSenderMessage smt_resp;
    GTESenderMessage1 gte_msg1;
};

// Round 2 message: Receiver -> Sender
struct Round2Message {
    GTEReceiverMessage gte_msg;
};

// Round 2 response: Sender -> Receiver
struct Round2Response {
    GTESenderMessage2 gte_msg2;
};

class OtPSICA {
public:
    /**
     * Initialize protocol parameters.
     * @param n Maximum set size
     * @param eta Number of sender's sets
     * @param tau Threshold value
     * @param beta Element bit-size (default 16)
     * @return Configured protocol parameters
     */
    static ProtocolParams InitParams(size_t n, size_t eta, size_t tau,
                                      size_t beta = 16);

    // ===== RECEIVER INTERFACE =====

    /**
     * Receiver Round 1: Generate encrypted polynomials.
     * @param X Receiver's set
     * @param params Protocol parameters
     * @return (ReceiverState, Round1Message)
     */
    static std::pair<ReceiverState, Round1Message>
    ReceiverRound1(const ElementVec& X, const ProtocolParams& params);

    /**
     * Receiver processes Round 1 response and generates Round 2 message.
     * @param state Receiver state from Round 1
     * @param resp Sender's Round 1 response
     * @return Round2Message
     */
    static Round2Message
    ReceiverProcessRound1(ReceiverState& state, const Round1Response& resp);

    /**
     * Receiver final: Decrypt and output results.
     * @param state Receiver state
     * @param resp Sender's Round 2 response
     * @return Vector of results: b_k for k in [1..eta]
     */
    static std::vector<bool>
    ReceiverFinal(const ReceiverState& state, const Round2Response& resp);

    // ===== SENDER INTERFACE =====

    /**
     * Sender Round 1: Process receiver's message, compute masked cardinalities.
     * @param msg Receiver's Round 1 message
     * @param sets Sender's sets Y^(1), ..., Y^(eta)
     * @param tau Threshold value
     * @param sec_params Security parameters
     * @return (SenderState, Round1Response)
     */
    static std::pair<SenderState, Round1Response>
    SenderRound1(const Round1Message& msg,
                 const std::vector<ElementVec>& sets,
                 size_t tau,
                 const SecurityParams& sec_params);

    /**
     * Sender Round 2: Process comparison request.
     * @param state Sender state from Round 1
     * @param msg Receiver's Round 2 message
     * @return Round2Response
     */
    static Round2Response
    SenderRound2(const SenderState& state, const Round2Message& msg);

    // ===== CONVENIENCE FUNCTIONS =====

    /**
     * Run the complete protocol (for testing).
     * @param X Receiver's set
     * @param Y_sets Sender's sets
     * @param tau Threshold
     * @param params Protocol parameters
     * @return Vector of results b_k
     */
    static std::vector<bool>
    RunProtocol(const ElementVec& X,
                const std::vector<ElementVec>& Y_sets,
                size_t tau,
                const ProtocolParams& params);
};

} // namespace otpsica

#endif // OTPSICA_PROTOCOL_H
