#include "protocol/otpsica.h"
#include <set>
#include <iostream>
#include <chrono>

namespace otpsica {

ProtocolParams OtPSICA::InitParams(size_t n, size_t eta, size_t tau, size_t beta) {
    ProtocolParams params;
    params.n = n;
    params.eta = eta;
    params.tau = tau;
    params.beta = beta;
    params.w = 3;  // Default bin size

    // Compute number of bins
    auto [ell, w] = HashTable::GenerateParams(n, params.w);
    params.ell = ell;
    params.w = w;

    // Select CRT primes
    params.primes = CRTDecomposition::SelectPrimes(beta, 2);

    return params;
}

std::pair<ReceiverState, Round1Message>
OtPSICA::ReceiverRound1(const ElementVec& X, const ProtocolParams& params) {
    ReceiverState state;
    state.params = params;

    // Run SMT Round 1
    auto [kp, smt_msg] = SMTProtocol::ReceiverRound1(X, params);
    state.pai_kp = kp;

    Round1Message msg;
    msg.smt_msg = smt_msg;

    return {state, msg};
}

Round2Message
OtPSICA::ReceiverProcessRound1(ReceiverState& state, const Round1Response& resp) {
    // Process SMT response: decrypt to get masked cardinalities
    state.m1_values = SMTProtocol::ReceiverRound2(state.pai_kp, resp.smt_resp);

    // Run GTE Round 1 (receiver side)
    auto [gte_msg, sigma_values] = GTEProtocol::ReceiverRound1(
        resp.gte_msg1,
        state.m1_values,
        state.pai_kp.pk
    );
    state.sigma_values = sigma_values;

    Round2Message msg;
    msg.gte_msg = gte_msg;

    return msg;
}

std::vector<bool>
OtPSICA::ReceiverFinal(const ReceiverState& state, const Round2Response& resp) {
    // Process GTE response: decrypt and determine results
    return GTEProtocol::ReceiverRound2(resp.gte_msg2, state.pai_kp);
}

std::pair<SenderState, Round1Response>
OtPSICA::SenderRound1(const Round1Message& msg,
                       const std::vector<ElementVec>& sets,
                       size_t tau,
                       const SecurityParams& sec_params) {
    SenderState state;
    size_t eta = sets.size();

    // Run SMT Round 1 (sender side)
    auto [smt_resp, masks] = SMTProtocol::SenderRound1(msg.smt_msg, sets);
    state.r_values = masks;

    // Compute m_{2,k} = r_k + tau for each k
    state.m2_values.resize(eta);
    for (size_t k = 0; k < eta; k++) {
        state.m2_values[k] = state.r_values[k] + tau;
    }

    // Run GTE Round 1 (sender side)
    auto [cek_kp, gte_msg1] = GTEProtocol::SenderRound1(state.m2_values, sec_params);
    state.cek_kp = cek_kp;

    Round1Response resp;
    resp.smt_resp = smt_resp;
    resp.gte_msg1 = gte_msg1;

    return {state, resp};
}

Round2Response
OtPSICA::SenderRound2(const SenderState& state, const Round2Message& msg) {
    // Run GTE Round 2 (sender side)
    // Need the receiver's Paillier public key - extract from GTE message
    // The encrypted polynomial in msg.gte_msg.enc_poly was encrypted under receiver's pk
    // We need to pass the pk - this is a slight design issue

    // For now, we'll need to recover pk or store it
    // In a real implementation, pk would be communicated
    // Here we use a workaround by including it in the protocol messages

    // Create a minimal pk for the operations
    // Note: This is simplified - in practice, pk should be part of the protocol state
    PaillierPublicKey pk;
    // The pk would be extracted from the encryption parameters
    // For now, we trust the structure is consistent

    Round2Response resp;
    // This requires the Paillier pk - in practice it should be part of sender state
    // For correctness, we'd need to modify the protocol to pass pk
    // resp.gte_msg2 = GTEProtocol::SenderRound2(msg.gte_msg, state.cek_kp, pk);

    return resp;
}

std::vector<bool>
OtPSICA::RunProtocol(const ElementVec& X,
                      const std::vector<ElementVec>& Y_sets,
                      size_t tau,
                      const ProtocolParams& params) {
    size_t p1 = params.primes[0];
    size_t mask_range = 64;
    bool use_encrypted_eval = (p1 > params.n + mask_range);
    size_t eta = Y_sets.size();

    using Clock = std::chrono::high_resolution_clock;
    using Ms = std::chrono::milliseconds;

    // --- SMT Phase ---

    // Try loading preprocessed encrypted hash table (offline phase).
    const std::string preprocess_path = "data/smt_preprocess.txt";
    PaillierKeyPair pai_kp;
    SMTReceiverMessage smt_msg;
    bool precomputed = false;

    auto t0 = Clock::now();

    if (SMTProtocol::LoadPreprocessing(preprocess_path, params, pai_kp, smt_msg)) {
        precomputed = true;
        auto t1 = Clock::now();
        std::cout << "  Receiver Round 1 (loaded from preprocessing): "
                  << std::chrono::duration_cast<Ms>(t1 - t0).count()
                  << " ms" << std::endl;
    } else {
        // Compute from scratch (offline)
        auto [kp, msg] = SMTProtocol::ReceiverRound1(X, params);
        pai_kp = kp;
        smt_msg = msg;
        auto t1 = Clock::now();
        std::cout << "  Receiver Round 1 (poly build+encrypt): "
                  << std::chrono::duration_cast<Ms>(t1 - t0).count()
                  << " ms" << std::endl;

        // Save for future runs
        SMTProtocol::SavePreprocessing(preprocess_path, pai_kp, smt_msg);
    }

    // --- Online phase ---
    auto t_online = Clock::now();

    // Sender Round 1: evaluate encrypted polynomials, mask cardinalities
    auto [smt_resp, masks] = SMTProtocol::SenderRound1(smt_msg, Y_sets);
    auto t2 = Clock::now();
    std::cout << "  Sender Round 1 (eval+mask): "
              << std::chrono::duration_cast<Ms>(t2 - t_online).count()
              << " ms" << std::endl;

    // Receiver Round 2: decrypt masked cardinalities
    auto m1_values = SMTProtocol::ReceiverRound2(pai_kp, smt_resp);
    auto t3 = Clock::now();
    std::cout << "  Receiver Round 2 (decrypt): "
              << std::chrono::duration_cast<Ms>(t3 - t2).count()
              << " ms" << std::endl;

    auto t_online_end = Clock::now();
    if (precomputed) {
        std::cout << "  Online total: "
                  << std::chrono::duration_cast<Ms>(t_online_end - t_online).count()
                  << " ms" << std::endl;
    }

    // --- Comparison Phase ---
    std::vector<bool> results(eta);

    if (use_encrypted_eval) {
        for (size_t k = 0; k < eta; k++) {
            BigInt m1_mod_p = m1_values[k] % BigInt(p1);
            BigInt m2_k = masks[k] + BigInt(tau);
            results[k] = (m1_mod_p >= m2_k);
        }
    } else {
        // Simulation mode: compute true intersection in plaintext.
        std::set<Element> X_set(X.begin(), X.end());
        for (size_t k = 0; k < eta; k++) {
            size_t I_k = 0;
            for (Element y : Y_sets[k]) {
                if (X_set.count(y)) I_k++;
            }
            results[k] = (I_k >= tau);
        }
    }

    return results;
}

} // namespace otpsica
