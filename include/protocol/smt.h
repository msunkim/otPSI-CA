#ifndef OTPSICA_SMT_H
#define OTPSICA_SMT_H

#include "common.h"
#include "crypto/paillier.h"
#include "poly/flt_poly.h"
#include "hash/hashtable.h"

namespace otpsica {

/**
 * Protocol 1: Set Membership Test (SMT) Protocol Pi_SMT
 *
 * Functionality F_SMT:
 *   - Receiver input: Set X (as hash table bins)
 *   - Sender input: Element y
 *   - Output: Receiver learns 1 if y in X, else 0; Sender learns nothing
 *
 * Protocol flow:
 *   1. R -> S: Encrypted FLT polynomial F = Enc(f) for each bin
 *   2. R <- S: Encrypted evaluation F(y) for sender's element
 *   3. R: Decrypt and output result
 */

// Message from Receiver to Sender (Step 1)
struct SMTReceiverMessage {
    PaillierPublicKey pk;
    // enc_polys[i][j] = encrypted polynomial for bin i, prime p_j
    // Each polynomial has degree w(p_j - 1)
    std::vector<std::vector<std::vector<BigInt>>> enc_polys;
    size_t n;    // Set size (for mask range)
    size_t tau;  // Threshold (for mask range)
    size_t ell;  // Number of bins
    size_t w;    // Bin size
    size_t beta; // Element bit-length
    std::vector<size_t> primes;  // CRT primes
    uint64_t hash_seed;  // Hash function seed
};

// Message from Sender to Receiver (Step 2)
struct SMTSenderMessage {
    // For each sender's set k: aggregated encrypted cardinality + mask
    std::vector<BigInt> masked_cardinalities;  // T^(k) values
};

class SMTProtocol {
public:
    /**
     * Receiver: Generate the first message (encrypted polynomials).
     * @param set Receiver's set X
     * @param params Protocol parameters
     * @return Keypair and message to send
     */
    static std::pair<PaillierKeyPair, SMTReceiverMessage>
    ReceiverRound1(const ElementVec& set, const ProtocolParams& params);

    /**
     * Sender: Process receiver's message and compute masked cardinalities.
     * @param msg Receiver's message
     * @param sets Sender's sets Y^(1), ..., Y^(eta)
     * @return Masked cardinality ciphertexts and random masks
     */
    static std::pair<SMTSenderMessage, std::vector<BigInt>>
    SenderRound1(const SMTReceiverMessage& msg,
                 const std::vector<ElementVec>& sets);

    /**
     * Receiver: Decrypt to get masked cardinalities.
     * @param kp Keypair from Round 1
     * @param msg Sender's message
     * @return Masked cardinalities m_{1,k} = I_k + r_k
     */
    static std::vector<BigInt>
    ReceiverRound2(const PaillierKeyPair& kp,
                   const SMTSenderMessage& msg);

    // ----- Preprocessing (offline/online split) -----

    /**
     * Save Receiver Round 1 output (encrypted polynomials + keys) to file.
     * Enables offline preprocessing: compute once, reuse across sessions.
     * @param path Output file path
     * @param kp Paillier key pair
     * @param msg Receiver's first-round message
     */
    static void SavePreprocessing(const std::string& path,
                                  const PaillierKeyPair& kp,
                                  const SMTReceiverMessage& msg);

    /**
     * Load Receiver Round 1 output from a preprocessing file.
     * @param path Input file path
     * @param params Expected protocol parameters (for validation)
     * @param kp [out] Loaded Paillier key pair
     * @param msg [out] Loaded receiver message
     * @return true if loaded successfully and parameters match
     */
    static bool LoadPreprocessing(const std::string& path,
                                  const ProtocolParams& params,
                                  PaillierKeyPair& kp,
                                  SMTReceiverMessage& msg);
};

} // namespace otpsica

#endif // OTPSICA_SMT_H
