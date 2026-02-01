#ifndef OTPSICA_GS19_H
#define OTPSICA_GS19_H

#include "common.h"
#include "crypto/paillier.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

namespace otpsica {

/**
 * Parameters for the GS19 PICT protocol.
 *
 * Based on: Ghosh, S. & Simkin, M. "The Communication Complexity of
 * Threshold Private Set Intersection." CRYPTO 2019, Section 5.
 *
 * Tests if |X ∩ Y| >= tau by checking if the symmetric set difference
 * |(X\Y) ∪ (Y\X)| <= 2t, where t = n - tau.
 */
struct GS19Params {
    size_t n;       // Set size
    size_t tau;     // Threshold: output "similar" if |X ∩ Y| >= tau
    size_t t;       // t = n - tau (half symmetric difference bound)
    size_t m;       // m = 2t + 1 (Hankel matrix dimension)
    BigInt q;       // Field prime: q > (4t^2 + 2t) * 2^kappa
    size_t kappa;   // Statistical security parameter (default 40)
};

/**
 * GS19 PICT (Private Intersection Cardinality Testing) protocol.
 *
 * Uses Hankel matrices constructed from power-sum polynomial encodings
 * of sets, encrypted under Paillier. Matrix singularity is tested via
 * the FINV protocol (KMWF07, Sections 4 & 5.1) which computes the
 * minimal polynomial of the matrix through linearly recurrent sequences
 * and Berlekamp-Massey.
 *
 * Computational cost per sender set: O(m^3) Paillier scalar multiplications
 * where m = 2(n - tau) + 1.
 */
class GS19PICT {
public:
    /**
     * Initialize protocol parameters.
     * @param n Set size
     * @param tau Threshold (|X ∩ Y| >= tau means "similar")
     * @param kappa Statistical security parameter (default 40)
     * @return GS19Params with computed t, m, q
     */
    static GS19Params InitParams(size_t n, size_t tau, size_t kappa = 40);

    /**
     * Alice's offline state after preprocessing.
     * Reusable across multiple sender sets.
     */
    struct AliceState {
        PaillierKeyPair kp;
        std::vector<BigInt> enc_diag;   // Enc(c_0)..Enc(c_{4t}), 4t+1 values
        BigInt u_val;                    // Evaluation point u (as ZZ)
        GS19Params params;
    };

    /**
     * Alice offline phase: Paillier keygen + power-sum encoding + encrypt Hankel diagonal.
     * @param X Alice's (receiver's) set
     * @param params Protocol parameters
     * @return AliceState with encrypted Hankel diagonal
     */
    static AliceState AliceOffline(const ElementVec& X, const GS19Params& params);

    /**
     * Run FINV for one sender set (Alice + Bob simulated on same machine).
     * @param state Alice's offline state
     * @param Y Sender's set
     * @return true if |X ∩ Y| >= tau
     */
    static bool RunForSender(const AliceState& state, const ElementVec& Y);

    /**
     * Run a limited number of FINV iterations (for benchmarking large params).
     * @param state Alice's offline state
     * @param Y Sender's set
     * @param max_iters Maximum number of FINV iterations to run
     * @param per_iter_ms Output: average time per iteration in ms
     * @return true (always, since we don't complete the test)
     */
    static bool RunForSenderPartial(const AliceState& state, const ElementVec& Y,
                                     size_t max_iters, double& per_iter_ms);

    /**
     * Full protocol for one pair of sets.
     */
    static bool RunProtocol(const ElementVec& X, const ElementVec& Y,
                            const GS19Params& params);

    /**
     * Multi-sender protocol: reuse Alice offline for all sender sets.
     */
    static std::vector<bool> RunProtocolMultiSender(
        const ElementVec& X, const std::vector<ElementVec>& Y_sets,
        const GS19Params& params);

    /**
     * Plaintext verification: directly compute Hankel matrix rank.
     * No encryption involved, used for correctness checking.
     */
    static bool PlaintextTest(const ElementVec& X, const ElementVec& Y,
                              const GS19Params& params);

private:
    /**
     * Compute power-sum polynomial evaluations.
     * p_S(u^i) = Σ_{a∈S} (u^a)^i for i = 0, ..., max_power
     * All arithmetic in Z_q.
     */
    static std::vector<BigInt> ComputePowerSums(
        const ElementVec& S, const BigInt& u,
        size_t max_power, const BigInt& q);

    /**
     * Encrypted Hankel matrix-vector multiplication.
     * Computes Enc(H_C · w) where H_C is stored as encrypted diagonal
     * (enc_diag[i+j] = Enc(H_C[i][j])) and w is plaintext.
     *
     * Cost: m^2 Paillier ScalarMul + m^2 Paillier Add operations.
     */
    static std::vector<BigInt> EncHankelVecMul(
        const PaillierPublicKey& pk,
        const std::vector<BigInt>& enc_diag,
        const std::vector<BigInt>& w_plain,
        size_t m);

    /** Find smallest prime > bound. */
    static BigInt FindPrime(const BigInt& bound);
};

} // namespace otpsica

#endif // OTPSICA_GS19_H
