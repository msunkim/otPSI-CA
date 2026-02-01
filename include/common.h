#ifndef OTPSICA_COMMON_H
#define OTPSICA_COMMON_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <vector>
#include <cstdint>
#include <random>
#include <set>
#include <algorithm>

namespace otpsica {

// Type aliases
using BigInt = NTL::ZZ;
using Element = uint64_t;
using ElementVec = std::vector<Element>;

// Security parameters
struct SecurityParams {
    size_t lambda = 128;        // Security parameter (bits)
    size_t N_bits = 2048;       // RSA modulus size
    size_t kappa = 256;         // DLog security parameter for CEK
};

// Protocol parameters
struct ProtocolParams {
    size_t n;                   // Maximum set size
    size_t eta;                 // Number of sender's sets
    size_t tau;                 // Threshold
    size_t beta = 16;           // Element bit-size
    size_t w = 3;               // Bin size
    size_t ell;                 // Number of bins (computed from n, w)
    std::vector<size_t> primes; // CRT primes {p_1, ..., p_nu}
};

// Random number generation
inline BigInt RandomBigInt(const BigInt& max) {
    BigInt r;
    NTL::RandomBnd(r, max);
    return r;
}

inline BigInt RandomBigIntNonZero(const BigInt& max) {
    BigInt r;
    do {
        NTL::RandomBnd(r, max);
    } while (NTL::IsZero(r));
    return r;
}

/**
 * Random set generator for protocol simulation and benchmarking.
 *
 * Generates sets of random elements with specified bit-length and size.
 * Useful for simulating receiver and sender sets with variable parameters.
 */
class RandomSetGenerator {
public:
    /**
     * Construct a generator with given seed.
     * @param seed Random seed (0 = use random device)
     */
    explicit RandomSetGenerator(uint64_t seed = 0) {
        if (seed == 0) {
            std::random_device rd;
            seed = rd();
        }
        rng_.seed(seed);
    }

    /**
     * Generate a random set of distinct elements.
     * @param beta Bit-length of elements (elements in range [1, 2^beta - 1])
     * @param size Number of elements to generate
     * @return Vector of distinct random elements
     */
    ElementVec Generate(size_t beta, size_t size) {
        ElementVec result;
        result.reserve(size);

        // Maximum element value: 2^beta - 1
        Element max_val = (beta >= 64) ? ~Element(0) : ((Element(1) << beta) - 1);
        std::uniform_int_distribution<Element> dist(1, max_val);

        std::set<Element> seen;
        while (result.size() < size) {
            Element e = dist(rng_);
            if (seen.find(e) == seen.end()) {
                result.push_back(e);
                seen.insert(e);
            }
        }

        return result;
    }

    /**
     * Generate multiple random sets (for sender).
     * @param beta Bit-length of elements
     * @param size Number of elements per set
     * @param num_sets Number of sets to generate
     * @return Vector of sets
     */
    std::vector<ElementVec> GenerateMultiple(size_t beta, size_t size, size_t num_sets) {
        std::vector<ElementVec> result(num_sets);
        for (size_t k = 0; k < num_sets; k++) {
            result[k] = Generate(beta, size);
        }
        return result;
    }

    /**
     * Generate sender sets with controlled intersection with receiver set.
     * @param receiver_set The receiver's set X
     * @param beta Bit-length of elements
     * @param size Number of elements per sender set
     * @param num_sets Number of sender sets
     * @param intersection_sizes Vector of desired intersection sizes for each set
     * @return Vector of sender sets with specified intersection sizes
     */
    std::vector<ElementVec> GenerateWithIntersection(
            const ElementVec& receiver_set,
            size_t beta,
            size_t size,
            size_t num_sets,
            const std::vector<size_t>& intersection_sizes) {

        std::vector<ElementVec> result(num_sets);
        Element max_val = (beta >= 64) ? ~Element(0) : ((Element(1) << beta) - 1);
        std::uniform_int_distribution<Element> dist(1, max_val);

        // Create a copy of receiver set for shuffling
        ElementVec rx_copy = receiver_set;

        // Set of all receiver elements for exclusion
        std::set<Element> rx_set(receiver_set.begin(), receiver_set.end());

        for (size_t k = 0; k < num_sets; k++) {
            std::set<Element> seen;
            ElementVec& Y = result[k];
            Y.reserve(size);

            // Determine intersection size for this set
            size_t int_size = (k < intersection_sizes.size())
                              ? intersection_sizes[k]
                              : 0;
            int_size = std::min(int_size, std::min(size, receiver_set.size()));

            // Shuffle and pick int_size elements from receiver set
            std::shuffle(rx_copy.begin(), rx_copy.end(), rng_);
            for (size_t i = 0; i < int_size; i++) {
                Y.push_back(rx_copy[i]);
                seen.insert(rx_copy[i]);
            }

            // Fill remaining with random elements NOT in receiver set
            while (Y.size() < size) {
                Element e = dist(rng_);
                if (seen.find(e) == seen.end() && rx_set.find(e) == rx_set.end()) {
                    Y.push_back(e);
                    seen.insert(e);
                }
            }
        }

        return result;
    }

    /**
     * Get the underlying RNG for custom operations.
     */
    std::mt19937_64& GetRng() { return rng_; }

private:
    std::mt19937_64 rng_;
};

} // namespace otpsica

#endif // OTPSICA_COMMON_H
