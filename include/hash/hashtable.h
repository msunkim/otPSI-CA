#ifndef OTPSICA_HASHTABLE_H
#define OTPSICA_HASHTABLE_H

#include "common.h"
#include <string>
#include <stdexcept>

#ifdef USE_OPENSSL
#include <openssl/sha.h>
#else
#include <cryptopp/sha.h>
#endif

namespace otpsica {

/**
 * Hash table for set representation.
 *
 * Structure: ell bins x w elements per bin
 * - ell: number of bins (determined by n and w to minimize overflow)
 * - w: maximum elements per bin (typically 3)
 *
 * The hash table enables efficient bin-wise polynomial operations
 * by grouping elements that hash to the same bin.
 */

class HashTable {
public:
    /**
     * Generate hash table parameters.
     * @param n Maximum number of elements
     * @param w Bin size (default 3)
     * @return Pair (ell, w) ensuring negligible overflow probability
     */
    static std::pair<size_t, size_t> GenerateParams(size_t n, size_t w = 3);

    /**
     * Construct a hash table with automatic bin calculation.
     * @param n Maximum set size
     * @param w Bin size
     * @param seed Random seed for hash function (optional)
     */
    HashTable(size_t n, size_t w = 3, uint64_t seed = 0);

    /**
     * Construct a hash table with explicit parameters.
     * Use this when both parties need identical hash table structure.
     * @param ell Number of bins
     * @param w Bin size
     * @param seed Hash function seed
     * @param explicit_ell Tag to distinguish from (n, w, seed) constructor
     */
    struct ExplicitParams {};
    HashTable(size_t ell, size_t w, uint64_t seed, ExplicitParams);

    /**
     * Insert an element into the hash table.
     * @param x Element to insert
     * @throws std::overflow_error if bin overflows
     */
    void Insert(Element x);

    /**
     * Check if an element is in the hash table.
     * @param x Element to check
     * @return true if x is in the table
     */
    bool Check(Element x) const;

    /**
     * Get the bin index for an element.
     */
    size_t GetBinIndex(Element x) const;

    /**
     * Get all elements in a bin.
     * @param bin_idx Bin index in [0, ell-1]
     * @return Vector of elements in the bin
     */
    const ElementVec& GetBin(size_t bin_idx) const;

    /**
     * Fill empty slots in all bins with random distinct elements.
     * This ensures each bin has exactly w elements (needed for FLT).
     * @param beta Element bit-size (elements drawn from [1, 2^beta])
     * @param party_offset Offset to make padding different per party (0=receiver, 1=sender)
     */
    void FillWithRandom(size_t beta, uint64_t party_offset = 0);

    /**
     * Get the number of bins.
     */
    size_t GetNumBins() const { return ell_; }

    /**
     * Get the bin capacity.
     */
    size_t GetBinSize() const { return w_; }

    /**
     * Get the hash function seed (for sharing between parties).
     */
    uint64_t GetSeed() const { return seed_; }

    /**
     * Set the hash function seed (to match another party's table).
     */
    void SetSeed(uint64_t seed);

private:
    size_t n_;                  // Maximum set size
    size_t ell_;                // Number of bins
    size_t w_;                  // Elements per bin
    uint64_t seed_;             // Hash function seed
    std::vector<ElementVec> bins_;  // The actual bins

    // Hash function: element -> bin index
    size_t Hash(Element x) const;
};

} // namespace otpsica

#endif // OTPSICA_HASHTABLE_H
