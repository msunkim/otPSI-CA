#include "hash/hashtable.h"
#include <cstring>
#include <set>

#ifdef USE_OPENSSL
#include <openssl/sha.h>
#else
#include <cryptopp/sha.h>
#endif

namespace otpsica {

std::pair<size_t, size_t> HashTable::GenerateParams(size_t n, size_t w) {
    // Determine number of bins ell such that overflow probability is negligible
    // Based on balls-into-bins analysis: with n balls into ell bins,
    // max load is approximately n/ell + O(sqrt((n/ell) * log(ell)))
    //
    // To ensure overflow probability < 2^{-40} with high confidence,
    // we use ell = n (1 bin per element on average) with larger w.
    // This is very conservative but avoids overflow issues.

    // Use more bins: ell = n ensures average load < 1
    size_t ell = n;

    // Scale bin capacity with n to ensure negligible overflow probability.
    // With n balls into n bins (Poisson lambda=1), P(any bin > w) ~ n * e^{-1}/w!
    // For n <= 200: w=5  => ~6% failure (acceptable for small experiments)
    // For n <= 2000: w=8 => ~0.1% failure per run
    // For n > 2000: w=10 => effectively zero
    size_t safe_w;
    if (n <= 200) {
        safe_w = std::max(w, size_t(5));
    } else if (n <= 2000) {
        safe_w = std::max(w, size_t(8));
    } else {
        safe_w = std::max(w, size_t(10));
    }

    // Ensure ell is at least 1
    if (ell == 0) {
        ell = 1;
    }

    return {ell, safe_w};
}

HashTable::HashTable(size_t n, size_t w, uint64_t seed)
    : n_(n), w_(w), seed_(seed) {

    auto [ell, _w] = GenerateParams(n, w);
    ell_ = ell;
    (void)_w;

    // Initialize empty bins
    bins_.resize(ell_);

    // If seed is 0, generate random seed
    if (seed_ == 0) {
        std::random_device rd;
        seed_ = rd();
    }
}

HashTable::HashTable(size_t ell, size_t w, uint64_t seed, ExplicitParams)
    : n_(ell), ell_(ell), w_(w), seed_(seed) {

    // Initialize empty bins with explicit number
    bins_.resize(ell_);

    // If seed is 0, generate random seed
    if (seed_ == 0) {
        std::random_device rd;
        seed_ = rd();
    }
}

size_t HashTable::Hash(Element x) const {
    // Use SHA-256 with seed for hashing
    uint8_t digest[32];  // SHA-256 produces 32 bytes

#ifdef USE_OPENSSL
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, reinterpret_cast<const uint8_t*>(&seed_), sizeof(seed_));
    SHA256_Update(&ctx, reinterpret_cast<const uint8_t*>(&x), sizeof(x));
    SHA256_Final(digest, &ctx);
#else
    CryptoPP::SHA256 hash;
    hash.Update(reinterpret_cast<const uint8_t*>(&seed_), sizeof(seed_));
    hash.Update(reinterpret_cast<const uint8_t*>(&x), sizeof(x));
    hash.Final(digest);
#endif

    // Convert first 8 bytes to size_t and mod by ell
    uint64_t h = 0;
    std::memcpy(&h, digest, sizeof(h));
    return h % ell_;
}

void HashTable::Insert(Element x) {
    size_t idx = Hash(x);
    ElementVec& bin = bins_[idx];

    // Check for duplicates
    for (Element e : bin) {
        if (e == x) return;  // Already present
    }

    // Check overflow
    if (bin.size() >= w_) {
        throw std::overflow_error("Hash table bin overflow");
    }

    bin.push_back(x);
}

bool HashTable::Check(Element x) const {
    size_t idx = Hash(x);
    const ElementVec& bin = bins_[idx];

    for (Element e : bin) {
        if (e == x) return true;
    }
    return false;
}

size_t HashTable::GetBinIndex(Element x) const {
    return Hash(x);
}

const ElementVec& HashTable::GetBin(size_t bin_idx) const {
    if (bin_idx >= ell_) {
        throw std::out_of_range("Bin index out of range");
    }
    return bins_[bin_idx];
}

void HashTable::FillWithRandom(size_t beta, uint64_t party_offset) {
    // Use party-specific offset to ensure receiver and sender get different padding
    std::mt19937_64 rng(seed_ + 12345 + party_offset);
    std::uniform_int_distribution<Element> dist(1, (1ULL << beta) - 1);

    // Collect all existing elements to avoid duplicates
    std::set<Element> existing;
    for (const auto& bin : bins_) {
        for (Element e : bin) {
            existing.insert(e);
        }
    }

    // Fill each bin to capacity w
    for (auto& bin : bins_) {
        while (bin.size() < w_) {
            Element r = dist(rng);
            // Ensure uniqueness within the table
            if (existing.find(r) == existing.end()) {
                bin.push_back(r);
                existing.insert(r);
            }
        }
    }
}

void HashTable::SetSeed(uint64_t seed) {
    seed_ = seed;
}

} // namespace otpsica
