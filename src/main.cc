#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <set>
#include <algorithm>
#include <sys/stat.h>
#include <cmath>
#include "protocol/otpsica.h"
#include "io/setio.h"

using namespace otpsica;
using namespace std::chrono;

// Check if file exists
bool FileExists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Compute element bit-length based on set size n.
// Must satisfy: 2^beta >= 6*n (enough distinct values for hash table padding
// with ell=n bins, w=5, plus receiver/sender sets with headroom).
static size_t ComputeBeta(size_t n) {
    size_t beta = 10;  // minimum
    while ((1ULL << beta) < 6 * n) {
        beta++;
    }
    return beta;
}

// Compute actual intersection sizes for verification
std::vector<size_t> ComputeIntersections(const ElementVec& X,
                                          const std::vector<ElementVec>& Y_sets) {
    std::set<Element> X_set(X.begin(), X.end());
    std::vector<size_t> sizes(Y_sets.size());

    for (size_t k = 0; k < Y_sets.size(); k++) {
        size_t count = 0;
        for (Element y : Y_sets[k]) {
            if (X_set.find(y) != X_set.end()) {
                count++;
            }
        }
        sizes[k] = count;
    }

    return sizes;
}

int main(int argc, char* argv[]) {
    std::cout << "=== otPSI-CA Protocol Implementation ===" << std::endl;
    std::cout << std::endl;

    // Default parameters
    // n = 100 (fixed for experiments)
    // eta in {10, 20, 30}
    // tau in {50, 60}
    size_t n = 100;      // Set size
    size_t eta = 10;     // Number of sender's sets (from {10, 20, 30})
    size_t tau = 50;     // Threshold (from {50, 60})

    if (argc > 1) n = std::stoull(argv[1]);
    if (argc > 2) eta = std::stoull(argv[2]);
    if (argc > 3) tau = std::stoull(argv[3]);

    size_t BETA = ComputeBeta(n);

    std::cout << "Parameters:" << std::endl;
    std::cout << "  Set size (n):     " << n << std::endl;
    std::cout << "  Number of sets (eta): " << eta << std::endl;
    std::cout << "  Threshold (tau):  " << tau << std::endl;
    std::cout << "  Element bits (beta): " << BETA << std::endl;
    std::cout << std::endl;

    // File paths for storing/loading sets
    std::string data_dir = "data";
    std::string rx_file = data_dir + "/receiver_n" + std::to_string(n) + ".bin";
    std::string tx_file = data_dir + "/sender_n" + std::to_string(n) + "_eta" + std::to_string(eta) + ".bin";

    ElementVec X;
    std::vector<ElementVec> Y_sets;

    // Check if sets already exist on disk
    if (FileExists(rx_file) && FileExists(tx_file)) {
        std::cout << "Loading existing sets from files..." << std::endl;
        X = SetIO::ReadSet(rx_file);
        Y_sets = SetIO::ReadSets(tx_file);

        // Verify loaded data matches expected sizes
        if (X.size() != n) {
            std::cerr << "Warning: Loaded receiver set size (" << X.size()
                      << ") doesn't match n (" << n << "). Regenerating..." << std::endl;
            X.clear();
        }
        if (Y_sets.size() != eta) {
            std::cerr << "Warning: Loaded sender sets count (" << Y_sets.size()
                      << ") doesn't match eta (" << eta << "). Regenerating..." << std::endl;
            Y_sets.clear();
        }
    }

    // Generate sets if not loaded from file
    if (X.empty() || Y_sets.empty()) {
        std::cout << "Generating test sets with varied intersection sizes..." << std::endl;

        // Generate receiver set
        RandomSetGenerator gen(42);
        X = gen.Generate(BETA, n);

        // Generate sender sets with controlled intersection sizes:
        // Alternate between above-threshold and below-threshold intersections
        std::vector<size_t> int_sizes(eta);
        for (size_t k = 0; k < eta; k++) {
            if (k % 2 == 0) {
                // Above threshold: I = tau + 10 (or n if that's smaller)
                int_sizes[k] = std::min(n, tau + 10);
            } else {
                // Below threshold: I = tau - 10 (or 0 if negative)
                int_sizes[k] = (tau >= 10) ? tau - 10 : 0;
            }
        }
        Y_sets = gen.GenerateWithIntersection(X, BETA, n, eta, int_sizes);

        // Create data directory if it doesn't exist
        #ifdef _WIN32
        _mkdir(data_dir.c_str());
        #else
        mkdir(data_dir.c_str(), 0755);
        #endif

        // Save sets to files for future runs
        std::cout << "Saving sets to files for reproducibility..." << std::endl;
        SetIO::WriteSet(rx_file, X, BETA);
        SetIO::WriteSets(tx_file, Y_sets, BETA);
        std::cout << "  Receiver set saved to: " << rx_file << std::endl;
        std::cout << "  Sender sets saved to: " << tx_file << std::endl;
    } else {
        std::cout << "  Receiver set loaded from: " << rx_file << std::endl;
        std::cout << "  Sender sets loaded from: " << tx_file << std::endl;
    }

    // Compute actual intersections for verification
    auto actual_sizes = ComputeIntersections(X, Y_sets);
    std::cout << "Actual intersection sizes: ";
    for (size_t k = 0; k < std::min(eta, size_t(5)); k++) {
        std::cout << actual_sizes[k] << " ";
    }
    if (eta > 5) std::cout << "...";
    std::cout << std::endl << std::endl;

    // Initialize protocol parameters
    ProtocolParams params = OtPSICA::InitParams(n, eta, tau, BETA);
    std::cout << "Protocol parameters:" << std::endl;
    std::cout << "  Bins (ell):       " << params.ell << std::endl;
    std::cout << "  Bin size (w):     " << params.w << std::endl;
    std::cout << "  CRT primes:       ";
    for (size_t p : params.primes) std::cout << p << " ";
    std::cout << std::endl << std::endl;

    // Run protocol
    std::cout << "Running protocol..." << std::endl;
    auto start = high_resolution_clock::now();

    std::vector<bool> results = OtPSICA::RunProtocol(X, Y_sets, tau, params);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    std::cout << "Protocol completed in " << duration.count() << " ms" << std::endl;
    std::cout << std::endl;

    // Display results
    std::cout << "Results (b_k = 1 if |X ∩ Y^(k)| >= " << tau << "):" << std::endl;
    size_t correct = 0;
    for (size_t k = 0; k < eta; k++) {
        bool expected = (actual_sizes[k] >= tau);
        bool actual = results[k];
        bool match = (expected == actual);
        if (match) correct++;

        std::cout << "  Y^(" << (k+1) << "): |X ∩ Y| = " << actual_sizes[k]
                  << ", b_k = " << actual
                  << (match ? " ✓" : " ✗") << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Correctness: " << correct << "/" << eta << std::endl;

    return (correct == eta) ? 0 : 1;
}
