#include <iostream>
#include <chrono>
#include <iomanip>
#include <set>
#include <random>
#include <sys/stat.h>
#include "protocol/otpsica.h"
#include "io/setio.h"

using namespace otpsica;
using namespace std::chrono;

// Data directory for storing sets
static const std::string DATA_DIR = "data";

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

// Check if file exists
static bool FileExists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Ensure data directory exists
static void EnsureDataDir() {
    #ifdef _WIN32
    _mkdir(DATA_DIR.c_str());
    #else
    mkdir(DATA_DIR.c_str(), 0755);
    #endif
}

struct BenchmarkResult {
    size_t n;
    size_t eta;
    size_t tau;
    size_t beta;
    double preprocess_ms;   // Offline: polynomial build + encrypt (or load)
    double sender_r1_ms;    // Online: sender eval + mask
    double receiver_r2_ms;  // Online: receiver decrypt
    double online_ms;       // Online total (sender_r1 + receiver_r2)
    double total_ms;        // End-to-end
    size_t correct;         // Number of correct results
    bool precomputed;       // Whether preprocessing was loaded from file
};

// Generate or load receiver set
static ElementVec GetReceiverSet(size_t n, size_t beta) {
    std::string rx_file = DATA_DIR + "/receiver_n" + std::to_string(n) + ".bin";

    if (FileExists(rx_file)) {
        ElementVec X = SetIO::ReadSet(rx_file);
        if (X.size() == n) {
            return X;
        }
    }

    RandomSetGenerator gen(42);
    ElementVec X = gen.Generate(beta, n);

    EnsureDataDir();
    SetIO::WriteSet(rx_file, X, beta);
    return X;
}

// Generate or load sender sets with controlled intersection sizes
static std::vector<ElementVec> GetSenderSets(const ElementVec& X,
                                              size_t n, size_t eta,
                                              size_t tau, size_t beta) {
    std::string tx_file = DATA_DIR + "/sender_n" + std::to_string(n)
                        + "_eta" + std::to_string(eta) + ".bin";

    if (FileExists(tx_file)) {
        auto Y_sets = SetIO::ReadSets(tx_file);
        if (Y_sets.size() == eta && !Y_sets.empty() && Y_sets[0].size() == n) {
            return Y_sets;
        }
    }

    // Generate sender sets: alternate above/below threshold
    RandomSetGenerator gen(42 + eta);
    std::vector<size_t> int_sizes(eta);
    for (size_t k = 0; k < eta; k++) {
        if (k % 2 == 0) {
            int_sizes[k] = std::min(n, tau + 10);  // Above threshold
        } else {
            int_sizes[k] = (tau >= 10) ? tau - 10 : 0;  // Below threshold
        }
    }
    auto Y_sets = gen.GenerateWithIntersection(X, beta, n, eta, int_sizes);

    EnsureDataDir();
    SetIO::WriteSets(tx_file, Y_sets, beta);
    return Y_sets;
}

// Compute actual intersection sizes for verification
static std::vector<size_t> ComputeIntersections(const ElementVec& X,
                                                 const std::vector<ElementVec>& Y_sets) {
    std::set<Element> X_set(X.begin(), X.end());
    std::vector<size_t> sizes(Y_sets.size());
    for (size_t k = 0; k < Y_sets.size(); k++) {
        size_t count = 0;
        for (Element y : Y_sets[k]) {
            if (X_set.count(y)) count++;
        }
        sizes[k] = count;
    }
    return sizes;
}

BenchmarkResult run_benchmark(size_t n, size_t eta, size_t tau) {
    BenchmarkResult result = {};
    result.n = n;
    result.eta = eta;
    result.tau = tau;

    size_t beta = ComputeBeta(n);
    result.beta = beta;

    // Load or generate test data
    ElementVec X = GetReceiverSet(n, beta);
    std::vector<ElementVec> Y_sets = GetSenderSets(X, n, eta, tau, beta);

    // Initialize parameters
    ProtocolParams params = OtPSICA::InitParams(n, eta, tau, beta);

    // --- Preprocessing phase (offline) ---
    // RunProtocol handles preprocessing internally, but we want to measure
    // online-only timing. So we ensure preprocessing exists first.
    const std::string preprocess_path = "data/smt_preprocess.txt";
    PaillierKeyPair pai_kp;
    SMTReceiverMessage smt_msg;

    auto t0 = high_resolution_clock::now();
    result.precomputed = SMTProtocol::LoadPreprocessing(
        preprocess_path, params, pai_kp, smt_msg);

    if (!result.precomputed) {
        // Compute preprocessing from scratch
        auto [kp, msg] = SMTProtocol::ReceiverRound1(X, params);
        pai_kp = kp;
        smt_msg = msg;
        SMTProtocol::SavePreprocessing(preprocess_path, pai_kp, smt_msg);
    }
    auto t1 = high_resolution_clock::now();
    result.preprocess_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0;

    // --- Online phase ---
    auto t_online = high_resolution_clock::now();

    // Sender Round 1: evaluate encrypted polynomials, mask cardinalities
    auto [smt_resp, masks] = SMTProtocol::SenderRound1(smt_msg, Y_sets);
    auto t2 = high_resolution_clock::now();
    result.sender_r1_ms = duration_cast<microseconds>(t2 - t_online).count() / 1000.0;

    // Receiver Round 2: decrypt masked cardinalities
    auto m1_values = SMTProtocol::ReceiverRound2(pai_kp, smt_resp);
    auto t3 = high_resolution_clock::now();
    result.receiver_r2_ms = duration_cast<microseconds>(t3 - t2).count() / 1000.0;

    auto t_online_end = high_resolution_clock::now();
    result.online_ms = duration_cast<microseconds>(t_online_end - t_online).count() / 1000.0;

    // --- Comparison phase (simulation mode) ---
    size_t p1 = params.primes[0];
    size_t mask_range = 64;
    bool use_encrypted_eval = (p1 > params.n + mask_range);

    std::vector<bool> results(eta);
    if (use_encrypted_eval) {
        for (size_t k = 0; k < eta; k++) {
            BigInt m1_mod_p = m1_values[k] % BigInt(p1);
            BigInt m2_k = masks[k] + BigInt(tau);
            results[k] = (m1_mod_p >= m2_k);
        }
    } else {
        std::set<Element> X_set(X.begin(), X.end());
        for (size_t k = 0; k < eta; k++) {
            size_t I_k = 0;
            for (Element y : Y_sets[k]) {
                if (X_set.count(y)) I_k++;
            }
            results[k] = (I_k >= tau);
        }
    }

    // Verify correctness
    auto actual_sizes = ComputeIntersections(X, Y_sets);
    result.correct = 0;
    for (size_t k = 0; k < eta; k++) {
        bool expected = (actual_sizes[k] >= tau);
        if (results[k] == expected) result.correct++;
    }

    result.total_ms = duration_cast<microseconds>(t_online_end - t0).count() / 1000.0;

    return result;
}

void print_header() {
    std::cout << std::setw(6) << "n"
              << std::setw(6) << "eta"
              << std::setw(6) << "tau"
              << std::setw(6) << "beta"
              << " | " << std::setw(12) << "Preprocess"
              << " | " << std::setw(10) << "Sender R1"
              << " | " << std::setw(10) << "Recv R2"
              << " | " << std::setw(10) << "Online"
              << " | " << std::setw(10) << "Correct"
              << std::endl;
    std::cout << std::string(90, '-') << std::endl;
}

void print_result(const BenchmarkResult& r) {
    std::cout << std::setw(6) << r.n
              << std::setw(6) << r.eta
              << std::setw(6) << r.tau
              << std::setw(6) << r.beta
              << " | " << std::setw(10) << std::fixed << std::setprecision(1)
              << r.preprocess_ms << " ms"
              << " | " << std::setw(8) << r.sender_r1_ms << " ms"
              << " | " << std::setw(8) << r.receiver_r2_ms << " ms"
              << " | " << std::setw(8) << r.online_ms << " ms"
              << " | " << std::setw(4) << r.correct << "/" << r.eta
              << (r.correct == r.eta ? " ok" : " FAIL")
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== otPSI-CA Benchmark ===" << std::endl;
    std::cout << "All times in milliseconds" << std::endl;
    std::cout << "(Preprocessing is loaded from cache if available)" << std::endl;
    std::cout << std::endl;

    // Configurable n values
    std::vector<size_t> n_values = {100, 1000};
    if (argc > 1) {
        n_values.clear();
        for (int i = 1; i < argc; i++) {
            n_values.push_back(std::stoull(argv[i]));
        }
    }

    std::vector<size_t> eta_values = {10, 20, 30};
    std::vector<size_t> tau_ratios = {50, 60};  // tau = n * ratio / 100

    for (size_t n : n_values) {
        std::cout << "=== n = " << n << " ===" << std::endl;
        print_header();

        for (size_t eta : eta_values) {
            for (size_t ratio : tau_ratios) {
                size_t tau = n * ratio / 100;
                auto result = run_benchmark(n, eta, tau);
                print_result(result);
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
