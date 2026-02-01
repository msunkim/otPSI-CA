#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <set>
#include <random>
#include <sys/stat.h>
#include "protocol/gs19.h"
#include "io/setio.h"

using namespace otpsica;
using namespace std::chrono;

static const std::string DATA_DIR = "data";

// Compute element bit-length based on set size n.
static size_t ComputeBeta(size_t n) {
    size_t beta = 10;
    while ((1ULL << beta) < 6 * n) {
        beta++;
    }
    return beta;
}

static bool FileExists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

static void EnsureDataDir() {
    mkdir(DATA_DIR.c_str(), 0755);
}

// Generate or load receiver set
static ElementVec GetReceiverSet(size_t n, size_t beta) {
    std::string rx_file = DATA_DIR + "/receiver_n" + std::to_string(n) + ".bin";
    if (FileExists(rx_file)) {
        ElementVec X = SetIO::ReadSet(rx_file);
        if (X.size() == n) return X;
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
    RandomSetGenerator gen(42 + eta);
    std::vector<size_t> int_sizes(eta);
    for (size_t k = 0; k < eta; k++) {
        if (k % 2 == 0) {
            int_sizes[k] = std::min(n, tau + 10);
        } else {
            int_sizes[k] = (tau >= 10) ? tau - 10 : 0;
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

// Paillier micro-benchmark
struct PaillierMicroBench {
    double scalar_mul_ms;
    double add_ms;
    double encrypt_ms;
    double decrypt_ms;
};

PaillierMicroBench run_paillier_microbench() {
    PaillierMicroBench result = {};
    const size_t NUM_OPS = 200;

    auto kp = Paillier::KeyGen(1024);

    // Prepare test data
    std::vector<BigInt> ciphertexts(NUM_OPS);
    std::vector<BigInt> scalars(NUM_OPS);
    std::vector<BigInt> plaintexts(NUM_OPS);

    for (size_t i = 0; i < NUM_OPS; i++) {
        plaintexts[i] = RandomBigInt(kp.pk.N);
        ciphertexts[i] = Paillier::Encrypt(kp.pk, plaintexts[i]);
        scalars[i] = RandomBigInt(kp.pk.N);
    }

    // ScalarMul
    {
        auto t0 = high_resolution_clock::now();
        for (size_t i = 0; i < NUM_OPS; i++) {
            Paillier::ScalarMul(kp.pk, ciphertexts[i], scalars[i]);
        }
        auto t1 = high_resolution_clock::now();
        result.scalar_mul_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0 / NUM_OPS;
    }

    // Add
    {
        auto t0 = high_resolution_clock::now();
        for (size_t i = 0; i < NUM_OPS; i++) {
            Paillier::Add(kp.pk, ciphertexts[i], ciphertexts[(i + 1) % NUM_OPS]);
        }
        auto t1 = high_resolution_clock::now();
        result.add_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0 / NUM_OPS;
    }

    // Encrypt
    {
        auto t0 = high_resolution_clock::now();
        for (size_t i = 0; i < NUM_OPS; i++) {
            Paillier::Encrypt(kp.pk, plaintexts[i]);
        }
        auto t1 = high_resolution_clock::now();
        result.encrypt_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0 / NUM_OPS;
    }

    // Decrypt
    {
        auto t0 = high_resolution_clock::now();
        for (size_t i = 0; i < NUM_OPS; i++) {
            Paillier::Decrypt(kp.sk, kp.pk, ciphertexts[i]);
        }
        auto t1 = high_resolution_clock::now();
        result.decrypt_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0 / NUM_OPS;
    }

    return result;
}

struct GS19BenchmarkResult {
    size_t n, tau, t, m;
    size_t eta;
    double alice_offline_ms;
    double bob_subtraction_ms;  // Average per sender
    double finv_per_iter_ms;
    double finv_total_ms;       // Measured or extrapolated
    double total_per_sender_ms;
    bool extrapolated;
    size_t correct;
};

GS19BenchmarkResult run_gs19_benchmark(size_t n, size_t eta, size_t tau,
                                        size_t max_finv_iters = 0) {
    GS19BenchmarkResult result = {};
    result.n = n;
    result.tau = tau;
    result.eta = eta;

    // Initialize parameters
    GS19Params params = GS19PICT::InitParams(n, tau);
    result.t = params.t;
    result.m = params.m;

    size_t beta = ComputeBeta(n);

    // Load or generate test data (same sets as otPSI-CA benchmark)
    ElementVec X = GetReceiverSet(n, beta);
    std::vector<ElementVec> Y_sets = GetSenderSets(X, n, eta, tau, beta);
    auto actual_sizes = ComputeIntersections(X, Y_sets);

    // --- Alice offline ---
    auto t0 = high_resolution_clock::now();
    auto alice_state = GS19PICT::AliceOffline(X, params);
    auto t1 = high_resolution_clock::now();
    result.alice_offline_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0;

    // Determine if we should run full or partial
    bool run_full = (max_finv_iters == 0);
    if (run_full) {
        // Full protocol for all sender sets
        result.correct = 0;
        double total_bob_sub = 0;
        double total_finv = 0;

        for (size_t k = 0; k < eta; k++) {
            // Measure Bob's subtraction (power sums + homomorphic sub)
            auto tb0 = high_resolution_clock::now();
            // We measure the full RunForSender which includes subtraction + FINV
            bool gs19_result = GS19PICT::RunForSender(alice_state, Y_sets[k]);
            auto tb1 = high_resolution_clock::now();
            double sender_ms = duration_cast<microseconds>(tb1 - tb0).count() / 1000.0;
            total_finv += sender_ms;

            bool expected = (actual_sizes[k] >= tau);
            if (gs19_result == expected) result.correct++;

            std::cout << "  Sender " << k << "/" << eta
                      << ": |X∩Y|=" << actual_sizes[k]
                      << " result=" << (gs19_result ? "similar" : "different")
                      << " expected=" << (expected ? "similar" : "different")
                      << " (" << std::fixed << std::setprecision(1) << sender_ms << " ms)"
                      << std::endl;
        }

        result.finv_total_ms = total_finv / eta;  // Average per sender
        result.total_per_sender_ms = result.finv_total_ms;
        result.finv_per_iter_ms = result.finv_total_ms / (2 * params.m);
        result.extrapolated = false;
    } else {
        // Partial: measure a few FINV iterations and extrapolate
        double per_iter_ms = 0;

        // Use first sender set for measurement
        GS19PICT::RunForSenderPartial(alice_state, Y_sets[0], max_finv_iters, per_iter_ms);

        result.finv_per_iter_ms = per_iter_ms;
        result.finv_total_ms = per_iter_ms * (2 * params.m);
        result.total_per_sender_ms = result.finv_total_ms;
        result.extrapolated = true;

        // Verify correctness with plaintext test
        result.correct = 0;
        for (size_t k = 0; k < eta; k++) {
            bool pt_result = GS19PICT::PlaintextTest(X, Y_sets[k], params);
            bool expected = (actual_sizes[k] >= tau);
            if (pt_result == expected) result.correct++;
        }
    }

    return result;
}

void print_gs19_header() {
    std::cout << std::setw(6) << "n"
              << std::setw(6) << "tau"
              << std::setw(6) << "t"
              << std::setw(6) << "m"
              << " | " << std::setw(12) << "Alice Off"
              << " | " << std::setw(12) << "FINV/iter"
              << " | " << std::setw(14) << "FINV total"
              << " | " << std::setw(14) << "Per-sender"
              << " | " << std::setw(10) << "Correct"
              << std::endl;
    std::cout << std::string(100, '-') << std::endl;
}

void print_gs19_result(const GS19BenchmarkResult& r) {
    std::cout << std::setw(6) << r.n
              << std::setw(6) << r.tau
              << std::setw(6) << r.t
              << std::setw(6) << r.m;

    auto format_time = [](double ms, bool extrap) -> std::string {
        std::ostringstream oss;
        if (ms < 1000) {
            oss << std::fixed << std::setprecision(1) << ms << " ms";
        } else if (ms < 60000) {
            oss << std::fixed << std::setprecision(1) << ms / 1000.0 << " s";
        } else if (ms < 3600000) {
            oss << std::fixed << std::setprecision(1) << ms / 60000.0 << " min";
        } else {
            oss << std::fixed << std::setprecision(1) << ms / 3600000.0 << " hrs";
        }
        if (extrap) oss << "*";
        return oss.str();
    };

    std::cout << " | " << std::setw(10) << format_time(r.alice_offline_ms, false)
              << " | " << std::setw(10) << format_time(r.finv_per_iter_ms, false)
              << " | " << std::setw(12) << format_time(r.finv_total_ms, r.extrapolated)
              << " | " << std::setw(12) << format_time(r.total_per_sender_ms, r.extrapolated)
              << " | " << std::setw(4) << r.correct << "/" << r.eta
              << (r.correct == r.eta ? " ok" : " FAIL")
              << (r.extrapolated ? "(pt)" : "")
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== GS19 PICT Benchmark (Baseline) ===" << std::endl;
    std::cout << "Ghosh & Simkin, CRYPTO 2019, Section 5" << std::endl;
    std::cout << "Paillier key: 1024-bit, kappa=40" << std::endl;
    std::cout << std::endl;

    // Paillier micro-benchmarks
    std::cout << "--- Paillier Micro-benchmarks (1024-bit, 200 ops) ---" << std::endl;
    auto micro = run_paillier_microbench();
    std::cout << "  ScalarMul: " << std::fixed << std::setprecision(3)
              << micro.scalar_mul_ms << " ms/op" << std::endl;
    std::cout << "  Add:       " << micro.add_ms << " ms/op" << std::endl;
    std::cout << "  Encrypt:   " << micro.encrypt_ms << " ms/op" << std::endl;
    std::cout << "  Decrypt:   " << micro.decrypt_ms << " ms/op" << std::endl;
    std::cout << std::endl;

    // Configurable n values
    std::vector<size_t> n_values = {100};
    if (argc > 1) {
        n_values.clear();
        for (int i = 1; i < argc; i++) {
            n_values.push_back(std::stoull(argv[i]));
        }
    }

    size_t eta = 10;
    std::vector<size_t> tau_ratios = {50, 60};  // tau = n * ratio / 100

    // Threshold for full vs partial FINV
    // For m > MAX_FULL_M, use partial measurement (3 iterations + extrapolate)
    // m=101 with 10 senders would take ~5 hours at 0.93ms/ScalarMul, so we
    // extrapolate for any m > 20 and verify correctness via plaintext test.
    const size_t MAX_FULL_M = 20;

    for (size_t n : n_values) {
        std::cout << "=== n = " << n << " ===" << std::endl;
        print_gs19_header();

        for (size_t ratio : tau_ratios) {
            size_t tau = n * ratio / 100;
            GS19Params params = GS19PICT::InitParams(n, tau);

            size_t max_iters = 0;  // 0 = full run
            if (params.m > MAX_FULL_M) {
                max_iters = 3;
                std::cout << "  (m=" << params.m
                          << " > " << MAX_FULL_M
                          << ": using " << max_iters
                          << " FINV iterations, extrapolating)" << std::endl;
            }

            auto result = run_gs19_benchmark(n, eta, tau, max_iters);
            print_gs19_result(result);
            std::cout << std::endl;
        }
    }

    std::cout << "* = extrapolated from per-iteration measurement" << std::endl;
    std::cout << "(pt) = correctness verified via plaintext test only" << std::endl;

    return 0;
}
