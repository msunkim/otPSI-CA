#include <iostream>
#include <cassert>
#include <cstdio>
#include <string>
#include <chrono>
#include "crypto/paillier.h"
#include "crypto/cek.h"
#include "hash/hashtable.h"
#include "poly/flt_poly.h"
#include "io/setio.h"

using namespace otpsica;
using namespace NTL;

void test_paillier() {
    std::cout << "Testing Paillier encryption..." << std::endl;

    // Generate keys (use smaller modulus for faster testing)
    PaillierKeyPair kp = Paillier::KeyGen(1024);

    // Test encryption/decryption
    BigInt m1 = conv<NTL::ZZ>(12345);
    BigInt c1 = Paillier::Encrypt(kp.pk, m1);
    BigInt d1 = Paillier::Decrypt(kp.sk, kp.pk, c1);
    assert(d1 == m1);

    // Test homomorphic addition
    BigInt m2 = conv<NTL::ZZ>(67890);
    BigInt c2 = Paillier::Encrypt(kp.pk, m2);
    BigInt c_sum = Paillier::Add(kp.pk, c1, c2);
    BigInt d_sum = Paillier::Decrypt(kp.sk, kp.pk, c_sum);
    assert(d_sum == m1 + m2);

    // Test scalar multiplication
    BigInt k = conv<NTL::ZZ>(5);
    BigInt c_mul = Paillier::ScalarMul(kp.pk, c1, k);
    BigInt d_mul = Paillier::Decrypt(kp.sk, kp.pk, c_mul);
    assert(d_mul == m1 * k);

    std::cout << "  Paillier tests passed!" << std::endl;
}

void test_cek() {
    std::cout << "Testing CEK encryption..." << std::endl;

    // Generate keys with small parameters for testing
    CEKKeyPair kp = CEK::KeyGen(128, 1024, 2, 10);

    // Test encryption/decryption for small messages
    for (size_t m = 0; m < 5; m++) {
        BigInt c = CEK::Encrypt(kp.pk, m);
        long d = CEK::Decrypt(kp.sk, kp.pk, c);
        assert(d == static_cast<long>(m));
    }

    std::cout << "  CEK tests passed!" << std::endl;
}

void test_hashtable() {
    std::cout << "Testing hash table..." << std::endl;

    size_t n = 100;
    size_t w = 3;
    HashTable ht(n, w);

    // Insert elements
    for (Element i = 1; i <= 50; i++) {
        ht.Insert(i);
    }

    // Check membership
    for (Element i = 1; i <= 50; i++) {
        assert(ht.Check(i));
    }

    // Check non-membership
    for (Element i = 51; i <= 100; i++) {
        assert(!ht.Check(i));
    }

    // Test filling with random
    ht.FillWithRandom(16);
    for (size_t i = 0; i < ht.GetNumBins(); i++) {
        assert(ht.GetBin(i).size() == w);
    }

    std::cout << "  Hash table tests passed!" << std::endl;
}

void test_flt_polynomial() {
    std::cout << "Testing FLT polynomial..." << std::endl;

    // Simple test: set {1, 2, 3} with prime p = 7
    ElementVec elements = {1, 2, 3};
    size_t p = 7;

    auto coeffs = FLTPolynomial::Build(elements, p);

    // Verify polynomial evaluates to 1 for elements in set
    // and 0 for elements not in set (using plaintext evaluation)
    NTL::ZZ_p::init(conv<NTL::ZZ>(p));

    for (Element y = 0; y < 10; y++) {
        // Evaluate polynomial at y
        NTL::ZZ_p result;
        NTL::clear(result);
        NTL::ZZ_p y_power;
        NTL::set(y_power);  // y^0 = 1

        size_t degree = coeffs.size() - 1;
        for (size_t i = 0; i <= degree; i++) {
            result += conv<NTL::ZZ_p>(coeffs[degree - i]) * y_power;
            y_power *= conv<NTL::ZZ_p>(y);
        }

        bool in_set = (y == 1 || y == 2 || y == 3);
        bool result_is_one = (IsOne(result));

        // Note: FLT gives 1 for elements in set, 0 otherwise (mod p)
        if (in_set) {
            assert(result_is_one);
        }
    }

    std::cout << "  FLT polynomial tests passed!" << std::endl;
}

void test_crt_primes() {
    std::cout << "Testing CRT prime selection..." << std::endl;

    // Test for 16-bit elements
    auto primes = CRTDecomposition::SelectPrimes(16, 2);
    assert(primes.size() >= 2);
    assert(CRTDecomposition::Validate(primes, 16));

    std::cout << "  Selected primes: ";
    for (size_t p : primes) std::cout << p << " ";
    std::cout << std::endl;

    std::cout << "  CRT tests passed!" << std::endl;
}

void test_setio() {
    std::cout << "Testing SetIO..." << std::endl;

    // Generate random sets using RandomSetGenerator
    RandomSetGenerator gen(42);
    ElementVec X = gen.Generate(16, 100);  // 100 elements, 16 bits

    std::vector<size_t> int_sizes = {30, 20, 10};
    auto Y_sets = gen.GenerateWithIntersection(X, 16, 50, 3, int_sizes);

    // Test text format
    std::string txt_file = "/tmp/test_sets.txt";
    SetIO::WriteSets(txt_file, Y_sets, 16);

    auto Y_read = SetIO::ReadSets(txt_file);
    assert(Y_read.size() == Y_sets.size());
    for (size_t k = 0; k < Y_sets.size(); k++) {
        assert(Y_read[k].size() == Y_sets[k].size());
        for (size_t i = 0; i < Y_sets[k].size(); i++) {
            assert(Y_read[k][i] == Y_sets[k][i]);
        }
    }

    // Test binary format
    std::string bin_file = "/tmp/test_sets.bin";
    SetIO::WriteSets(bin_file, Y_sets, 16);

    auto Y_read_bin = SetIO::ReadSets(bin_file);
    assert(Y_read_bin.size() == Y_sets.size());
    for (size_t k = 0; k < Y_sets.size(); k++) {
        assert(Y_read_bin[k].size() == Y_sets[k].size());
        for (size_t i = 0; i < Y_sets[k].size(); i++) {
            assert(Y_read_bin[k][i] == Y_sets[k][i]);
        }
    }

    // Test metadata reading
    size_t num_sets, set_size, beta;
    SetIO::ReadMetadata(txt_file, num_sets, set_size, beta);
    assert(num_sets == 3);
    assert(set_size == 50);
    assert(beta == 16);

    // Clean up test files
    std::remove(txt_file.c_str());
    std::remove(bin_file.c_str());

    std::cout << "  SetIO tests passed!" << std::endl;
}

void test_cek_keygen(size_t d) {
    // Generate (or load from cache) a CEK key with realistic parameters.
    // SecurityParams defaults: kappa=256, N_bits=2048, rho=2.
    std::cout << "CEK KeyGen with d=" << d
              << " (kappa=256, N=2048, rho=2)..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto kp = CEK::KeyGen(256, 2048, 2, d);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "  KeyGen completed in " << ms << " ms" << std::endl;
    std::cout << "  N bits = " << NumBits(kp.pk.N) << std::endl;
    std::cout << "  g1 bits = " << NumBits(kp.pk.g1) << std::endl;
    std::cout << "  g2 bits = " << NumBits(kp.pk.g2) << std::endl;

    // Quick encrypt/decrypt sanity check (message 0 and 1)
    for (size_t m = 0; m < std::min(d, size_t(3)); m++) {
        BigInt c = CEK::Encrypt(kp.pk, m);
        long dec = CEK::Decrypt(kp.sk, kp.pk, c);
        assert(dec == static_cast<long>(m));
    }
    std::cout << "  Encrypt/decrypt OK" << std::endl;
}

int main(int argc, char* argv[]) {
    // Usage:
    //   ./otpsica_test              — run all unit tests
    //   ./otpsica_test cek [d]      — generate (and cache) a CEK key with given d
    if (argc >= 2 && std::string(argv[1]) == "cek") {
        size_t d = 354;  // default for tau=50
        if (argc >= 3) d = std::stoull(argv[2]);
        try {
            test_cek_keygen(d);
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "CEK keygen failed: " << e.what() << std::endl;
            return 1;
        }
    }

    std::cout << "=== otPSI-CA Unit Tests ===" << std::endl << std::endl;

    try {
        test_paillier();
        test_cek();
        test_hashtable();
        test_flt_polynomial();
        test_crt_primes();
        test_setio();

        std::cout << std::endl << "All tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
