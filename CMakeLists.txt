cmake_minimum_required(VERSION 3.16)
project(otPSICA VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Find required libraries

# NTL library (for big integer and polynomial arithmetic)
find_library(NTL_LIB ntl REQUIRED)
find_path(NTL_INCLUDE NTL/ZZ.h REQUIRED)

# GMP library (required by NTL)
find_library(GMP_LIB gmp REQUIRED)

# Crypto++ library (for hashing)
# Try pkg-config first, then fall back to find_library
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(CRYPTOPP QUIET libcryptopp)
endif()

if(NOT CRYPTOPP_FOUND)
    # Try finding Crypto++ directly
    find_library(CRYPTOPP_LIBRARIES NAMES cryptopp crypto++ cryptlib
                 PATHS /usr/lib /usr/local/lib /opt/homebrew/lib)
    find_path(CRYPTOPP_INCLUDE_DIRS cryptopp/sha.h
              PATHS /usr/include /usr/local/include /opt/homebrew/include)

    if(CRYPTOPP_LIBRARIES AND CRYPTOPP_INCLUDE_DIRS)
        set(CRYPTOPP_FOUND TRUE)
    endif()
endif()

if(NOT CRYPTOPP_FOUND)
    message(WARNING "Crypto++ not found. Using OpenSSL for hashing instead.")
    find_package(OpenSSL REQUIRED)
    set(USE_OPENSSL TRUE)
    add_definitions(-DUSE_OPENSSL)
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${NTL_INCLUDE}
)

if(CRYPTOPP_FOUND)
    include_directories(${CRYPTOPP_INCLUDE_DIRS})
elseif(USE_OPENSSL)
    include_directories(${OPENSSL_INCLUDE_DIR})
endif()

# Source files
set(CRYPTO_SOURCES
    src/crypto/paillier.cc
    src/crypto/cek.cc
)

set(HASH_SOURCES
    src/hash/hashtable.cc
)

set(POLY_SOURCES
    src/poly/flt_poly.cc
)

set(PROTOCOL_SOURCES
    src/protocol/smt.cc
    src/protocol/gte.cc
    src/protocol/otpsica.cc
    src/protocol/gs19.cc
)

set(IO_SOURCES
    src/io/setio.cc
)

# Library
add_library(otpsica STATIC
    ${CRYPTO_SOURCES}
    ${HASH_SOURCES}
    ${POLY_SOURCES}
    ${PROTOCOL_SOURCES}
    ${IO_SOURCES}
)

target_link_libraries(otpsica
    ${NTL_LIB}
    ${GMP_LIB}
)

if(CRYPTOPP_FOUND)
    target_link_libraries(otpsica ${CRYPTOPP_LIBRARIES})
elseif(USE_OPENSSL)
    target_link_libraries(otpsica OpenSSL::Crypto)
endif()

# Main executable
add_executable(otpsica_main src/main.cc)
target_link_libraries(otpsica_main otpsica)

# Test executable
add_executable(otpsica_test tests/test_main.cc)
target_link_libraries(otpsica_test otpsica)

# Benchmark executable
add_executable(otpsica_bench tests/benchmark.cc)
target_link_libraries(otpsica_bench otpsica)

# GS19 baseline benchmark executable
add_executable(gs19_bench tests/benchmark_gs19.cc)
target_link_libraries(gs19_bench otpsica)

# Install
install(TARGETS otpsica otpsica_main
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include)
