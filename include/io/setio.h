#ifndef OTPSICA_SETIO_H
#define OTPSICA_SETIO_H

#include "common.h"
#include <string>
#include <fstream>
#include <stdexcept>

namespace otpsica {

/**
 * File I/O for element sets.
 *
 * Supports two formats:
 * - Text format (.txt): Human-readable, one element per line
 * - Binary format (.bin): Compact binary representation
 *
 * Text format for single set:
 *   # Set: <size> elements, <beta> bits
 *   <element_1>
 *   <element_2>
 *   ...
 *
 * Text format for multiple sets:
 *   # Sets: <num_sets> sets, <size> elements each, <beta> bits
 *   # Set 0
 *   <element_1>
 *   ...
 *   # Set 1
 *   <element_1>
 *   ...
 *
 * Binary format:
 *   [4 bytes: magic "SETS"]
 *   [4 bytes: version]
 *   [8 bytes: num_sets]
 *   [8 bytes: elements_per_set]
 *   [8 bytes: beta]
 *   [data: num_sets * elements_per_set * 8 bytes]
 */
class SetIO {
public:
    // File format enumeration
    enum class Format {
        Text,
        Binary,
        Auto  // Detect from file extension
    };

    /**
     * Write a single set to file.
     * @param filename Output file path
     * @param elements The set to write
     * @param beta Bit-length of elements (for metadata)
     * @param format File format (default: auto-detect from extension)
     */
    static void WriteSet(const std::string& filename,
                         const ElementVec& elements,
                         size_t beta = 0,
                         Format format = Format::Auto);

    /**
     * Write multiple sets to file.
     * @param filename Output file path
     * @param sets Vector of sets to write
     * @param beta Bit-length of elements (for metadata)
     * @param format File format (default: auto-detect from extension)
     */
    static void WriteSets(const std::string& filename,
                          const std::vector<ElementVec>& sets,
                          size_t beta = 0,
                          Format format = Format::Auto);

    /**
     * Read a single set from file.
     * @param filename Input file path
     * @param format File format (default: auto-detect from extension)
     * @return The set of elements
     */
    static ElementVec ReadSet(const std::string& filename,
                              Format format = Format::Auto);

    /**
     * Read multiple sets from file.
     * @param filename Input file path
     * @param format File format (default: auto-detect from extension)
     * @return Vector of sets
     */
    static std::vector<ElementVec> ReadSets(const std::string& filename,
                                            Format format = Format::Auto);

    /**
     * Read set metadata without loading full data.
     * @param filename Input file path
     * @param[out] num_sets Number of sets in file
     * @param[out] set_size Number of elements per set
     * @param[out] beta Bit-length of elements
     */
    static void ReadMetadata(const std::string& filename,
                             size_t& num_sets,
                             size_t& set_size,
                             size_t& beta);

private:
    // Text format I/O
    static void WriteSetText(std::ofstream& out,
                             const ElementVec& elements,
                             size_t beta);

    static void WriteSetsText(std::ofstream& out,
                              const std::vector<ElementVec>& sets,
                              size_t beta);

    static ElementVec ReadSetText(std::ifstream& in);

    static std::vector<ElementVec> ReadSetsText(std::ifstream& in);

    // Binary format I/O
    static void WriteSetsBinary(std::ofstream& out,
                                const std::vector<ElementVec>& sets,
                                size_t beta);

    static std::vector<ElementVec> ReadSetsBinary(std::ifstream& in);

    // Utility functions
    static Format DetectFormat(const std::string& filename);
    static bool HasExtension(const std::string& filename, const std::string& ext);

    // Binary format constants
    static constexpr uint32_t MAGIC = 0x53455453;  // "SETS"
    static constexpr uint32_t VERSION = 1;
};

} // namespace otpsica

#endif // OTPSICA_SETIO_H
