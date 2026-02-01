#include "io/setio.h"
#include <sstream>
#include <cstring>

namespace otpsica {

// Utility: check file extension
bool SetIO::HasExtension(const std::string& filename, const std::string& ext) {
    if (filename.length() < ext.length()) return false;
    return filename.compare(filename.length() - ext.length(), ext.length(), ext) == 0;
}

// Detect format from file extension
SetIO::Format SetIO::DetectFormat(const std::string& filename) {
    if (HasExtension(filename, ".bin") || HasExtension(filename, ".dat")) {
        return Format::Binary;
    }
    return Format::Text;
}

// Write single set (delegates to WriteSets)
void SetIO::WriteSet(const std::string& filename,
                     const ElementVec& elements,
                     size_t beta,
                     Format format) {
    std::vector<ElementVec> sets = {elements};
    WriteSets(filename, sets, beta, format);
}

// Write multiple sets
void SetIO::WriteSets(const std::string& filename,
                      const std::vector<ElementVec>& sets,
                      size_t beta,
                      Format format) {
    if (format == Format::Auto) {
        format = DetectFormat(filename);
    }

    std::ofstream out;
    if (format == Format::Binary) {
        out.open(filename, std::ios::binary);
    } else {
        out.open(filename);
    }

    if (!out.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }

    if (format == Format::Binary) {
        WriteSetsBinary(out, sets, beta);
    } else {
        WriteSetsText(out, sets, beta);
    }

    out.close();
}

// Read single set
ElementVec SetIO::ReadSet(const std::string& filename, Format format) {
    auto sets = ReadSets(filename, format);
    if (sets.empty()) {
        return ElementVec();
    }
    return sets[0];
}

// Read multiple sets
std::vector<ElementVec> SetIO::ReadSets(const std::string& filename, Format format) {
    if (format == Format::Auto) {
        format = DetectFormat(filename);
    }

    std::ifstream in;
    if (format == Format::Binary) {
        in.open(filename, std::ios::binary);
    } else {
        in.open(filename);
    }

    if (!in.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    std::vector<ElementVec> result;
    if (format == Format::Binary) {
        result = ReadSetsBinary(in);
    } else {
        result = ReadSetsText(in);
    }

    in.close();
    return result;
}

// Read metadata only
void SetIO::ReadMetadata(const std::string& filename,
                         size_t& num_sets,
                         size_t& set_size,
                         size_t& beta) {
    Format format = DetectFormat(filename);

    std::ifstream in;
    if (format == Format::Binary) {
        in.open(filename, std::ios::binary);
    } else {
        in.open(filename);
    }

    if (!in.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    if (format == Format::Binary) {
        // Read binary header
        uint32_t magic, version;
        uint64_t ns, ss, b;

        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        if (magic != MAGIC) {
            throw std::runtime_error("Invalid binary file format (bad magic)");
        }

        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&ns), sizeof(ns));
        in.read(reinterpret_cast<char*>(&ss), sizeof(ss));
        in.read(reinterpret_cast<char*>(&b), sizeof(b));

        num_sets = static_cast<size_t>(ns);
        set_size = static_cast<size_t>(ss);
        beta = static_cast<size_t>(b);
    } else {
        // Parse text header
        std::string line;
        num_sets = 0;
        set_size = 0;
        beta = 0;

        while (std::getline(in, line)) {
            if (line.empty()) continue;

            if (line[0] == '#') {
                // Parse header comments
                if (line.find("Sets:") != std::string::npos) {
                    // Multi-set header: "# Sets: <num_sets> sets, <size> elements each, <beta> bits"
                    std::sscanf(line.c_str(), "# Sets: %zu sets, %zu elements each, %zu bits",
                               &num_sets, &set_size, &beta);
                } else if (line.find("Set:") != std::string::npos && num_sets == 0) {
                    // Single set header: "# Set: <size> elements, <beta> bits"
                    std::sscanf(line.c_str(), "# Set: %zu elements, %zu bits",
                               &set_size, &beta);
                    num_sets = 1;
                }
            } else {
                // First data line reached, stop reading
                break;
            }
        }
    }

    in.close();
}

// Text format: write single set
void SetIO::WriteSetText(std::ofstream& out,
                         const ElementVec& elements,
                         size_t beta) {
    out << "# Set: " << elements.size() << " elements, " << beta << " bits\n";
    for (Element e : elements) {
        out << e << "\n";
    }
}

// Text format: write multiple sets
void SetIO::WriteSetsText(std::ofstream& out,
                          const std::vector<ElementVec>& sets,
                          size_t beta) {
    if (sets.empty()) {
        out << "# Sets: 0 sets, 0 elements each, " << beta << " bits\n";
        return;
    }

    size_t max_size = 0;
    for (const auto& s : sets) {
        if (s.size() > max_size) max_size = s.size();
    }

    out << "# Sets: " << sets.size() << " sets, " << max_size << " elements each, " << beta << " bits\n";

    for (size_t k = 0; k < sets.size(); k++) {
        out << "# Set " << k << "\n";
        for (Element e : sets[k]) {
            out << e << "\n";
        }
    }
}

// Text format: read sets
std::vector<ElementVec> SetIO::ReadSetsText(std::ifstream& in) {
    std::vector<ElementVec> result;
    ElementVec current_set;
    bool in_set = false;

    std::string line;
    while (std::getline(in, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Comment line
            if (line.find("# Set ") != std::string::npos || line.find("# Set:") != std::string::npos) {
                // Start of a new set
                if (in_set && !current_set.empty()) {
                    result.push_back(std::move(current_set));
                    current_set.clear();
                }
                in_set = true;
            }
            // Skip other comment lines
            continue;
        }

        // Data line: parse element
        Element e;
        std::istringstream iss(line);
        if (iss >> e) {
            current_set.push_back(e);
            in_set = true;
        }
    }

    // Add last set if not empty
    if (!current_set.empty()) {
        result.push_back(std::move(current_set));
    }

    return result;
}

// Text format: read single set (for backward compatibility)
ElementVec SetIO::ReadSetText(std::ifstream& in) {
    auto sets = ReadSetsText(in);
    if (sets.empty()) return ElementVec();
    return sets[0];
}

// Binary format: write sets
void SetIO::WriteSetsBinary(std::ofstream& out,
                            const std::vector<ElementVec>& sets,
                            size_t beta) {
    // Write header
    uint32_t magic = MAGIC;
    uint32_t version = VERSION;
    uint64_t num_sets = sets.size();

    // Find max set size (all sets padded to this size with 0)
    uint64_t max_size = 0;
    for (const auto& s : sets) {
        if (s.size() > max_size) max_size = s.size();
    }
    uint64_t beta64 = beta;

    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&num_sets), sizeof(num_sets));
    out.write(reinterpret_cast<const char*>(&max_size), sizeof(max_size));
    out.write(reinterpret_cast<const char*>(&beta64), sizeof(beta64));

    // Write actual sizes for each set
    for (const auto& s : sets) {
        uint64_t actual_size = s.size();
        out.write(reinterpret_cast<const char*>(&actual_size), sizeof(actual_size));
    }

    // Write data (each set padded to max_size)
    for (const auto& s : sets) {
        for (Element e : s) {
            out.write(reinterpret_cast<const char*>(&e), sizeof(e));
        }
        // Padding with zeros
        Element zero = 0;
        for (size_t i = s.size(); i < max_size; i++) {
            out.write(reinterpret_cast<const char*>(&zero), sizeof(zero));
        }
    }
}

// Binary format: read sets
std::vector<ElementVec> SetIO::ReadSetsBinary(std::ifstream& in) {
    // Read header
    uint32_t magic, version;
    uint64_t num_sets, max_size, beta;

    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != MAGIC) {
        throw std::runtime_error("Invalid binary file format (bad magic)");
    }

    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != VERSION) {
        throw std::runtime_error("Unsupported binary file version");
    }

    in.read(reinterpret_cast<char*>(&num_sets), sizeof(num_sets));
    in.read(reinterpret_cast<char*>(&max_size), sizeof(max_size));
    in.read(reinterpret_cast<char*>(&beta), sizeof(beta));
    (void)beta;  // Not used for reading, just metadata

    // Read actual sizes
    std::vector<uint64_t> actual_sizes(num_sets);
    for (size_t k = 0; k < num_sets; k++) {
        in.read(reinterpret_cast<char*>(&actual_sizes[k]), sizeof(uint64_t));
    }

    // Read data
    std::vector<ElementVec> result(num_sets);
    for (size_t k = 0; k < num_sets; k++) {
        result[k].resize(actual_sizes[k]);

        // Read actual elements
        for (size_t i = 0; i < actual_sizes[k]; i++) {
            in.read(reinterpret_cast<char*>(&result[k][i]), sizeof(Element));
        }

        // Skip padding
        Element dummy;
        for (size_t i = actual_sizes[k]; i < max_size; i++) {
            in.read(reinterpret_cast<char*>(&dummy), sizeof(dummy));
        }
    }

    return result;
}

} // namespace otpsica
