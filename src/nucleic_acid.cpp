//
// Created by tiho on 12/6/22.
//

#include "../include/nucleic_acid.h"

NucleicAcid::NucleicAcid(
        const std::string& name,
        const std::string& data)
        : NucleicAcid(
        name.c_str(), name.size(),
        data.c_str(), data.size()) {}

NucleicAcid::NucleicAcid(
        const char* name, std::uint32_t name_len,
        const char* data, std::uint32_t data_len)
        : id(num_objects++),
          name(name, name_len),
          deflated_data(),
          block_quality(),
          inflated_len(data_len),
          is_reverse_complement(0) {
    deflated_data.reserve(data_len / 32. + .999);
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < data_len; ++i) {
        std::uint64_t c = kNucleotideCoder[static_cast<std::uint8_t>(data[i])];
        if (c == 255ULL) {
            throw std::invalid_argument(
                    "[biosoup::NucleicAcid::NucleicAcid] error: not a nucleotide");
        }
        block |= c << ((i << 1) & 63);
        if (((i + 1) & 31) == 0 || i == data_len - 1) {
            deflated_data.emplace_back(block);
            block = 0;
        }
    }
}

NucleicAcid::NucleicAcid(
        const char* name, std::uint32_t name_len,
        const char* data, std::uint32_t data_len,
        const char* quality, std::uint32_t quality_len)
        : NucleicAcid(
        name, name_len,
        data, data_len) {
    block_quality.reserve(quality_len / 64. + .999);
    for (std::uint32_t i = 0; i < quality_len; i += 64) {
        std::uint32_t j = std::min(i + 64, quality_len);
        auto block = std::accumulate(
                quality + i,
                quality + j,
                0,
                [] (const std::uint32_t& sum, const char& q) -> std::uint32_t {
                    return sum + (q - '!');
                });
        block_quality.emplace_back(block / (j - i));
    }
}

std::uint64_t NucleicAcid::Code(std::uint32_t i) const {
    std::uint64_t x = 0;
    if (is_reverse_complement) {
        i = inflated_len - i - 1;
        x = 3;
    }
    return ((deflated_data[i >> 5] >> ((i << 1) & 63)) & 3) ^ x;
}

std::uint8_t NucleicAcid::Score(std::uint32_t i) const {
    if (is_reverse_complement) {
        i = inflated_len - i - 1;
    }
    return block_quality[i >> 6];
}

std::string NucleicAcid::InflateData(std::uint32_t i, std::uint32_t len) const {
    if (i >= inflated_len) {
        return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
        dst += kNucleotideDecoder[Code(i)];
    }
    return dst;
}

std::string NucleicAcid::InflateQuality(std::uint32_t i, std::uint32_t len) const {  // NOLINT
    if (block_quality.empty() || i >= inflated_len) {
        return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
        dst += Score(i) + '!';
    }
    return dst;
}

void NucleicAcid::ReverseAndComplement() {   // Watson-Crick base pairing
    is_reverse_complement ^= 1;
}