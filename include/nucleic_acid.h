//
// Created by tiho on 12/6/22.
//

#ifndef BYTEPAIRRAM_NUCLEIC_ACID_H
#define BYTEPAIRRAM_NUCLEIC_ACID_H

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <numeric>
#include <string>
#include <stdexcept>
#include <vector>

constexpr static std::uint8_t kNucleotideCoder[] = {
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255,   0, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255,   0,   1 ,  1,   0, 255, 255,   2,
        3, 255, 255,   2, 255,   1,   0, 255,
        255, 255,   0,   1,   3,   3,   2,   0,
        255,   3, 255, 255, 255, 255, 255, 255,
        255,   0,   1,   1,   0, 255, 255,   2,
        3, 255, 255,   2, 255,   1,   0, 255,
        255, 255,   0,   1,   3,   3,   2,   0,
        255,   3, 255, 255, 255, 255, 255, 255
};

constexpr static char kNucleotideDecoder[] = {
        'A', 'C', 'G', 'T'
};

class NucleicAcid {
public:
    NucleicAcid() = default;

    NucleicAcid(
            const std::string& name,
            const std::string& data);

    NucleicAcid(
            const char* name, std::uint32_t name_len,
            const char* data, std::uint32_t data_len);

    NucleicAcid(
            const std::string& name,
            const std::string& data,
            const std::string& quality);

    NucleicAcid(
            const char* name, std::uint32_t name_len,
            const char* data, std::uint32_t data_len,
            const char* quality, std::uint32_t quality_len);

    NucleicAcid(const NucleicAcid&) = default;
    NucleicAcid& operator=(const NucleicAcid&) = default;

    NucleicAcid(NucleicAcid&&) = default;
    NucleicAcid& operator=(NucleicAcid&&) = default;

    ~NucleicAcid() = default;

    std::uint64_t Code(std::uint32_t i) const;

    std::uint8_t Score(std::uint32_t i) const;

    std::string InflateData(std::uint32_t i = 0, std::uint32_t len = -1) const;

    std::string InflateQuality(std::uint32_t i = 0, std::uint32_t len = -1) const;

    void ReverseAndComplement();

    static std::atomic<std::uint32_t> num_objects;
    //static std::uint32_t num_objects;

    std::uint32_t id;  // (optional) initialize num_objects to 0
    std::string name;
    std::vector<std::uint64_t> deflated_data;
    std::vector<std::uint8_t> block_quality;  // (optional) Phred quality scores
    std::uint32_t inflated_len;
    bool is_reverse_complement;
};

#endif //BYTEPAIRRAM_NUCLEIC_ACID_H
