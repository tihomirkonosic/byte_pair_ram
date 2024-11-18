//
// Created by tiho on 12/7/22.
//

#ifndef BYTEPAIRRAM_KMER_H
#define BYTEPAIRRAM_KMER_H

#include <cstdint>
#include <vector>

struct Kmer {
public:
    Kmer() = default;
    Kmer(std::uint64_t position) : position(position) {}

    std::uint64_t position;
    std::vector<char> nucleotides;
};

#endif //BYTEPAIRRAM_KMER_H
