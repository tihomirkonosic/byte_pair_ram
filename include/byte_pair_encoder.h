//
// Created by tiho on 12/6/22.
//

#ifndef BYTEPAIRRAM_BYTE_PAIR_ENCODER_H
#define BYTEPAIRRAM_BYTE_PAIR_ENCODER_H

#include <stdio.h>
#include <memory>
#include "nucleic_acid.h"
#include "Kmer.h"

#define BLOCKSIZE 5000 /* Maximum block size */
#define HASHSIZE 4096 /* Size of hash table */
#define MAXCHARS 200 /* Char set per block */
#define THRESHOLD 3 /* Minimum pair count */

class BytePairEncoder {
public:
    std::vector<Kmer> Encode(const std::unique_ptr<NucleicAcid> &sequence, std::uint32_t kmerLength);
    int Lookup (unsigned char a, unsigned char b);
    int Fileread (FILE *input);
    void Filewrite (FILE *);
    void Compress (FILE *, FILE *);
    void Expand (FILE *input, FILE *output);

private:
    unsigned char buffer[BLOCKSIZE]; /* Data block */
    unsigned char leftcode[256]; /* Pair table */
    unsigned char rightcode[256]; /* Pair table */
    unsigned char left[HASHSIZE]; /* Hash table */
    unsigned char right[HASHSIZE]; /* Hash table */
    unsigned char count[HASHSIZE]; /* Pair count */
    int size; /* Size of current data block */
};

#endif //BYTEPAIRRAM_BYTE_PAIR_ENCODER_H
