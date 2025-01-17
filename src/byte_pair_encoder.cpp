//
// Created by tiho on 12/6/22.
//

#include "../include/byte_pair_encoder.h"

std::vector<Kmer> BytePairEncoder::Encode(const std::unique_ptr<NucleicAcid> &sequence, std::uint32_t kmerLength) {
    std::vector<Kmer> dst;

    for (std::uint32_t i = 0; i < sequence->inflated_len; ++i) {
        std::uint64_t c = sequence->Code(i);
        char ch = kNucleotideDecoder[c];

        if (i <= sequence->inflated_len - kmerLength)
            dst.emplace_back(i + 1);

        for (std::int32_t k = 0; k < kmerLength; ++k)
            if (i >= k)
                dst[i - k].nucleotides.emplace_back(ch);
    }

    int leftch, rightch, code, oldsize;
    int index, r, w, best, done = 0;

    // Compress each data block until end of file
    while (!done) {
        //done = Fileread(infile);
        code = 256;

        for (;;) {
            // Get next unused char for pair code
            for (code--; code >= 0; code--)
                if (code == leftcode[code] && !rightcode[code])
                    break;

            // Must quit if no unused chars left
            if (code < 0)
                break;

            // Find most frequent pair of chars
            for (best = 2, index = 0; index < HASHSIZE; index++)
                if (count[index] > best) {
                    best = count[index];
                    leftch = left[index];
                    rightch = right[index];
                }

            // Done if no more compression possible
            if (best < THRESHOLD)
                break;

            // Replace pairs in data, adjust pair counts
            oldsize = size - 1;

            for (w = 0, r = 0; r < oldsize; r++) {
                if (buffer[r] == leftch && buffer[r + 1] == rightch) {
                    if (r > 0) {
                        index = Lookup(buffer[w - 1], leftch);
                        if (count[index] > 1) --count[index];
                        index = Lookup(buffer[w - 1], code);
                        if (count[index] < 255) ++count[index];
                    }

                    if (r < oldsize - 1) {
                        index = Lookup(rightch, buffer[r + 2]);
                        if (count[index] > 1) --count[index];
                        index = Lookup(code, buffer[r + 2]);
                        if (count[index] < 255) ++count[index];
                    }
                    buffer[w++] = code;
                    r++;
                    size--;
                } else {
                    buffer[w++] = buffer[r];
                }
            }

            buffer[w] = buffer[r];

            // Add to pair substitution table
            leftcode[code] = leftch;
            rightcode[code] = rightch;

            // Delete pair from hash table
            index = Lookup(leftch, rightch);
            count[index] = 1;
        }
    }

    return dst;
}

// Return index of character pair in hash table
// Deleted nodes have count of 1 for hashing
int BytePairEncoder::Lookup(unsigned char a, unsigned char b) {
    int index;
    // Compute hash key from both characters
    index = (a ^ (b << 5)) & (HASHSIZE - 1);

    // Search for pair or first empty slot
    while ((left[index] != a || right[index] != b) && count[index] != 0)
        index = (index + 1) & (HASHSIZE - 1);

    // Store pair in table
    left[index] = a;
    right[index] = b;

    return index;
}

// Read next block from input file into buffer
int BytePairEncoder::Fileread(FILE *input) {
    int c, index, used = 0;

    // Reset hash table and pair table
    for (c = 0; c < HASHSIZE; c++)
        count[c] = 0;

    for (c = 0; c < 256; c++) {
        leftcode[c] = c;
        rightcode[c] = 0;
    }

    size = 0;

    // Read data until full or few unused chars
    while (size < BLOCKSIZE && used < MAXCHARS && (c = getc(input)) != EOF) {
        if (size > 0) {
            index = Lookup(buffer[size - 1], c);
            if (count[index] < 255) ++count[index];
        }
        buffer[size++] = c;

        // Use rightcode to flag data chars found
        if (!rightcode[c]) {
            rightcode[c] = 1;
            used++;
        }
    }

    return c == EOF;
}

// Write each pair table and data block to output
void BytePairEncoder::Filewrite(FILE *output) {
    int i, len, c = 0;

    // For each character 0..255
    while (c < 256) {
        if (c == leftcode[c]) {
            // If not a pair code, count run of literals
            len = 1;
            c++;
            while (len < 127 && c < 256 && c == leftcode[c]) {
                len++;
                c++;
            }

            putc(len + 127, output);
            len = 0;
            if (c == 256)
                break;
        } else {
            // Else count run of pair codes
            len = 0;
            c++;
            while (len < 127 && c < 256 && c != leftcode[c] ||
                   len < 125 && c < 254 && c + 1 != leftcode[c + 1]) {
                len++;
                c++;
            }

            putc(len, output);
            c -= len + 1;
        }

        // Write range of pairs to output
        for (i = 0; i <= len; i++) {
            putc(leftcode[c], output);

            if (c != leftcode[c])
                putc(rightcode[c], output);

            c++;
        }
    }

    // Write size bytes and compressed data block
    putc(size / 256, output);
    putc(size % 256, output);
    fwrite(buffer, size, 1, output);
}


// Compress from input file to output file
void BytePairEncoder::Compress(FILE *infile, FILE *outfile) {
    int leftch, rightch, code, oldsize;
    int index, r, w, best, done = 0;

    // Compress each data block until end of file
    while (!done) {
        done = Fileread(infile);
        code = 256;

        for (;;) {
            // Get next unused char for pair code
            for (code--; code >= 0; code--)
                if (code == leftcode[code] && !rightcode[code])
                    break;

            // Must quit if no unused chars left
            if (code < 0)
                break;

            // Find most frequent pair of chars
            for (best = 2, index = 0; index < HASHSIZE; index++)
                if (count[index] > best) {
                    best = count[index];
                    leftch = left[index];
                    rightch = right[index];
                }

            // Done if no more compression possible
            if (best < THRESHOLD)
                break;

            // Replace pairs in data, adjust pair counts
            oldsize = size - 1;

            for (w = 0, r = 0; r < oldsize; r++) {
                if (buffer[r] == leftch && buffer[r + 1] == rightch) {
                    if (r > 0) {
                        index = Lookup(buffer[w - 1], leftch);
                        if (count[index] > 1) --count[index];
                        index = Lookup(buffer[w - 1], code);
                        if (count[index] < 255) ++count[index];
                    }

                    if (r < oldsize - 1) {
                        index = Lookup(rightch, buffer[r + 2]);
                        if (count[index] > 1) --count[index];
                        index = Lookup(code, buffer[r + 2]);
                        if (count[index] < 255) ++count[index];
                    }
                    buffer[w++] = code;
                    r++;
                    size--;
                } else {
                    buffer[w++] = buffer[r];
                }
            }

            buffer[w] = buffer[r];

            // Add to pair substitution table
            leftcode[code] = leftch;
            rightcode[code] = rightch;

            // Delete pair from hash table
            index = Lookup(leftch, rightch);
            count[index] = 1;
        }

        Filewrite(outfile);
    }
}

// Decompress data from input to output
void BytePairEncoder::Expand(FILE *input, FILE *output) {
    unsigned char left[256], right[256], stack[30];
    short int c, count, i, size;

    // Unpack each block until end of file
    while ((count = getc(input)) != EOF) {
        // Set left to itself as literal flag
        for (i = 0; i < 256; i++)
            left[i] = i;

        // Read pair table
        for (c = 0;;) {
            // Skip range of literal bytes
            if (count > 127) {
                c += count - 127;
                count = 0;
            }
            if (c == 256) break;

            // Read pairs, skip right if literal
            for (i = 0; i <= count; i++, c++) {
                left[c] = getc(input);
                if (c != left[c])
                    right[c] = getc(input);
            }
            if (c == 256) break;
            count = getc(input);
        }

        // Calculate packed data block size
        size = 256 * getc(input) + getc(input);

        // Unpack data block
        for (i = 0;;) {
            // Pop byte from stack or read byte
            if (i)
                c = stack[--i];
            else {
                if (!size--) break;
                c = getc(input);
            }
            // Output byte or push pair on stack
            if (c == left[c])
                putc(c, output);
            else {
                stack[i++] = right[c];
                stack[i++] = left[c];
            }
        }
    }
}
