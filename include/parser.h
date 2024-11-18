//
// Created by Tiho on 6.12.2022..
//

#ifndef BYTEPAIRRAM_PARSER_H
#define BYTEPAIRRAM_PARSER_H

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>
#include "nucleic_acid.h"

class Parser {
public:
    Parser(const Parser&) = delete;
    Parser& operator=(const Parser&) = delete;

    Parser(Parser&&) = delete;
    Parser& operator=(Parser&&) = delete;

    virtual ~Parser() {}

    static std::unique_ptr<Parser> Create(const std::string& path);

    void Reset();

    std::vector<std::unique_ptr<NucleicAcid>> Parse(std::uint64_t bytes, bool shorten_names = true);

protected:
    Parser(gzFile file, std::uint32_t storage_size);

    const std::vector<char>& buffer() const {
        return buffer_;
    }

    std::uint32_t buffer_ptr() const {
        return buffer_ptr_;
    }

    std::uint32_t buffer_bytes() const {
        return buffer_bytes_;
    }

    const std::vector<char>& storage() const {
        return storage_;
    }

    std::uint32_t storage_ptr() const {
        return storage_ptr_;
    }

    bool Read();

    void Store(std::uint32_t count, bool strip = false);

    void Terminate(std::uint32_t i);

    void Clear();

    static std::uint32_t RightStrip(const char* str, std::uint32_t str_len);

    static std::uint32_t Shorten(const char* str, std::uint32_t str_len);

private:
    std::unique_ptr<gzFile_s, int(*)(gzFile)> file_;
    std::vector<char> buffer_;
    std::uint32_t buffer_ptr_;
    std::uint32_t buffer_bytes_;
    std::vector<char> storage_;
    std::uint32_t storage_ptr_;
};

#endif //BYTEPAIRRAM_PARSER_H
