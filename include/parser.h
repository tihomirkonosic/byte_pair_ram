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

    void Reset() {
        gzseek(file_.get(), 0, SEEK_SET);
        buffer_ptr_ = 0;
        buffer_bytes_ = 0;
    }

    std::vector<std::unique_ptr<NucleicAcid>> Parse(std::uint64_t bytes, bool shorten_names = true) {
        std::vector<std::unique_ptr<NucleicAcid>> dst;
        std::uint64_t parsed_bytes = 0;
        std::uint32_t data_ptr = 0;

        auto create_T = [&] () -> void {
            if (data_ptr == 0) {
                throw std::invalid_argument(
                        "[bioparser::FastaParser] error: invalid file format");
            }

            auto name_len = shorten_names ?
                            this->Shorten(this->storage().data(), data_ptr) :
                            this->RightStrip(this->storage().data(), data_ptr);

            auto data_len = this->storage_ptr() - data_ptr;

            if (name_len == 0 || this->storage()[0] != '>' || data_len == 0) {
                throw std::invalid_argument(
                        "[bioparser::FastaParser] error: invalid file format");
            }

            dst.emplace_back(std::unique_ptr<NucleicAcid>(new NucleicAcid(
                    static_cast<const char*>(this->storage().data() + 1), name_len - 1,
                    static_cast<const char*>(this->storage().data() + data_ptr), data_len)));  // NOLINT

            parsed_bytes += this->storage_ptr();
            data_ptr = 0;
            this->Clear();
        };

        bool is_eof = false;
        bool is_name = true;

        while (true) {
            auto buffer_ptr = this->buffer_ptr();
            for (; buffer_ptr < this->buffer_bytes(); ++buffer_ptr) {
                auto c = this->buffer()[buffer_ptr];
                if (c == '\n') {
                    this->Store(buffer_ptr - this->buffer_ptr(), !is_name);
                    if (is_name) {
                        data_ptr = this->storage_ptr();
                        is_name = false;
                    }
                } else if (!is_name && c == '>') {
                    is_name = true;
                    create_T();
                    if (parsed_bytes >= bytes) {
                        return dst;
                    }
                }
            }
            if (this->buffer_ptr() < buffer_ptr) {
                this->Store(buffer_ptr - this->buffer_ptr(), !is_name);
            }

            if (is_eof) {
                break;
            }
            is_eof = this->Read();
        }

        if (this->storage_ptr() != 0) {
            create_T();
        }

        return dst;
    }

protected:
    Parser(gzFile file, std::uint32_t storage_size)
            : file_(file, gzclose),
              buffer_(65536, 0),  // 64 kB
              buffer_ptr_(0),
              buffer_bytes_(0),
              storage_(storage_size, 0),
              storage_ptr_(0) {}

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

    bool Read() {
        buffer_ptr_ = 0;
        buffer_bytes_ = gzread(file_.get(), buffer_.data(), buffer_.size());
        return buffer_bytes_ < buffer_.size();
    }

    void Store(std::uint32_t count, bool strip = false) {
        if (buffer_ptr_ + count > buffer_.size()) {
            throw std::invalid_argument(
                    "[bioparser::Parser::Store] error: buffer overflow");
        }
        if (storage_ptr_ + count > storage_.size()) {
            storage_.resize(2 * storage_.size());
        }
        std::memcpy(&storage_[storage_ptr_], &buffer_[buffer_ptr_], count);
        storage_ptr_ += strip ? RightStrip(&storage_[storage_ptr_], count) : count;
        buffer_ptr_ += count + 1;  // ignore sought character
    }

    void Terminate(std::uint32_t i) {
        storage_[i] = '\0';
    }

    void Clear() {
        storage_ptr_ = 0;
    }

    static std::uint32_t RightStrip(const char* str, std::uint32_t str_len) {
        while (str_len > 0 && std::isspace(str[str_len - 1])) {
            --str_len;
        }
        return str_len;
    }

    static std::uint32_t Shorten(const char* str, std::uint32_t str_len) {
        for (std::uint32_t i = 0; i < str_len; ++i) {
            if (std::isspace(str[i])) {
                return i;
            }
        }
        return str_len;
    }

private:
    std::unique_ptr<gzFile_s, int(*)(gzFile)> file_;
    std::vector<char> buffer_;
    std::uint32_t buffer_ptr_;
    std::uint32_t buffer_bytes_;
    std::vector<char> storage_;
    std::uint32_t storage_ptr_;
};

#endif //BYTEPAIRRAM_PARSER_H
