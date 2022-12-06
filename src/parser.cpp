//
// Created by tiho on 12/6/22.
//

#include "../include/parser.h"

std::atomic<std::uint32_t> NucleicAcid::num_objects{0};

std::unique_ptr<Parser> Parser::Create(const std::string& path) {
    auto file = gzopen(path.c_str(), "r");
    if (file == nullptr) {
        throw std::invalid_argument(
                "[bioparser::Parser::Create] error: unable to open file " + path);
    }
    return std::unique_ptr<Parser>(new Parser(file, 4194304));
}