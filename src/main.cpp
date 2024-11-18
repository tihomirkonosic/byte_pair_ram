#include <iostream>
#include "../include/parser.h"
#include "../include/Kmer.h"
#include "../include/byte_pair_encoder.h"

void Help() {
    std::cout <<
              "usage: BytePairRam <target>\n";
}

int main(int argc, char * argv[]) {
    if (argc == 1) {
        Help();
        return 0;
    }

    auto tparser = Parser::Create(argv[1]);
    if (tparser == nullptr) {
        return 1;
    }

    std::vector<std::unique_ptr<NucleicAcid>> targets;
    try {
        targets = tparser->Parse(1ULL << 32);
    } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    if (targets.empty()) {
        return 1;
    }

    BytePairEncoder encoder;
    std::vector<Kmer> kmers = encoder.Encode(targets.front(), 20);

    std::cout << "Kmers: " << kmers.size() << std::endl;
}
