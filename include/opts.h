//
// Created by Fatemeh Almodaresi on 2019-11-03.
//

#ifndef BVOPERATORS_OPTS_H
#define BVOPERATORS_OPTS_H

#include <cstdint>

namespace BVOperators {
    constexpr char wvIdxFileName[] = "wv.bin";
    constexpr char idxInfoFileName[] = "idxInfo.bin";
    constexpr char rankStatFileName[] = "rank.stat";
    constexpr char selectStatFileName[] = "select.stat";
}

enum Operation {
    acc, rnk, sel
};
struct Opts {
    uint64_t minBVSize{10000000};
    uint64_t maxBVSize{1000000000};
    uint64_t jumpSize{10000000};
    std::string prefix = "console";
    std::string inputFile;
    Operation operation = Operation::acc;
};

#endif //BVOPERATORS_OPTS_H
