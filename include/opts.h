//
// Created by Fatemeh Almodaresi on 2019-11-03.
//

#ifndef BVOPERATORS_OPTS_H
#define BVOPERATORS_OPTS_H

#include <cstdint>

struct Opts {
    uint64_t minBVSize{10000};
    uint64_t maxBVSize{1000000};
    uint64_t jumpSize{100000};
};

#endif //BVOPERATORS_OPTS_H
