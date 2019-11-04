//
// Created by Fatemeh Almodaresi on 2019-11-02.
//

#ifndef CMSC858_HW1_FATRANK_H
#define CMSC858_HW1_FATRANK_H

#include "opts.h"
#include "compact_vector/compact_vector.hpp"


class FatRank {
public:
    FatRank(uint64_t size) {
        cvec.resize(size);
    }

    void doIt(compact::vector<uint64_t, 1> &cvec);
private:
    compact::vector<uint64_t, 1> cvec;
};

#endif //CMSC858_HW1_FATRANK_H
