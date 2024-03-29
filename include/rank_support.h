//
// Created by Fatemeh Almodaresi on 2019-11-02.
//

#ifndef BVOPERATORS_RANKSUPPORT_H
#define BVOPERATORS_RANKSUPPORT_H

#include <ostream>

#include "opts.h"
#include "compact_vector/compact_vector.hpp"


class Rank_support {
public:
    explicit Rank_support(compact::vector<uint64_t, 1> &cvecIn) :
    cvec(cvecIn) , Rs(0) , Rb (0) , Rp(0) {
        construct();
    }

    uint64_t operator()(uint64_t i) {return rank1(i);}
    uint64_t rank1(uint64_t i);
    uint64_t rank0(uint64_t i);
    uint64_t overhead();

    uint64_t getBvSize() const;
    uint64_t getSetIdxLessEqual(uint64_t i, uint64_t v);

private:
    uint64_t bvSize;
    uint32_t s = 0;
    uint32_t b = 0;
    uint32_t p = 0;
    uint32_t blocksPerSuperBlock = 0;
    compact::vector<uint64_t, 1> &cvec;
    compact::vector<uint64_t> Rs;
    compact::vector<uint64_t> Rb;
    compact::vector<uint64_t> Rp;
    void construct();
};

#endif //BVOPERATORS_RANKSUPPORT_H
