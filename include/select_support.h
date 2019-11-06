//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#ifndef BVOPERATORS_SELECTSUPPORT_H
#define BVOPERATORS_SELECTSUPPORT_H

#include "rank_support.h"

class Select_support {
public:
    explicit Select_support(Rank_support &rank_supportIn): r(rank_supportIn) {}
    int64_t operator()(uint64_t i) {return select1(i);}
    int64_t select1(uint64_t i);
    int64_t select0(uint64_t i);
    int64_t recursiveSelect(uint64_t s, uint64_t e, uint64_t g);
    int64_t recursiveSelect0(uint64_t s, uint64_t e, uint64_t g);
    uint64_t overhead();
private:
    Rank_support &r;
};


#endif //BVOPERATORS_SELECTSUPPORT_H
