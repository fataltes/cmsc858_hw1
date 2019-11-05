//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#ifndef BVOPERATORS_SELECTSUPPORT_H
#define BVOPERATORS_SELECTSUPPORT_H

#include "rank_support.h"

class Select_support {
public:
    Select_support(Rank_support &rank_supportIn): r(rank_supportIn) {}
    uint64_t operator()(uint64_t i) {return select1(i);}
    uint64_t select1(uint64_t i);
    uint64_t select0(uint64_t i);
    uint64_t recursiveSelect(uint64_t s, uint64_t e, uint64_t g);
    uint64_t overhead();
private:
    Rank_support &r;
};


#endif //BVOPERATORS_SELECTSUPPORT_H
