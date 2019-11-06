//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#ifndef BVOPERATORS_WAVELETTREE_H
#define BVOPERATORS_WAVELETTREE_H

#include <string>

#include "select_support.h"

class WaveletTree {

public:
    WaveletTree(Opts &opts);

    bool initializingWVTree(std::string &fileName);
    bool serialize();

private:
    std::string indexPrefix;
    uint64_t seqLen{0};
    uint32_t charLen{0};
    std::map<char, uint64_t> chars;
    compact::vector<uint64_t, 1> *wv;
    std::vector<uint64_t> spos;
    bool construct(std::string &fileName);
    void insertIntoWVRecursively(uint64_t c, uint64_t level);
};


#endif //BVOPERATORS_WAVELETTREE_H
