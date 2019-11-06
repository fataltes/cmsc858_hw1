//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#ifndef BVOPERATORS_WAVELETTREE_H
#define BVOPERATORS_WAVELETTREE_H

#include <string>

#include "select_support.h"

class WaveletTree {

public:
    explicit WaveletTree(std::string &inputFile, bool loadIndex=false);

    bool initializeWVTree(std::string &fileName);
    bool serialize(std::string &prefix);

    char access(uint64_t idx);
    int rank(char c, uint64_t idx);
    int select(char c, uint64_t idx);

private:
    std::string indexPrefix;
    uint64_t seqLen{0};
    uint32_t charLen{0};
    std::map<char, uint64_t> chars;
    std::vector<char> inverseChars;
    compact::vector<uint64_t, 1> *wv;
    std::vector<uint64_t> spos;
    std::vector<uint64_t> srank;
    bool construct(std::string &fileName);
    void insertIntoWVRecursively(uint64_t c, uint64_t level);
    bool loadIdx(std::string &prefix);
    Rank_support *r;
    Select_support *s;
};


#endif //BVOPERATORS_WAVELETTREE_H
