//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#ifndef BVOPERATORS_WAVELETTREE_H
#define BVOPERATORS_WAVELETTREE_H

#include <string>

#include "select_support.h"

class WaveletTree {

public:
    std::map<char, uint64_t> chars;

    explicit WaveletTree(std::string &inputFile, bool loadIndex=false);

    virtual ~WaveletTree();

    bool initializeWVTree(std::string &fileName);
    bool serialize(std::string &prefix);

    char access(uint64_t idx);
    int64_t rank(char c, uint64_t idx);
    int64_t select(char c, uint64_t idx);

    uint64_t size() {
        return wv.bytes();
    }

    bool constructRankSupport();

private:
    std::string indexPrefix;
    uint64_t seqLen{0};
    uint32_t charLen{0};
    std::vector<char> inverseChars;
    compact::vector<uint64_t, 1> wv;
    std::vector<uint64_t> spos;
    std::vector<uint64_t> srank;
    bool construct(std::string &fileName);
    void insertIntoWVRecursively(uint64_t c, uint64_t level);
    bool loadIdx(std::string &prefix);
    Rank_support *r;
    Select_support *s;
};


#endif //BVOPERATORS_WAVELETTREE_H
