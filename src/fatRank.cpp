//
// Created by Fatemeh Almodaresi on 2019-11-02.
//
#include "fatRank.h"
#include <vector>
#include <cmath>

void FatRank::doIt(compact::vector<uint64_t, 1> &cvec) {
    uint64_t wrdCnt = 64;
    uint64_t n = cvec.size();
    auto s = static_cast<uint32_t>(std::pow(std::log2(n), 2) / 2);
    auto b = static_cast<uint32_t>(std::log2(n) / 2);
    auto p = static_cast<uint32_t>(std::ceil(std::log2(b)));
//    std::vector<uint64_t> Rs, Rb;
    compact::vector<uint64_t> Rs(s, static_cast<uint64_t>(std::ceil(n/(double)s)));
    Rs.clear_mem();
    compact::vector<uint64_t> Rb(b, static_cast<uint64_t>(std::ceil(n/(double)b)));
    Rb.clear_mem();
    compact::vector<uint64_t> Rp(p, static_cast<uint64_t>(std::ceil(n/(double)p)));
    uint64_t i{0}, j{0}, nextBlock{0}, nextSuperBlock{0}, bitCnt{0}, accumOneCnt{0}, superBlockOneCnt{0}, RsIdx{0}, RbIdx{0}, RpIdx{0};

    // Filling Rs
    while (i < n) {
        bitCnt=std::min(wrdCnt, nextSuperBlock-i);
        auto wrd = cvec.get_int(i, bitCnt);
        auto ones = __builtin_popcount(wrd);
        accumOneCnt += ones;
        i += bitCnt;
        if (i == nextSuperBlock) {
            Rs[RsIdx++] = accumOneCnt;
            j++;
            nextSuperBlock = std::min(n-1, j*s);
            superBlockOneCnt = 0;
        }
    }

    // Filling Rb and Rp
    i = 0;
    while (i < n) {
        bitCnt=std::min(b, static_cast<uint32_t >(nextSuperBlock-i));
        auto wrd = cvec.get_int(i, bitCnt);
        for (auto pCnt = b-1; pCnt >= 0; pCnt--) {
            Rp[RpIdx++] = __builtin_popcount(wrd >> pCnt);
        }
        superBlockOneCnt += __builtin_popcount(wrd);
        i += bitCnt;
        if (i == nextBlock) {
            Rb[RbIdx++] = superBlockOneCnt;
        }
        if (i == nextSuperBlock) {
            superBlockOneCnt = 0;
        }
    }
}

int benchMarkRank(Opts& opts) {}
int benchMarkSelect(Opts& opts) {}
int benchMarkWaveletTrees(Opts& opts) {}
