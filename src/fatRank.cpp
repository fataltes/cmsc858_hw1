//
// Created by Fatemeh Almodaresi on 2019-11-02.
//
#include "fatRank.h"
#include <vector>
#include <cmath>

/**
 *
 * @param cvec
 */
void FatRank::construct() {
    uint64_t n = cvec.size();
    std::cerr << "cvec size: " << n << "\n";
    s = static_cast<uint32_t>(std::pow(std::log2(n), 2) / 2);
    std::cerr << "s = " << s << " ";
    /*if (s < 8) {
        s = 8;
    } else if (s < 16) {
        s = 16;
    } else if (s < 32) {
        s = 32;
    } else if (s < 64) {
        s = 64;
    }
    std::cerr << "word-rounded s = " << s << "\n";
     */
    b = static_cast<uint32_t>(std::log2(n)/2);
    std::cerr << "b = " << b << " ";
    p = static_cast<uint32_t>(std::ceil(std::log2(b)));
    std::cerr << "p = " << p << "\n";

//    std::vector<uint64_t> Rs, Rb;
    auto sWidth = static_cast<uint32_t >(std::ceil(std::log2(n)));
    Rs = new compact::vector<uint64_t>(sWidth, static_cast<uint64_t>(std::ceil(n/(double)s)));
    Rs->clear_mem();
    auto bWidth = static_cast<uint32_t >(std::ceil(std::log2(sWidth)));
    Rb = new compact::vector<uint64_t>(bWidth, static_cast<uint64_t>(std::ceil(n/(double)b)));
    Rb->clear_mem();
    auto pWidth = static_cast<uint32_t >(std::ceil(std::log2(bWidth)));
    Rp = new compact::vector<uint64_t>(pWidth, static_cast<uint64_t>(std::pow(2, b) * b));
    Rp->clear_mem();
    std::cerr << "Rs.size(): " << Rs->size() << " Rb.size(): " << Rb->size() << " Rp.size(): " << Rp->size() << "\n";

    // Filling Rp
    for (auto i = 0; i < Rp->size(); i++) {
        uint32_t accumRank = 0;
        for (auto j = 0; j < b; j++) {
            accumRank += (i >> j) & 1;
//            std::cerr << i << " " << j << " " << accumRank << "\n";
            (*Rp)[i*b + j] = accumRank;
        }
    }

    // Filling Rs and Rb
    uint32_t i{0}, bitCnt{0}, accumOneCnt{0}, superBlockOneCnt{0}, blockOneCnt{0}, endOfBlock{s}, RsIdx{0}, RbIdx{0};

    while (i < n) {
        if (i % s == 0) {
            (*Rs)[RsIdx++] = accumOneCnt;
            superBlockOneCnt = 0;
            endOfBlock = RsIdx*s;
        }
        bitCnt=endOfBlock-i > 0 and endOfBlock-i < b? endOfBlock - i : b;
        auto wrd = cvec.get_int(i, bitCnt);
        blockOneCnt = __builtin_popcount(wrd);
        (*Rb)[RbIdx++] = superBlockOneCnt;
        superBlockOneCnt += blockOneCnt;
        accumOneCnt += blockOneCnt;
        i+=bitCnt;
    }
}

uint64_t FatRank::rank1(uint64_t i) {
    uint64_t val = 0;
    if (i >= cvec.size()) {
        std::cerr << "ERROR! Index requested is out of range. Index: " << i << " bitvector size: " << cvec.size() << "\n";
        std::exit(5);
    }
    auto pType = cvec.get_int((i/b)*b, b);
    val = (*Rs)[i/s] + (*Rb)[i/b] + (*Rp)[pType*b+i%b];
//    std::cerr << i << ":rs=" << (*Rs)[i/s] << " rb=" << (*Rb)[i/b] << " rp=" << (*Rp)[pType*b+i%b] << " pType=" << pType << "\n";
    return val;
}

uint64_t FatRank::rank0(uint64_t i) {
    return i+1 - rank1(i);
}

uint64_t FatRank::overhead() {
    return Rs->bits() + Rb->bits() + Rp->bits();
}

int benchMarkRank(Opts& opts) {
    compact::vector<uint64_t, 1> cvec(256);
    cvec.clear_mem();
    std::cerr << "Started benchmarking rank .. " << cvec.size() <<"\n";
    for (uint64_t i = 0; i < cvec.size(); i+=10) {
        cvec[i] = 1;
    }
    std::cerr << "Done filling out the bv.\n";
    FatRank r(cvec);
    std::cerr << cvec.size() << "\n";
    for (uint64_t i = 0; i < cvec.size(); i++) {
        std::cerr << i << ":" << r(i) << "\n";
    }
    return EXIT_SUCCESS;
}
int benchMarkSelect(Opts& opts) {}
int benchMarkWaveletTrees(Opts& opts) {}
