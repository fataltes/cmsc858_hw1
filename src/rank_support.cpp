//
// Created by Fatemeh Almodaresi on 2019-11-02.
//
#include "rank_support.h"
#include <vector>
#include <cmath>

/**
 *
 * @param cvec
 */
void Rank_support::construct() {
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
    blocksPerSuperBlock = static_cast<uint32_t >(std::ceil(s/(double)b));
    Rb = new compact::vector<uint64_t>(bWidth, static_cast<uint64_t>(Rs->size()*blocksPerSuperBlock));
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
        (*Rb)[RbIdx] = superBlockOneCnt;
        RbIdx++;
        bitCnt=endOfBlock-i > 0 and endOfBlock-i < b? endOfBlock - i : b;
        auto wrd = cvec.get_int(i, bitCnt);
        blockOneCnt = __builtin_popcount(wrd);
        superBlockOneCnt += blockOneCnt;
        accumOneCnt += blockOneCnt;
        i+=bitCnt;
    }
}

uint64_t Rank_support::rank1(uint64_t i) {
    uint64_t val = 0;
    if (i >= cvec.size()) {
        std::cerr << "ERROR! Index requested is out of range. Index: " << i << " bitvector size: " << cvec.size() << "\n";
        std::exit(5);
    }
    uint64_t RsIdx = i/s;
    uint64_t RbIdx = RsIdx*blocksPerSuperBlock + (i%s)/b;

    auto blockStartOnBV = i-(i%s)%b;
    auto nextBlockStart = std::min(cvec.size(),(RsIdx+1)*s);
    auto bits = nextBlockStart - blockStartOnBV != 0 and nextBlockStart - blockStartOnBV < b ? nextBlockStart - blockStartOnBV : b;
    auto pType = cvec.get_int(blockStartOnBV, bits);
    uint64_t RpIdx = pType*b+(i%s)%b;
    val = (*Rs)[RsIdx] + (*Rb)[RbIdx] + (*Rp)[RpIdx];
//    std::cerr << i << ":rs=" << RsIdx << "," << (*Rs)[RsIdx] << " rb=" << RbIdx << "," << (*Rb)[RbIdx]
//    << " rp=" << RpIdx << "," << (*Rp)[RpIdx] << " blockStart=" << blockStartOnBV << " bits=" << bits << "\n";
    return val;
}

uint64_t Rank_support::rank0(uint64_t i) {
    return i+1 - rank1(i);
}

uint64_t Rank_support::overhead() {
    return Rs->bits() + Rb->bits() + Rp->bits();
}

int benchMarkRank(Opts& opts) {
    compact::vector<uint64_t, 1> cvec(211);
    cvec.clear_mem();
    std::cerr << "Started benchmarking rank .. " << cvec.size() <<"\n";
    for (uint64_t i = 0; i < cvec.size(); i+=10) {
        cvec[i] = 1;
    }
    std::cerr << "Done filling out the bv.\n";
    Rank_support r(cvec);
    std::cerr << cvec.size() << "\n";
    auto prevRank = 0;
    for (uint64_t i = 0; i < cvec.size(); i++) {
        auto rank = r(i);
        if (rank != prevRank) {
            std::cerr << i << ":" << r(i) << "\n";
        }
        prevRank = rank;
    }
    return EXIT_SUCCESS;
}
int benchMarkSelect(Opts& opts) {}
int benchMarkWaveletTrees(Opts& opts) {}
