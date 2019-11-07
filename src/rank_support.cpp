//
// Created by Fatemeh Almodaresi on 2019-11-02.
//
#include "rank_support.h"
#include <vector>
#include <cmath>
#include <chrono>

/**
 *
 * @param cvec
 */
void Rank_support::construct() {
    bvSize = cvec.size();
    std::cerr << "cvec size: " << bvSize << "\n";
    s = static_cast<uint32_t>(std::pow(std::log2(bvSize), 2) / 2);
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
    std::cerr << "word-rounded s = " << s << "\bvSize";
     */
    b = static_cast<uint32_t>(std::log2(bvSize)/2);
    std::cerr << "b = " << b << " ";
    p = static_cast<uint32_t>(std::ceil(std::log2(b)));
    std::cerr << "p = " << p << "\n";

    auto sWidth = static_cast<uint32_t >(std::ceil(std::log2(bvSize)));
    Rs.set_m_bits(sWidth);
    Rs.resize(static_cast<uint64_t>(std::ceil(bvSize/(double)s)));
    Rs.clear_mem();
    auto bWidth = static_cast<uint32_t >(std::ceil(std::log2(sWidth)));
    blocksPerSuperBlock = static_cast<uint32_t >(std::ceil(s/(double)b));
    Rb.set_m_bits(bWidth);
    Rb.resize(static_cast<uint64_t>(Rs.size()*blocksPerSuperBlock));
    Rb.clear_mem();
    // this b+1 was fucking tricky!!!!
    // In this special case the value you put in each cell could actually be EXACTLY EQUAL to b
    // So here we don't want to represent b numbers but present numbers from 0 to fucking b inclusive
    // In cases b is a power of 2 for that one last value b we need log2(b+1) bits!!! rather than log2(b) bits
    auto pWidth = static_cast<uint32_t >(std::ceil(std::log2(b+1)));
    Rp.set_m_bits(pWidth);
    Rp.resize(static_cast<uint64_t>(std::pow(2, b) * b));
    Rp.clear_mem();

    // Filling Rp
    for (auto i = 0; i < Rp.size()/b; i++) {
        uint32_t accumRank = 0;
        for (auto j = 0; j < b; j++) {
            accumRank += (i >> j) & 1;
            Rp[i*b + j] = accumRank;
        }
    }
    // Filling Rs and Rb
    uint32_t i{0}, bitCnt{0}, accumOneCnt{0}, superBlockOneCnt{0}, blockOneCnt{0}, endOfBlock{s}, RsIdx{0}, RbIdx{0};

    while (i < bvSize) {
        if (i % s == 0) {
            Rs[RsIdx++] = accumOneCnt;
            superBlockOneCnt = 0;
            endOfBlock = RsIdx*s;
        }
        Rb[RbIdx] = superBlockOneCnt;
        RbIdx++;
        bitCnt=endOfBlock-i > 0 and endOfBlock-i < b? endOfBlock - i : b;
        bitCnt = std::min(static_cast<uint32_t >(cvec.size()-i), bitCnt);
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
        return BVOperators::INVALID;
//        std::exit(5);
    }
    uint64_t RsIdx = i/s;
    uint64_t RbIdx = RsIdx*blocksPerSuperBlock + (i%s)/b;

    auto blockStartOnBV = i-(i%s)%b;
    auto nextBlockStart = std::min(cvec.size(),(RsIdx+1)*s);
    auto bits = nextBlockStart - blockStartOnBV != 0 and nextBlockStart - blockStartOnBV < b ? nextBlockStart - blockStartOnBV : b;
    auto pType = cvec.get_int(blockStartOnBV, bits);
    uint64_t RpIdx = pType*b+(i%s)%b;
    if (RsIdx >= Rs.size() or RbIdx >= Rb.size() or RpIdx >= Rp.size()) {
        std::cerr << "Error! rank=" << i << " RsIdx=" << RsIdx << " RbIdx=" << RbIdx << " RpIdx" << RpIdx << "\n";
        std::exit(3);
    }
    val = Rs[RsIdx] + Rb[RbIdx] + Rp[RpIdx];
//    std::cerr << "\nRs=" << RsIdx << " " << Rs[RsIdx] << " Rb=" << RbIdx << " " << Rb[RbIdx]
//    << " Rp=" << RpIdx << " " << Rp[RpIdx] << " (i%s)/b=" << (i%s)/b << " (i%s)%b=" << (i%s)%b <<
//    " ptype=" << pType << "\n";
    return val;
}

uint64_t Rank_support::rank0(uint64_t i) {
    return i+1 - rank1(i);
}

uint64_t Rank_support::overhead() {
    return (Rs.bytes() + Rb.bytes() + Rp.bytes())*8;
}

uint64_t Rank_support::getBvSize() const {
    return bvSize;
}

uint64_t Rank_support::getSetIdxLessEqual(uint64_t i, uint64_t v) {
    auto idx = i;
    while (true) {
        if (!idx) return 0;
        auto start = static_cast<uint64_t >(std::max(0, static_cast<int>(idx)-63));
        auto bits = idx-start+1;
        auto wrd = cvec.get_int(start, bits);
        for (auto shiftVal = bits-1; shiftVal >= 0; shiftVal--) {
            if ( (wrd >> shiftVal & 1) == v) return idx;
            idx--;
        }
    }
}

int benchmarkRank(Opts &opts) {
    std::ofstream o;
    if (opts.prefix != "console") {
        std::string filename = opts.prefix + "/" + BVOperators::rankStatFileName;
        o.open(filename, std::ios::out);
        o << "bv_size\trank_size\tavg_rank_time\n";
    } else {
        std::cout << "bv_size\trank_size\tavg_rank_time\n";
    }
    for (auto v = opts.minBVSize; v < opts.maxBVSize; v+=opts.jumpSize) {
//        std::cerr << "V" << v << "\n";
        compact::vector<uint64_t, 1> cvec(v);
        cvec.clear_mem();
        for (uint64_t i = 0; i < cvec.size(); i+=10) {
            cvec[i] = 1;
        }
        Rank_support r(cvec);
        auto start = std::chrono::high_resolution_clock::now();
        for (uint32_t t = 0; t < cvec.size(); t++) {
            r(t);
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        if (opts.prefix == "console") {
            std::cout << v << "\t" << r.overhead() << "\t" << elapsed.count() / cvec.size() << "\n";
        } else {
            o << v << "\t" << r.overhead() << "\t" << elapsed.count() / cvec.size() << "\n";
        }
    }
    if (opts.prefix != "console") {
        o.close();
    }

 /*   compact::vector<uint64_t, 1> cvec(1100);
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
            std::cerr << i << ":" << rank << "\n";
        }
        prevRank = rank;
    }*/
    return EXIT_SUCCESS;
}
