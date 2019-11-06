//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#include <select_support.h>

#include "select_support.h"

uint64_t Select_support::select1(uint64_t i) {
    if (i == 0) {
        std::cerr << "Warning! select works on values greater than 0.\n";
        std::exit(5);
    }
    auto n = r.getBvSize();
    uint64_t s{0}, e{n};
    return recursiveSelect(s, e, i);
}

uint64_t Select_support::recursiveSelect(uint64_t s, uint64_t e, uint64_t g) {
    auto m = (e+s)/2;
    if (m == r.getBvSize()) {
        std::cerr << "Warning: Select input > total # of 1s in the bv. returning the bv size.\n";
        return r.getBvSize();
    }
    auto rank = r(m);
    if (rank == g) {
        return r.getSetIdxLessEqual(m);
    } else if (rank > g) {
        recursiveSelect(s, m-1, g);
    } else {
        recursiveSelect(m+1, e, g);
    }
}


uint64_t Select_support::select0(uint64_t i) {
    return 0;
}

uint64_t Select_support::overhead() {
    return 0;
}

int benchmarkSelect(Opts &opts) {
    compact::vector<uint64_t, 1> cvec(211);
    cvec.clear_mem();
    std::cerr << "Started benchmarking rank .. " << cvec.size() <<"\n";
    for (uint64_t i = 0; i < cvec.size(); i+=10) {
        cvec[i] = 1;
    }
    std::cerr << "Done filling out the bv.\n";
    std::cerr << cvec.size() << "\n";
    Rank_support r(cvec);
    Select_support s(r);
    for (uint64_t i = 1; i < cvec.size()/10 + 5; i++) {
        auto select = s(i);
        std::cerr << i << ":" << select << "\n";
    }
    return EXIT_SUCCESS;
}
