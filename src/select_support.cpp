//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#include <chrono>

#include "select_support.h"

int64_t Select_support::select1(uint64_t i) {
    if (i == 0) {
        std::cerr << "Error! select works on values greater than 0.\n";
        std::exit(5);
    }
    if (r(r.getBvSize()-1) < i) return -1;
    auto n = r.getBvSize();
    uint64_t s{0}, e{n};
    return recursiveSelect(s, e, i);
}

int64_t Select_support::select0(uint64_t i) {
    if (i == 0) {
        std::cerr << "Error! select works on values greater than 0.\n";
        std::exit(5);
    }
    if (r.rank0(r.getBvSize()-1) < i) return -1;
    auto n = r.getBvSize();
    uint64_t s{0}, e{n};
    return recursiveSelect0(s, e, i);
}

int64_t Select_support::recursiveSelect(uint64_t s, uint64_t e, uint64_t g) {
    auto m = (e+s)/2;
    if (m == r.getBvSize()) {
//        std::cerr << "Warning: Select input > total # of 1s in the bv. returning the bv size.\n";
        return BVOperators::INVALID;
    }
    if (s > e) {
//        std::cerr << g << " not found\n";
        return BVOperators::INVALID;
    }
    auto rank = r(m);
    if (rank == g) {
        return r.getSetIdxLessEqual(m, 1);
    } else if (rank > g) {
        recursiveSelect(s, m-1, g);
    } else {
        recursiveSelect(m+1, e, g);
    }
}

int64_t Select_support::recursiveSelect0(uint64_t s, uint64_t e, uint64_t g) {
    auto m = (e+s)/2;
    if (m == r.getBvSize()) {
//        std::cerr << "Warning: Select input > total # of 1s in the bv. returning the bv size.\n";
        return BVOperators::INVALID;
    }
    if (s > e) {
//        std::cerr << g << " not found\n";
        return BVOperators::INVALID;
    }

    auto rank = r.rank0(m);
    if (rank == g) {
        return r.getSetIdxLessEqual(m, 0);
    } else if (rank > g) {
        recursiveSelect(s, m-1, g);
    } else {
        recursiveSelect(m+1, e, g);
    }
}

uint64_t Select_support::overhead() {
    return 0;
}

int benchmarkSelect(Opts &opts) {
    std::ofstream o;
    if (opts.prefix != "console") {
        std::string filename = opts.prefix + "/" + BVOperators::selectStatFileName;
        o.open(filename, std::ios::out);
        o << "bv_size\tselect_size\tavg_select_time\n";
    } else {
        std::cout << "bv_size\tselect_size\tavg_select_time\n";
    }
    for (auto v = opts.minBVSize; v < opts.maxBVSize; v+=opts.jumpSize) {
        std::cerr << "V" << v << "\n";
        compact::vector<uint64_t, 1> cvec(v);
        cvec.clear_mem();
        for (uint64_t i = 0; i < cvec.size(); i+=10) {
            cvec[i] = 1;
        }
        Rank_support r(cvec);
        Select_support s(r);
        auto start = std::chrono::high_resolution_clock::now();
        for (uint32_t t = 1; t < cvec.size()/10; t++) {
            s(t);
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
    /*compact::vector<uint64_t, 1> cvec(40000000);
    cvec.clear_mem();
    std::cerr << "Started benchmarking rank .. " << cvec.size() <<"\n";
    for (uint64_t i = 0; i < cvec.size(); i+=11) {
        cvec[i] = 1;
    }
    std::cerr << "Done filling out the bv.\n";
    std::cerr << cvec.size() << "\n";
    Rank_support r(cvec);
    Select_support s(r);
    for (uint64_t i = 1; i < cvec.size()/10 + 5; i++) {
        auto select = s(i);
        std::cerr << i << ":" << select << "\n";
    }*/
    return EXIT_SUCCESS;
}
