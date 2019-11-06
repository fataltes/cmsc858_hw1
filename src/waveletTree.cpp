//
// Created by Fatemeh Almodaresi on 2019-11-05.
//

#include <vector>
#include <map>
#include <waveletTree.h>
#include <cmath>
#include <compact_vector/compact_vector.hpp>

#include "waveletTree.h"

WaveletTree::WaveletTree(std::string &inputFile, bool loadIndex) {
    if (loadIndex) {
        loadWVTree(inputFile);
    } else {
        initializingWVTree(inputFile);
        wv = new compact::vector<uint64_t, 1>(static_cast<uint64_t >(std::ceil(std::log2(chars.size())) * seqLen));
        wv->clear_mem();
        construct(inputFile);
    }
}

bool WaveletTree::loadWVTree(std::string &prefix) {
    std::cout << "loading the index from " + prefix + "..\n";
    std::ifstream idxInfo(prefix + "/" + BVOperators::idxInfoFileName,
                          std::ios::binary | std::ios::in);
    idxInfo.read(reinterpret_cast<char*>(&seqLen), sizeof(seqLen));
    uint32_t charCnt{0};
    idxInfo.read(reinterpret_cast<char*>(&charCnt), sizeof(charCnt));
    inverseChars.resize(charCnt);
    for (uint32_t i = 0; i < charCnt; i++) {
        char c;
        idxInfo.read(&c, sizeof(c));
        chars[c] = i;
        inverseChars[i] = c;
    }
    charLen = static_cast<uint32_t >(std::ceil(std::log2(charCnt)));
    idxInfo.close();
    std::cerr << "seqLen=" << seqLen << " , charCnt=" << charCnt << " , charBits=" << charLen << "\n";
    std::string fileName = prefix + "/" + BVOperators::wvIdxFileName;
    wv = new compact::vector<uint64_t, 1>(10);
    wv->deserialize(fileName, false);
    r = new Rank_support(*wv);
    s = new Select_support(*r);


    // construct the sPos and sRank vectors
    srank.resize(static_cast<uint64_t >(pow(2, charLen) - 1));
    spos.resize(srank.size());
    std::cerr << "seqLen=" << seqLen << " , srank.size()=" << srank.size() << "\n";
    spos[1] = seqLen;
    srank[1] = r->rank1(spos[1]-1);
    for (uint64_t i = 2; i < srank.size(); i++) {
        uint64_t par = (i-1)/2;
        spos[i] = spos[par] + seqLen;
        if (i % 2  == 0) {
            spos[i] += (spos[par+1]-srank[par+1])-(spos[par]-srank[par]);
        }
//        std::cerr << i << " " << spos[i] << "\n";
        srank[i] = r->rank1(spos[i]-1);
    }
    std::cerr << "Index loaded successfully.\n";
    return true;
}

bool WaveletTree::initializingWVTree(std::string &fileName) {
    uint64_t bufferSize = 100000;
    std::vector<char> buffer(bufferSize, 0);
    std::ifstream in( fileName, std::ios::binary | std::ios::ate);
    auto fileLen = in.tellg();
    if (!fileLen) return false;
    seqLen = static_cast<uint64_t >(fileLen)-1;
    std::cerr << "seqLen = " << seqLen << "\n";
    in.seekg(0);
    uint64_t cursor = bufferSize;
    while (cursor < seqLen) {
        in.read(buffer.data(), bufferSize);
        for (auto c : buffer) {
            if (chars.find(c) == chars.end()) {
                chars[c] = 0;
            }
            chars[c]++;
        }
        cursor += bufferSize;
    }
    // last piece of the file
    bufferSize = seqLen - (cursor-bufferSize);
    buffer.resize(bufferSize);
    in.read(buffer.data(), bufferSize);
    for (auto c : buffer) {
        if (chars.find(c) == chars.end()) {
            chars[c] = 0;
        }
        chars[c]++;
    }
    in.close();

    charLen = static_cast<uint32_t >(std::ceil(std::log2(chars.size())));
    std::cerr << "charCnt=" << chars.size() << " , charBits=" << charLen << "\n";
    spos.resize(static_cast<uint64_t >(pow(2, charLen+1) - 1));
    std::cerr << "spos.size()=" << spos.size() << "\n";
    auto lowestLevelStart = static_cast<uint64_t >(spos.size()-pow(2, charLen));
    uint64_t idx{0};
    for (auto &v : chars) {
        spos[lowestLevelStart+idx] = v.second;
//        std::cerr << "idx=" << idx << " " << lowestLevelStart+idx << " " << spos[lowestLevelStart+idx] << "\n";
        v.second = idx++;
    }
    // we're done with reading the file at this point

    for (int i = lowestLevelStart - 1; i >= 0; i--) {
//        std::cerr << i << "->" << 2*i+1 << "," << 2*i+2 << "\n";
        spos[i] = spos[2*i+1] + spos[2*i+2];
    }
    uint64_t s{0}, len{0};
    for (auto l = 0; l < log2(spos.size()); l++) {
        s = s + len;
        uint64_t prevCnt=l*seqLen, currCnt = spos[s];
        len = (uint64_t)std::pow(2, l);
        for (auto i = s; i < s + len; i++) {
            spos[i] = prevCnt;
            prevCnt += currCnt;
            currCnt = spos[i + 1];
        }
    }
    std::cerr << "Wavelet tree initialized.\n";

    /*for (auto i = 0; i < spos.size(); i++) {
        std::cerr << i << ":" << spos[i] << "\n";
    }*/
    return true;
}

bool WaveletTree::construct(std::string &fileName) {
    std::cerr << "Filling the wavelet tree..\n";
    uint64_t bufferSize = 100000;
    std::vector<char> buffer(bufferSize, 0);
    std::ifstream in(fileName, std::ios::binary | std::ios::in);
    uint64_t cursor = bufferSize;
    while (cursor < seqLen) {
        in.read(buffer.data(), bufferSize);
        for (auto c : buffer) {
//            std::cerr << chars[c] << "\n";
            insertIntoWVRecursively(chars[c], 0);
        }
        cursor += bufferSize;
    }
    // last piece of the file
    bufferSize = seqLen - (cursor-bufferSize);
    buffer.resize(bufferSize);
    in.read(buffer.data(), bufferSize);
    for (auto c : buffer) {
//        std::cerr << chars[c] << "\n";
        insertIntoWVRecursively(chars[c], 0);
    }
    in.close();
    std::cerr << "Wavelet tree fully constructed.\n";
}

void WaveletTree::insertIntoWVRecursively(uint64_t c, uint64_t level) {
    auto rowStartIdx = static_cast<uint64_t >(std::pow(2, level))-1;
    auto bucket = c >> (charLen - level);
    auto idx = rowStartIdx + bucket;
    (*wv)[spos[idx]] = c >> (charLen - level - 1) & 1;
    spos[idx]++;
    if (level+1 == charLen) return;
    insertIntoWVRecursively(c, level+1);
}

bool WaveletTree::serialize(std::string &prefix) {
    if (prefix == "console") {
        std::cout << "writing the wv to console\n";
        for (auto c = 0; c < charLen; c++) {
            for (auto i = 0; i < seqLen; i++) {
                std::cout << ((*wv)[c*seqLen+i]? "1 ":"0 ");
            }
            std::cout << "\n";
        }
    } else {
        std::cout << "serializing the wavelet to " + prefix + "\n";
        std::ofstream idxInfo(prefix + "/" + BVOperators::idxInfoFileName,
                std::ios::binary | std::ios::out);
        idxInfo.write(reinterpret_cast<char*>(&seqLen), sizeof(seqLen));
        auto charCnt = static_cast<uint32_t >(chars.size());
        idxInfo.write(reinterpret_cast<char*>(&charCnt), sizeof(charCnt));
        for (auto &kv : chars) {
            idxInfo.write(&kv.first, sizeof(kv.first));
        }
        idxInfo.close();
        std::ofstream wvIdx(prefix + "/" + BVOperators::wvIdxFileName,
                std::ios::binary | std::ios::out);
        wv->serialize(wvIdx);
        wvIdx.close();
    }
    return true;
}

char WaveletTree::access(uint64_t idx) {
    uint64_t charId = 0;
    uint64_t blockIdx{0}, blockStart{0};
    uint64_t offset = idx, vIdx{blockStart+offset}, v{};
    for (uint64_t level = 0; level < charLen; level++) {
        // load value at current level
        v = (*wv)[vIdx];
        charId = (charId << 1) | v;
        //find the next level index
        if (v) {
            offset = (r->rank1(vIdx) - 1) - srank[blockIdx];
        } else {
            offset = (r->rank0(vIdx) - 1) - (spos[blockIdx]-srank[blockIdx]);
        }
        blockIdx = blockIdx*2+v+1;
        blockStart = spos[blockIdx];
        vIdx = blockStart+offset;

    }
    return inverseChars[charId];
}

uint64_t WaveletTree::rank(char c, uint64_t idx) {
    uint64_t charId = 0;
    uint64_t blockStart = 0;
    uint64_t rnk = idx;
    for (uint64_t level = 0; level < charLen; level++) {
        auto blockIdx = static_cast<uint64_t >(std::pow(2, level)-1+charId);
        blockStart = spos[blockIdx];
        auto vIdx = level*seqLen+blockStart+rnk;
        uint16_t v = (*wv)[vIdx];
        charId = (charId << 1) | v;
        if (v) {
            rnk = r->rank1(vIdx) - srank[blockIdx];
        } else {
            rnk = vIdx + 1 - r->rank1(vIdx) - (blockIdx + 1 - srank[blockIdx]);
        }
    }
    return rnk;
}

uint64_t WaveletTree::select(char c, uint64_t idx) {
    return wv->size();
}






int constructWaveletTree(Opts &opts) {
    WaveletTree wv(opts.inputFile);
    wv.serialize(opts.prefix);
    WaveletTree wv2(opts.prefix, true);
    opts.prefix = "console";
    wv2.serialize(opts.prefix);
}

int operateOnWaveletTree(Opts &opts) {
    WaveletTree wv(opts.prefix, true);
    std::string console = "console";
    wv.serialize(console);
    std::ifstream query(opts.inputFile, std::ios::in);
    uint64_t idx;
    char c;
    if (opts.operation == Operation::acc) {
        std::cout << "Results of ACCESS operation:\n";
        while (query.good()) {
            query >> idx;
            if (query.good())
                std::cout << idx << ":" << wv.access(idx) << "\n";
        }
    } else if (opts.operation == Operation::rnk) {
        std::cout << "Results of RANK operation:\n";
        while (query.good()) {
            query >> c >> idx;
            if (query.good())
                std::cout << idx << ":" << wv.rank(c, idx) << "\n";
        }
    } else if (opts.operation == Operation::sel) {
        std::cout << "Results of SELECT operation:\n";
        while (query.good()) {
            query >> c >> idx;
            if (query.good())
                std::cout << idx << ":" << wv.select(c, idx) << "\n";
        }
    }
}
