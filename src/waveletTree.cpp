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
        loadIdx(inputFile);
    } else {
        initializeWVTree(inputFile);
        wv = new compact::vector<uint64_t, 1>(static_cast<uint64_t >(std::ceil(std::log2(chars.size())) * seqLen));
        wv->clear_mem();
        construct(inputFile);
    }
}

bool WaveletTree::loadIdx(std::string &prefix) {
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
    wv = new compact::vector<uint64_t, 1>(0);
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
        srank[i] = r->rank1(spos[i]-1);
    }
    std::cerr << "Index loaded successfully.\n";
    return true;
}

bool WaveletTree::initializeWVTree(std::string &fileName) {
    uint64_t bufferSize = 100000;
    std::vector<char> buffer(bufferSize, 0);
    std::ifstream in( fileName, std::ios::binary | std::ios::ate);
    auto fileLen = in.tellg();
    if (!fileLen)  {
        std::cerr << "ERROR! Sequence file is either empty or corrupted.\n";
        std::exit(3);
    }
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
        v.second = idx++;
    }
    // we're done with reading the file at this point

    for (int64_t i = lowestLevelStart - 1; i >= 0; i--) {
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
            insertIntoWVRecursively(chars[c], 0);
        }
        cursor += bufferSize;
    }
    // last piece of the file
    bufferSize = seqLen - (cursor-bufferSize);
    buffer.resize(bufferSize);
    in.read(buffer.data(), bufferSize);
    for (auto c : buffer) {
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

/**
 *
 * @param idx between 0 and n, length of the original sequence (inclusive)
 * @return the character at index idx
 */
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

/**
 * A Top to Bottom bit-rank call over the levels of the tree
 *
 * @param c a character in the alphabet defined for wavelet tree
 * @param idx between 0 and n, length of the sequence that wavelet tree was constructed over (both inclusive)
 * @return the rank of character c up to and including position idx
 */
int64_t WaveletTree::rank(char c, uint64_t idx) {
    if (chars.find(c) == chars.end()) {
        std::cerr << "Character " << c << " invalid. skipping...\n";
        return BVOperators::INVALID;
    }
    uint64_t charId = chars[c];
    uint64_t blockIdx{0}, blockStart{0}, vIdx{idx}, v{0};
    int64_t offset{0};
    for (uint64_t level = 0; level < charLen; level++) {
        // load value at current level
        v = (charId >> (charLen-level-1)) & 1;
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
    return offset+1; // offset is 0-based but rank should give the count and should be 1-based
}

/**
 * A Bottom to Top bit-select call over the levels of the tree
 *
 * @param c a character in the alphabet defined for wavelet tree
 * @param idx greater than 1 and less than the length of the sequence
 * @return INVALID (-1) if total number of occurrences of the character c is less than idx
 * returns the index i of the sequence that character c has occurred idx times up to and including that
 */
int64_t WaveletTree::select(char c, uint64_t idx) {
    if (chars.find(c) == chars.end()) {
        std::cerr << "Character " << c << " invalid. skipping.\n";
        return BVOperators::INVALID;
    }
    if (idx == 0) {
        std::cerr << "Error! select works on values greater than 0.\n";
        std::exit(5);
    }
    int64_t childPos = 0;
    uint64_t charId = chars[c];
//    std::cerr << "char " << c << " charId=" << charId << "\n";
    for (int64_t level = charLen-1; level >= 0; level--) {
        uint64_t rIdx = charId >> (charLen-level-1);
        uint64_t lastBit = rIdx & 1;
        auto blockIdx = static_cast<uint64_t >(std::pow(2, level)-1 + (rIdx >> 1));
        auto start = spos[blockIdx];
        // if there is no next block, then end is end of the wvt bitvector
        auto end = blockIdx+1 == spos.size()? r->getBvSize() : spos[blockIdx+1];
//        std::cerr << "level=" << level << " idx=" << idx << " rIdx=" << rIdx << " lastBit=" << lastBit
//        << " rIdx >> 1=" << (rIdx >> 1) << " blockIdx=" << blockIdx
//        << " start=" << start << " end=" << end;
        if (lastBit) {
//            std::cerr << " 1=" << idx << " " << srank[blockIdx] << " " << idx + srank[blockIdx];
            childPos = s->recursiveSelect(start, end, idx + srank[blockIdx]);
        } else {
//            std::cerr << " 0=" << idx << " " << spos[blockIdx] << " " << srank[blockIdx] << " " << idx + (spos[blockIdx] - srank[blockIdx]);
            childPos = s->recursiveSelect0(start, end, idx + (spos[blockIdx] - srank[blockIdx]));
        }
//        std::cerr << "\n childPos=" << childPos << "\n";
        if (childPos == BVOperators::INVALID) {
            return BVOperators::INVALID;
        }
        idx = childPos - spos[blockIdx] + 1; // +1 is for converting the index to count ("i"th 1/0 to ("i"+1) 1s/0s)
    }
    if (idx == 0) {
        return BVOperators::INVALID;
    }
   return idx-1;
}



/***
 * Main functions that will be called from the main main!
 */

/**
 * Constructs The wavelet tree
 *
 * @param opts
 * @return
 */
int constructWaveletTree(Opts &opts) {
    WaveletTree wv(opts.inputFile);
    wv.serialize(opts.prefix);
    WaveletTree wv2(opts.prefix, true);
    opts.prefix = "console";
    wv2.serialize(opts.prefix);
}

/**
 * Operates select, rank, or access of the queries in the input file over the wavelet tree index
 *
 * @param opts
 * @return
 */
int operateOnWaveletTree(Opts &opts) {
    WaveletTree wv(opts.prefix, true);
//    std::string console = "console";
//    wv.serialize(console);
    std::ifstream query(opts.inputFile, std::ios::in);
    uint64_t idx;
    char c;
    if (opts.operation == Operation::acc) {
        std::cout << "\n\n\nResults of ACCESS operation:\n\n";
        while (query.good()) {
            query >> idx;
            if (query.good())
                std::cout << idx << ":" << wv.access(idx) << "\n";
        }
    } else if (opts.operation == Operation::rnk) {
        std::cout << "\n\n\nResults of RANK operation:\n\n";
        while (query.good()) {
            query >> c >> idx;
            if (query.good()) {
                auto res = wv.rank(c, idx);
                if (res >= 0) {
                    std::cout << "character " << c << " happens " << res << " times up to index " << idx << "\n";
                }
            }
        }
    } else if (opts.operation == Operation::sel) {
        std::cout << "\n\nResults of SELECT operation:\n";
        while (query.good()) {
            query >> c >> idx;
            if (query.good()) {
                auto res = wv.select(c, idx);
                if (res >= 0) {
                    std::cout << idx << (idx % 10 == 1?"st ":(idx % 10 == 2)?"nd ":"th ") << c << " happens at index " << res << "\n";
                } else {
                    std::cout << "character " << c << " occurs < " << idx << " times\n";
                }
            }
        }
    }
}
