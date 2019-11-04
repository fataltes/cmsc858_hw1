//
// Created by Fatemeh Almodaresi on 2019-11-02.
//


#include "clipp.h"
//#include "spdlog/spdlog.h"
//#include "spdlog/fmt/ostr.h"
//#include "spdlog/fmt/fmt.h"
#include "opts.h"
#include "fatRank.h"

//#include "CLI/Timer.hpp"
using namespace clipp;

int benchMarkRank(Opts& opts);
int benchMarkSelect(Opts& opts);
int benchMarkWaveletTrees(Opts& opts);

int main(int argc, char* argv[])  {
    (void) argc;
    Opts opts;
    enum class mode {help, rank, select, wv};
    mode selected = mode::help;

    // Rank Mode, prepare the rank performance distribution over bv size
    auto rankMode = (
            command("rank").set(selected, mode::rank),
                    (option("-s", "--startSize") & value("minBVSize", opts.minBVSize)) % "The bv size to start benchmarking with (default:10,000)",
                    (option("-e", "--endSize") & value("maxBVSize", opts.maxBVSize)) % "The max bv size to benchmark (default:1,000,000)",
                    (option("-j", "--jumpSize") & value("deltaSize", opts.jumpSize)) % "The bv size to start with (default: 100,000)"
    );
    auto selectMode = (
            command("select").set(selected, mode::select),
                    (option("-s", "--startSize") & value("minBVSize", opts.minBVSize)) % "The bv size to start benchmarking with (default:10,000)",
                    (option("-e", "--endSize") & value("maxBVSize", opts.maxBVSize)) % "The max bv size to benchmark (default:1,000,000)",
                    (option("-j", "--jumpSize") & value("deltaSize", opts.jumpSize)) % "The bv size to start with (default: 100,000)"
    );
    auto wvMode = (
            command("wavelet").set(selected, mode::wv),
                    (option("-s", "--startSize") & value("minBVSize", opts.minBVSize)) % "The bv size to start benchmarking with (default:10,000)",
                    (option("-e", "--endSize") & value("maxBVSize", opts.maxBVSize)) % "The max bv size to benchmark (default:1,000,000)",
                    (option("-j", "--jumpSize") & value("deltaSize", opts.jumpSize)) % "The bv size to start with (default: 100,000)"
    );

    //Multithreaded console logger(with color support)
//    auto console = spdlog::stderr_color_mt("console");

    bool showHelp = false;
    auto cli = (
            (rankMode | selectMode | wvMode | command("help").set(selected,mode::help) |
            option("--help", "-h").set(showHelp, true) % "show help"
            ));

    decltype(parse(argc, argv, cli)) res;
    try {
        res = parse(argc, argv, cli);
        if (showHelp) {
            std::cout << make_man_page(cli, "bvOperate");
            return 0;
        }
    } catch (std::exception &e) {
        std::cout << "\n\nparsing command line failed with exception: " << e.what() << "\n";
        std::cout << "\n\n";
        std::cout << make_man_page(cli, "bvOperate");
        return 1;
    }

    if(res) {
        switch(selected) {
            case mode::rank: benchMarkRank(opts);  break;
            case mode::select: benchMarkSelect(opts);  break;
            case mode::wv: benchMarkWaveletTrees(opts); break;
            case mode::help: std::cout << make_man_page(cli, "bvOperate"); break;
        }
    }

}