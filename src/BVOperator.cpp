//
// Created by Fatemeh Almodaresi on 2019-11-02.
//


#include "clipp.h"
//#include "spdlog/spdlog.h"
//#include "spdlog/fmt/ostr.h"
//#include "spdlog/fmt/fmt.h"
#include "opts.h"
#include "rank_support.h"

//#include "CLI/Timer.hpp"
using namespace clipp;

int benchmarkRank(Opts &opts);
int benchmarkSelect(Opts &opts);
int constructWaveletTree(Opts &opts);
int operateOnWaveletTree(Opts &opts);

int main(int argc, char* argv[])  {
    (void) argc;
    Opts opts;
    enum class mode {help, rank, select, wv_construct, wv_operations};
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
    auto wvConstructMode = (
            command("wv_build").set(selected, mode::wv_construct),
                    (option("-i", "--input_file") & value("inputFile", opts.inputFile)) % "The file containing the input sequence",
                    (option("-p", "--index_prefix") & value("indexPrefix", opts.prefix)) % "The directory to store the index (default:cout in console)"
    );
    std::string type = "access";
    auto wvOperationMode = (
            command("wv").set(selected, mode::wv_operations),
                    (option("-t", "--type") & value("OptType", type)) % "options:(access, rank, select), default:access",
                    (option("-i", "--input_file") & value("inputFile", opts.inputFile)) % "The file containing the input queries",
                    (option("-p", "--index_prefix") & value("indexPrefix", opts.prefix)) % "The parent directory of the index"
    );


    //Multithreaded console logger(with color support)
//    auto console = spdlog::stderr_color_mt("console");

    bool showHelp = false;
    auto cli = (
            (rankMode | selectMode |
            wvConstructMode | wvOperationMode |
                    command("help").set(selected,mode::help) |
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
            case mode::rank:
                benchmarkRank(opts); std::cerr << "Done\n"; break;
            case mode::select:
                benchmarkSelect(opts);  break;
            case mode::wv_construct:
                constructWaveletTree(opts); break;
            case mode::wv_operations:
                if (type == "access") {
                    opts.operation = Operation::acc;
                } else if (type == "rank") {
                    opts.operation = Operation::rnk;
                } else if (type == "select") {
                    opts.operation = Operation::sel;
                } else {
                    std::cerr << "Undefined Operation: " << type << "\n";
                    std::cerr << "Please choose from these 3 options: (access, rank, select)\n";
                    std::exit(5);
                }
                operateOnWaveletTree(opts); break;
            case mode::help: std::cout << make_man_page(cli, "bvOperate"); break;
        }
    }
    return EXIT_SUCCESS;
}