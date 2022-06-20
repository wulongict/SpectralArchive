//
// Created by wulong on 11/21/17.
//


#include "CSocketServer.h"
#include "CSpectralArchive.h"
#include "../psmvalidator/dependentlibs/msbasic/CDebugMode.h"
#include "../psmvalidator/dependentlibs/msbasic/commonfun.h"
#include "dependentlibs/spectralIndex/CKMeans.h"
#include <fstream>

using namespace std;

boost::program_options::variables_map getParam(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description patternLR("Command line only"), configs_infile("Config file");
    patternLR.add_options()
            ("config", po::value<string>(),
             "the config file")
            ("help,h", po::bool_switch()->default_value(false), "print help information");

    configs_infile.add_options()
            ("mzxmlfiles,m", po::value<string>()->default_value("/data/wulong/data/honeybee/all_mzXML.txt"),
             "the text file containing all the mzXML files")
            ("pepxmls,p", po::value<string>()->default_value("/data/wulong/data/honeybee/all_pepxml.txt"),
             "the text file containing all the pepxml files by X!Tandem or Comet")
            ("indexfile,i", po::value<string>()->default_value("/data/wulong/data/honeybee/faissX.index"),
             "if index does not exist, create index file using [mzxmlfiles], otherwise, load and use it directly ")
            ("rebuild,r", po::bool_switch()->default_value(false),
             "build indexfile again! Ignore and overwrite the existing index file!")
            ("update", po::bool_switch()->default_value(false),
             "update index with new mzXML and new pepxml. default false.")
            ("updaterawdata", po::value<string>()->default_value(""),
             "update index with new mzXML")
            ("updategt", po::value<string>()->default_value(""),
             "update peptide annotation table with input file (pepxml)")
            ("indexstrs", po::value<string>()->default_value(
                    "IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16"),
             "the string for creating indices.")
            ("updategtlist", po::value<string>()->default_value(""),
             "update peptide annotation table with input file list (pepxmllist)")
            ("use_gpu,g", po::bool_switch()->default_value(false),
             "use gpu or not")

            ("removeprecursor", po::bool_switch()->default_value(false),
             "remove precursor from spectrum, clear peaks in range [mz - 17, mz + 3] ")
            ("useflankingbins", po::bool_switch()->default_value(false),
             "add peak to flanking bins with half of the intensity ")
            ("inputsource", po::value<string>()->default_value("cmd"),
             "the source of query: 1. cmd; 2. socket; 3. stdin. The cmd will only search once! "
             "while the socket and stdin will open a while loop and keep accepting input;")
            ("queryindex", po::value<int>()->default_value(1000),
             "input absolute queryindex from mzxmlfiles. This is for check the groundtruth. One should always find itself with any query; ")
            ("topn", po::value<int>()->default_value(10),
             "return top N of the search result. Default N=10")
            ("port", po::value<int>()->default_value(8700),
             "The port which socket to listen to. Default: port=8700")
            ("maxconnection", po::value<int>()->default_value(10),
             "The port which socket to listen to. Default: port=10; "
             "check /proc/sys/net/core/somaxconn for maximal number of connections allowed on local computer")
            ("centroidinit", po::value<int>()->default_value(0),
             "the method to generate initial centroids (used in kmeans)")
            ("indexstring", po::value<string>()->default_value("IVF100,PQ8"),
             "the string for index_factory() to build index");

    patternLR.add(configs_infile);
    boost::program_options::variables_map vm;
    po::positional_options_description p;
    p.add("inputfile", -1);
    po::store(po::command_line_parser(argc, argv).options(patternLR).positional(p).run(), vm);

    if (vm.count("config")) {
        std::ifstream ifs{vm["config"].as<string>().c_str()};
        if (ifs) {
            po::store(parse_config_file(ifs, configs_infile), vm);
        }
    }

    po::notify(vm);

    if (vm.at("help").as<bool>() or argc == 1) {
        cout << "This is a commandline tool for building an searchable index of millions of spectra. "
                "To make this possible, we create a compressed version of spectra vectors, where each "
                "spectra vector comes from binning of all those peaks. " << endl << endl
             << "The tool itself can be used as a website. It uses socket to process requests from browser. "
                "A spectral network is created based on k Nearest Neighbours." << endl << endl;
        //cout << "Build: " << __DATE__ << " " << __TIME__ << endl;
        cout << patternLR << endl;
        cout << "Help info finish" << endl;
        exit(0);
    }
    displayParam(vm);

    return vm;
}


void displayTitle() {
    spdlog::get("A")->info("\n"
                           "-------------------------------------------------\n"
                           "#       Spectral Archive Search Engine          #\n"
                           "#       Build Date: {} {}        #\n"
                           "#     Developed in Henry Lam's Group @ HKUST    #\n"
                           "-------------------------------------------------", __DATE__, __TIME__);
}

int main(int argc, char *argv[]) {
    try {
        initlog("spectral_clustering.log", "A");
        spdlog::set_level(spdlog::level::debug);
        displayTitle();

        auto vm = getParam(argc, argv);
        string indexfilename = vm.at("indexfile").as<string>();
        string mzXMLList = vm.at("mzxmlfiles").as<string>();
        string pepxmls = vm.at("pepxmls").as<string>();
        bool rebuild = vm.at("rebuild").as<bool>();
        bool use_gpu = vm.at("use_gpu").as<bool>();
        string indexFactoryStr = vm.at("indexstring").as<string>();
        string inputsource = vm.at("inputsource").as<string>();

        int topn = vm.at("topn").as<int>();
        bool update_index = vm.at("update").as<bool>();
        string new_experimental_data = vm.at("updaterawdata").as<string>();
        string new_search_result = vm.at("updategt").as<string>();
        string new_search_result_list = vm.at("updategtlist").as<string>();

        int port_listening = vm["port"].as<int>();
        int max_connection_allowed = vm["maxconnection"].as<int>();

        bool removeprecursor = vm["removeprecursor"].as<bool>();
        bool useflankingbins = vm["useflankingbins"].as<bool>();
        int initcenoption = vm.at("centroidinit").as<int>();
        string indexstrs = vm.at("indexstrs").as<string>();

        CPQParam option;
        option._option = static_cast<CPQParam::PQ_OPTIONS>(initcenoption);
        int tol = 15;
        int minpeaknum = 3;
        int topPeakNum = 50, seed = 0;
        CSpectralArchive archive(mzXMLList, pepxmls, indexfilename, removeprecursor, useflankingbins, tol, minpeaknum,
                                 false, option, indexstrs, use_gpu, rebuild, seed, topPeakNum, false, false, false);

        if (update_index) {
            archive.update(new_experimental_data, new_search_result, new_search_result_list, "");
            spdlog::get("A")->info("Index updating finished!");

        } else {

            if (inputsource == "socket") {
                spdlog::get("A")->info("Starting index server via socket.. ");
                archive.setnProbe(256);
                CSocketServer(max_connection_allowed, port_listening, archive, topn);
            }

            cout << "Releasing gpu, spectable and index " << endl;
        }
    }
    catch (const exception &ex) {
        cout << "Program exit due to exception: " << ex.what() << endl;
    }
    return 0;
}

