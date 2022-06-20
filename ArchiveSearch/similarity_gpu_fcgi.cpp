//
// Created by wulong on 11/14/18.
//

#include "dependentlibs/msbasic/commonfun.h"
#include "dependentlibs/msbasic/CDebugMode.h"
#include "../librarymsms/CMzFileReader.h"
#include "FastCgiInterface.h"
#include "dependentlibs/spectralIndex/CKMeans.h"
#include "CSpectralArchive.h"
#include <thread>

using namespace std;

boost::program_options::variables_map getParam(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description patternLR("Command line only"), configs_infile("Config file");
    patternLR.add_options()
            ("config", po::value<string>(),
             "the config file")
            ("help,h", po::bool_switch()->default_value(false), "print help information");
    configs_infile.add_options()
            ("mzxmlfiles,m", po::value<string>()->default_value(""),
             "the text file containing all the mzXML files")
            ("pepxmls,p", po::value<string>()->default_value(""),
             "the text file containing all the pepxml files by X!Tandem or Comet")
            ("indexfile,i", po::value<string>()->default_value("./"),
             "if index does not exist, create index file using [mzxmlfiles], otherwise, load and use it directly ")
            ("rebuild,r", po::bool_switch()->default_value(false),
             "build indexfile again! Ignore and overwrite the existing index file!")
             ("verbose,v", po::bool_switch()->default_value(false),
             "print information for debugging . ")
            ("update", po::bool_switch()->default_value(false),
             "update index with new mzXML and new pepxml. default false.")
            ("updaterawdata", po::value<string>()->default_value(""),
             "update index with new mzXML")
            ("updaterawdatalist", po::value<string>()->default_value(""),
             "update index with new mzXMLlist, this is the file with many datafiles")
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
            ("recall_true_neighbor", po::bool_switch()->default_value(false),
             "calculate the recall of true neighbor (defined as most similar spectrum by dot product score) ")
            ("recallTrueNeighborTopK", po::value<int>()->default_value(10),
             " Consider top K true neighbors; default 10  ")
            ("recallTrueNeighborMinDP", po::value<double>()->default_value(0.6),
             "minimum dot product allowed for filter true nearest neighbors; default 0.6  ")
            ("minDpOfNeighborRecordedInSqlDB", po::value<double>()->default_value(0.5),
             "minimum dot product allowed for filter true nearest neighbors; default 0.5 ")
            ("newimp", po::bool_switch()->default_value(false),
             "use our newly implemented index ")
            ("debug_cpu_gpu", po::bool_switch()->default_value(false),
             "compare cpu and gpu result")
            ("inputsource", po::value<string>()->default_value("cmd"),
             "the source of query: 1. cmd; 2. socket; 3. stdin. The cmd will only search once! "
             "while the socket and stdin will open a while loop and keep accepting input;")
            ("queryindex", po::value<int>()->default_value(1000),
             "input absolute queryindex from mzxmlfiles. This is for check the groundtruth. One should always find itself with any query; ")
            ("numprobe", po::value<int>()->default_value(256),
             "number of clusters to be inspected. default: 256 ")
            ("numPvalueCalculator", po::value<int>()->default_value(1),
             "the different pvalues to be calculated. Recommended 1, or 3. default: 1 ")
            ("saveBackgroundScore", po::value<bool>()->default_value(false),
             "only evaulated the top hit; default: false ")
             ("skipBackgroundScoreCalc", po::value<bool>()->default_value(false),
                 "skip background score calculation on command line search mode")
            ("minPeakNum", po::value<int>()->default_value(6),
             " The value is used in visualization  ")
            ("tolerance", po::value<int>()->default_value(1),
             "1 for high-res MS2 and 15 for low-res MS2. default: 1 ")
            ("tophit", po::value<bool>()->default_value(false),
             "only evaulated the top hit; default: false ")
             ("search_in_range_of_archive", po::value<bool>()->default_value(false),
             "search for spectrum in archive with range: first to last; default: false ")
            ("createFileNameBlackList", po::value<bool>()->default_value(false),
             "use filename to build black list in search result, avoid using node to annotate itself.  default: false ")
             // the logic is not easy to follow.
            ("first", po::value<long>()->default_value(0),
             "the first query in data file ")
            ("bgspecseed", po::value<int>()->default_value(0),
             "the seed for generate random background spectra list for pvalue calculation ")
            ("last", po::value<long>()->default_value(-1),
             "the last query in data file; -1 means all ")
            ("centroidinit", po::value<int>()->default_value(0),
             "the method to generate initial centroids (used in kmeans)")
            ("searchBatchSize", po::value<int>()->default_value(2000),
             "the size of batch during search. Default 2000")
            ("topn", po::value<int>()->default_value(10),
             "return top N of the search result. Default N=10")
            ("port", po::value<int>()->default_value(8700),
             "The port which socket to listen to. Default: port=8700")
            ("maxconnection", po::value<int>()->default_value(10),
             "The port which socket to listen to. Default: port=10; "
             "check /proc/sys/net/core/somaxconn for maximal number of connections allowed on local computer")

            ("datafile", po::value<string>()->default_value(""),
             "the datafile for search; file formats supported: mzXML, mzML")
            ("validation", po::value<string>()->default_value(""),
             "the search file from another search engine; Supported file: interact-X.ipro.pep.xml (iProphet) or interact-X.pep.xml (PeptideProphet)")
            // todo: the option below should be removed!
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

// NOTE: How to use the tool
// curl -d "Chris WU" localhost:8001/id
// The order matters; Start nginx first, then start fcgi file
// sudo nginx -c ~/bitbucket/codejam/tools/Release/spectral_archive_nginx_fcgi.conf -s reload
// spawn-fcgi -p 8000 -n -- chriswu_fcgi
// cd ~/bitbucket/codejam/tools/Release
// spawn-fcgi -p 8710 -n -- fastcgi_similarity_gpu.fcgi --config ./similarityNetTool.conf
// usage: https://www.systutorials.com/docs/linux/man/1-spawn-fcgi/
// get a mixture spectra:
// http://spec.ust.hk:8709/id/3704218?TOPN=30

// run the web appliaction: ./scripts/webinterface_8710.bash

int main(int argc, char *argv[]) {
    try {
        initlog("spectral_clustering.log", "A");
        spdlog::set_level(spdlog::level::debug);
        displayTitle();
        spdlog::get("A")->info("CMD: {}", argToStr(argc, argv));

        auto vm = getParam(argc, argv);
        string indexfilename = vm.at("indexfile").as<string>();
        string mzXMLList = vm.at("mzxmlfiles").as<string>();
        string pepxmls = vm.at("pepxmls").as<string>();
        bool rebuild = vm.at("rebuild").as<bool>();
        bool verbose = vm.at("verbose").as<bool>();
        bool use_gpu = vm.at("use_gpu").as<bool>();
        string indexFactoryStr = vm.at("indexstring").as<string>();
        string inputsource = vm.at("inputsource").as<string>();
        int topn = vm.at("topn").as<int>();
        int numPvalueCalculator = vm.at("numPvalueCalculator").as<int>();
        bool saveBackgroundScore = vm.at("saveBackgroundScore").as<bool>();
        int bgspecseed = vm.at("bgspecseed").as<int>();
        int numprobe = vm.at("numprobe").as<int>();
        bool update_index = vm.at("update").as<bool>();
        string new_experimental_data = vm.at("updaterawdata").as<string>();
        string new_experimental_datalist = vm.at("updaterawdatalist").as<string>();
        string new_search_result = vm.at("updategt").as<string>();
        string new_search_result_list = vm.at("updategtlist").as<string>();
        string datafile = vm.at("datafile").as<string>();
        string searchfile = vm.at("validation").as<string>();
        bool newImp = vm.at("newimp").as<bool>();
        bool tophit = vm.at("tophit").as<bool>();
        bool search_range = vm.at("search_in_range_of_archive").as<bool>();
        long first = vm.at("first").as<long>();
        long last = vm.at("last").as<long>();

        int port_listening = vm["port"].as<int>();
        int max_connection_allowed = vm["maxconnection"].as<int>();

        bool removeprecursor = vm["removeprecursor"].as<bool>();
        bool useflankingbins = vm["useflankingbins"].as<bool>();

        bool debug_cpu_gpu = vm["debug_cpu_gpu"].as<bool>();
        bool skipBackgroundScoreCalc=vm["skipBackgroundScoreCalc"].as<bool>();
        int initcenoption = vm.at("centroidinit").as<int>();
        string indexstrs = vm.at("indexstrs").as<string>();
        int tolerance = vm["tolerance"].as<int>();
        int minPeakNum = vm["minPeakNum"].as<int>();
        int searchBatchSize = vm["searchBatchSize"].as<int>();
        bool recallTrueNeighbor = vm["recall_true_neighbor"].as<bool>();
        int recallTNNtopK = vm["recallTrueNeighborTopK"].as<int>();
        double recallTNNMinDP = vm["recallTrueNeighborMinDP"].as<double>();
        double minDpOfNeighborRecordedInSqlDB = vm["minDpOfNeighborRecordedInSqlDB"].as<double>();
        bool createFileNameBlackList = vm["createFileNameBlackList"].as<bool>();


        CPQParam option;
        option._option = static_cast<CPQParam::PQ_OPTIONS>(initcenoption);
        const int topPeakNum = 50;

        CSpectralArchive archive(mzXMLList, pepxmls, indexfilename, removeprecursor, useflankingbins, tolerance,
                                 minPeakNum, newImp, option, indexstrs, use_gpu, rebuild, bgspecseed, topPeakNum,
                                 createFileNameBlackList, saveBackgroundScore, verbose);
        // July 2 2019:  fixed a bug for GPU index, the setnprobe and toGPU are not interchangable. First to setnprobe();
        archive.setnProbe(numprobe);

        if (update_index) {
            archive.update(new_experimental_data, new_search_result, new_search_result_list, new_experimental_datalist);
            spdlog::get("A")->info("Index updating finished!");
        } else {
            if (inputsource == "socket") {
                //nginx
                spdlog::get("A")->info("Starting index server via socket.. ");
                shared_ptr<CSocketServerSummary> socketSummary = make_shared<CSocketServerSummary>();
                CFastCGIServer fastcgiserver(max_connection_allowed, port_listening, archive, topn,0, socketSummary);
                fastcgiserver.startFastCGIServer();

            }
            else if (inputsource == "msocket") {
                //nginx
                spdlog::get("A")->info("Starting index server via socket.. ");
                vector<shared_ptr<CFastCGIServer>> multiServer;
                int threadNum = 10;

                shared_ptr<CSocketServerSummary> socketSummary = make_shared<CSocketServerSummary>();
                for(int i = 0; i < threadNum; i ++){
                    multiServer.push_back(make_shared<CFastCGIServer>(max_connection_allowed, port_listening, archive, topn, i,socketSummary));
                }
                // create threads
                vector<thread> threads;

                for(int i = 0; i < threadNum; i ++){
                    threads.push_back(std::thread(std::bind(&CFastCGIServer::startFastCGIServer, multiServer[i])));
                }

                for(int i = 0; i < threadNum; i ++){
                    threads[i].join();
                }

            }
            else {
                if (search_range){
                    archive.searchNeighborsWithin(minDpOfNeighborRecordedInSqlDB, first, last);
                }
                else if (datafile.empty()) {

                    cout << "datafile should not be empty" << endl;
                    throw runtime_error("--datafile should not be empty!");
                }
                else {
                    double mzTol = 2 * tolerance * 2000.0 / 65535;
                    CMzFileReader mzfile(datafile, false, false, true, mzTol, minPeakNum, verbose);

                    archive.searchMzFileInBatch(mzfile, first, last, searchfile, topn, numPvalueCalculator,
                                                recallTrueNeighbor,
                                                searchBatchSize, bgspecseed, recallTNNtopK, recallTNNMinDP,
                                                skipBackgroundScoreCalc, useflankingbins);
                }
            }
        }
        spdlog::get("A")->info("done");
        return 0;
    }
    catch (const exception &ex) {
        cout << "[Error] " << ex.what() << " program will exit!" << endl;
        spdlog::get("A")->info("program exit with error: {}", ex.what());
        return -1;
    }

}
