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
#include "CTimerSummary.h"
#include "Util.h"

using namespace std;

boost::program_options::variables_map getParam(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description genericOptions("Command line only"), basicOptions("Config file"), hiddenOptions(
            "advancedOptions");
    genericOptions.add_options()
            ("config", po::value<string>(),
             "the config file")
            ("help,h", po::bool_switch()->default_value(false), "print help information");
    basicOptions.add_options()
            ("verbose,v", po::bool_switch()->default_value(false),
             "print information for debugging . ")
            ("archivename,n", po::value<string>()->default_value(""),
             "the name of the spectral archive")
            ("mzxmlfiles,m", po::value<string>()->default_value(""),
             "the text file containing all the mzXML files")

             // -----------options for build archive-------------
            ("indexstrs", po::value<string>()->default_value(
                     "IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16;IVF256,PQ16"),
             "the string for creating indices.")
            ("indexshuffleseeds", po::value<string>()->default_value("default"),
             "the string for creating indices with shuffle seeds. if the value is not default,"
             " a customized seeds list can be: customized:a;b;c;d;e;f when there are six indexes. "
             "The characters a-f should be integers. default equivalent to customized:0;1;2;3;4;5.")

             // option to launch the archive.
            ("inputsource", po::value<string>()->default_value("cmd"),
             "the source of query: 1. socket; 2. msocket; 3. cmd. The cmd will only search the given datafile, "
             "while the socket and msocket will keep accepting queries from client side")
            ("numprobe", po::value<int>()->default_value(256),
             "number of buckets to look into, the web API allows user to change this value in each query. default: 256 ")
            ("use_gpu,g", po::bool_switch()->default_value(false),
             "use gpu or not")

             // options for update existing archive
            ("update", po::bool_switch()->default_value(false),
             "update index with new mzXML and new pepxml. default false.")
            ("updaterawdata", po::value<string>()->default_value(""),
             "update index with new mzXML")
            ("updaterawdatalist", po::value<string>()->default_value(""),
             "update index with new mzXMLlist, this is the file with many datafiles")
            ("updategt", po::value<string>()->default_value(""),
             "update peptide annotation table with input file (pepxml)")
             ("updategtlist", po::value<string>()->default_value(""),
                     "update peptide annotation table with input file list (pepxmllist)")

             // options to shrink existing archive.
            ("popDataFileNum", po::value<int>()->default_value(0),
             "number of data files to be removed from tail of the archive. For nonzero value n, the last n "
             "data files will be removed. default: 0 ")
            ("shrinkto", po::value<int>()->default_value(-1),
             "shrink the archive to a smaller size. specify the number of spectra. Data files will be removed "
             "until the archive is smaller than this size. -1 means no shrink value. default: -1 ")


             // options for searching against existing archive
            ("datafile", po::value<string>()->default_value(""),
             "the datafile for search; file formats supported: mzXML, mzML")
            ("searchBatchSize", po::value<int>()->default_value(2000),
             "the size of batch during search. Default 2000")
            ("topn", po::value<int>()->default_value(10),
             "return top N of the search result. Default N=10")
             ("first", po::value<long>()->default_value(0),
                                 "the first query in data file ")
            ("last", po::value<long>()->default_value(-1),
             "the last query in data file; -1 means all ");
    hiddenOptions.add_options()
            ("pepxmls,p", po::value<string>()->default_value(""),
             "the text file containing all the pepxml files by X!Tandem or Comet")
            ("indexfile,i", po::value<string>()->default_value(""),
             "[DEPRECATED] if index does not exist, create index file using [mzxmlfiles], otherwise, load and use it directly ")
            ("rebuild,r", po::bool_switch()->default_value(false),
             "build indexfile again! Ignore and overwrite the existing index file!")

            // default options use for building archive. should always be true.
            ("removeprecursor", po::bool_switch()->default_value(false),
             "remove precursor from spectrum, clear peaks in range [mz - a, mz + b] ")
            ("useflankingbins", po::bool_switch()->default_value(false),
             "add peak to flanking bins with half of the intensity ")

             // --------------- TNN recall related options --------------------
            ("recall_true_neighbor", po::bool_switch()->default_value(false),
             "calculate the recall of true neighbor (defined as most similar spectrum by dot product score) ")
            ("recallTrueNeighborTopK", po::value<int>()->default_value(10),
             " Consider top K true neighbors; default 10  ")
             ("recallTrueNeighborMinPeakNumInQuery", po::value<int>()->default_value(10),
             " Skip query spectra with very few peaks.  default 10  ")
            ("recallTrueNeighborMinDP", po::value<double>()->default_value(0.6),
             "minimum dot product allowed for filter true nearest neighbors; default 0.6  ")

             // search function implemented inside archive.
            ("search_in_range_of_archive", po::value<bool>()->default_value(false),
             "search for spectrum in archive with range: first to last; default: false ")
            ("minDpOfNeighborRecordedInSqlDB", po::value<double>()->default_value(0.5),
             "minimum dot product allowed for filter true nearest neighbors; default 0.5 ")

            ("debug_cpu_gpu", po::bool_switch()->default_value(false),
             "compare cpu and gpu result")

            ("queryindex", po::value<int>()->default_value(1000),
             "input absolute queryindex from mzxmlfiles. This is for check the groundtruth. One should always find itself with any query; ")


            // 10k background score related options
            ("numPvalueCalculator", po::value<int>()->default_value(1),
             "the different pvalues to be calculated. Recommended 1, or 3. default: 1 ")
            ("saveBackgroundScore", po::value<bool>()->default_value(false),
             "only evaulated the top hit; default: false ")
             ("skipBackgroundScoreCalc", po::value<bool>()->default_value(false),
                 "skip background score calculation on command line search mode")
            ("bgspecseed", po::value<int>()->default_value(0),
             "the seed for generate random background spectra list for pvalue calculation ")

             // search results processing
             // 1. results visualization options
             ("plotHistogram", po::value<bool>()->default_value(false),
             "plot histogram of calibrated ssm score")
             // 2. keep only top hit for filter step analysis
            ("tophit", po::value<bool>()->default_value(true),
             "only evaulated the top hit; default: false ")

             // ---------------------- the parameters should not be changed.
            ("minPeakNum", po::value<int>()->default_value(6),
             " The value is used in visualization  ")
            ("tolerance", po::value<int>()->default_value(1),
             "1 for high-res MS2 and 15 for low-res MS2. default: 1 ")




            ("createFileNameBlackList", po::value<bool>()->default_value(false),
             "use filename to build black list in search result, avoid using node to annotate itself.  default: false ")

            // options related to my implementation of ivfpq; not used
            ("newimp", po::bool_switch()->default_value(false),
             "use our newly implemented index ")
            ("centroidinit", po::value<int>()->default_value(0),
             "the method to generate initial centroids (used in kmeans)")


             // the following options should not be used any more. to be removed.
            ("port", po::value<int>()->default_value(8700),
             "The port which socket to listen to. Default: port=8700")
            ("maxconnection", po::value<int>()->default_value(10),
             "The port which socket to listen to. Default: port=10; "
             "check /proc/sys/net/core/somaxconn for maximal number of connections allowed on local computer")
            ("indexstring", po::value<string>()->default_value("IVF100,PQ8"),
             "the string for index_factory() to build index")


            ("validation", po::value<string>()->default_value(""),
             "the search file from another search engine; Supported file: interact-X.ipro.pep.xml "
             "(iProphet) or interact-X.pep.xml (PeptideProphet)");

    po::options_description cmdline_options;
    cmdline_options.add(genericOptions).add(basicOptions).add(hiddenOptions);
//
    po::options_description config_file_options;
    config_file_options.add(basicOptions).add(hiddenOptions);
//
    po::options_description visible("Allowed options");
    visible.add(genericOptions).add(basicOptions);

//    genericOptions.add(basicOptions);
    boost::program_options::variables_map vm;
    boost::program_options::variables_map visiable_vm;
    po::positional_options_description p;
    p.add("inputfile", -1);
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

    po::store(po::command_line_parser(argc, argv).options(visible).run(), visiable_vm);

    if (vm.count("config")) {
        std::ifstream ifs{vm["config"].as<string>().c_str()};
        if (ifs) {
            // po::store will not overwrite existing options.
            po::store(parse_config_file(ifs, config_file_options), vm);
            po::store(parse_config_file(ifs, config_file_options), visiable_vm);
        }
    }

    po::notify(vm);
    po::notify(visiable_vm);
//    cout  <<  vm.at("inputsource").as<string>() << endl;
    if (vm.at("help").as<bool>() or argc == 1) {
        cout << "Spectroscape is a open source tool for building an searchable spectral archive. It can be used in both command line and also as a fastcgi service together with a web interface."
                << endl << endl;

        cout << visible << endl;
        cout << "Help info finish" << endl;
        exit(0);
    }
    displayParam(visiable_vm);
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



// run the web appliaction: ./scripts/webinterface_8710.bash


int main(int argc, char *argv[]) {

    try {
        SimpleTimer st(false);
        initlog("spectral_clustering.log", "A");
        spdlog::set_level(spdlog::level::debug);
        displayTitle();
        spdlog::get("A")->info("CMD: {}", argToStr(argc, argv));

        auto vm = getParam(argc, argv);
        string indexfilename = vm.at("indexfile").as<string>();
        string mzXMLList = vm.at("mzxmlfiles").as<string>();
        string archivename = vm.at("archivename").as<string>();
        string pepxmls = vm.at("pepxmls").as<string>();
        bool rebuild = vm.at("rebuild").as<bool>();
        bool verbose = vm.at("verbose").as<bool>();
        bool use_gpu = vm.at("use_gpu").as<bool>();
        string indexFactoryStr = vm.at("indexstring").as<string>();
        string inputsource = vm.at("inputsource").as<string>();
        int topn = vm.at("topn").as<int>();
        int popDataFileNum = vm.at("popDataFileNum").as<int>();
        int shrinkto = vm.at("shrinkto").as<int>();

        int numPvalueCalculator = vm.at("numPvalueCalculator").as<int>();
        bool saveBackgroundScore = vm.at("saveBackgroundScore").as<bool>();
        bool plotHistogram = vm.at("plotHistogram").as<bool>();
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
        string indexshuffleseeds = vm.at("indexshuffleseeds").as<string>();
        int tolerance = vm["tolerance"].as<int>();
        int minPeakNum = vm["minPeakNum"].as<int>();
        int searchBatchSize = vm["searchBatchSize"].as<int>();
        bool recallTrueNeighbor = vm["recall_true_neighbor"].as<bool>();
        int recallTNNtopK = vm["recallTrueNeighborTopK"].as<int>();
        int minPeakNumInSpec = vm["recallTrueNeighborMinPeakNumInQuery"].as<int>();
        double recallTNNMinDP = vm["recallTrueNeighborMinDP"].as<double>();
        double minDpOfNeighborRecordedInSqlDB = vm["minDpOfNeighborRecordedInSqlDB"].as<double>();
        bool createFileNameBlackList = vm["createFileNameBlackList"].as<bool>();


        CPQParam option;
        option._option = static_cast<CPQParam::PQ_OPTIONS>(initcenoption);
        const int topPeakNum = 50;

        // pass seeds in.
        CSpectralArchive archive(mzXMLList, pepxmls, indexfilename, removeprecursor, useflankingbins, tolerance,
                                 minPeakNum, newImp, option, indexstrs, use_gpu, rebuild, bgspecseed, topPeakNum,
                                 createFileNameBlackList, saveBackgroundScore, verbose, archivename, indexshuffleseeds);
        archive.setRecallTNN(recallTrueNeighbor, recallTNNtopK, recallTNNMinDP, minPeakNumInSpec); // the parameter is set here. 
        // July 2 2019:  fixed a bug for GPU index, the setnprobe and toGPU are not interchangable. First to setnprobe();
        archive.setnProbe(numprobe);
        int removeFileNumToShrink = archive.getNumOfFilesToRemoveForShrinkingArchiveTo(shrinkto);

        popDataFileNum = popDataFileNum > removeFileNumToShrink ? popDataFileNum: removeFileNumToShrink;
        if (popDataFileNum > 0 ){
            // remove data get top priority. after remove, program will exit.
            // remove information of the last file
            archive.remove(popDataFileNum);
        }
        else if (update_index) {
            archive.update(new_experimental_data, new_search_result, new_search_result_list, new_experimental_datalist);
            spdlog::get("A")->info("Index updated!");
        }
        else  {
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
                // all the other cases goes here,
                if (search_range){
                    archive.searchNeighborsWithin(minDpOfNeighborRecordedInSqlDB, first, last);
                }
                else if (datafile.empty()) {

                    cout << "datafile should not be empty" << endl;
                    throw runtime_error("--datafile should not be empty!");
                }
                else {

                    double mzTol = 2 * tolerance * 2000.0 / 65535;
                    CMzFileReader mzfile(datafile, false, false, removeprecursor, mzTol, minPeakNum, verbose);
                    mzfile.init(false, false, removeprecursor); // added back.

                    archive.searchMzFileInBatch(mzfile, first, last, searchfile, topn, numPvalueCalculator,
                                                recallTrueNeighbor,
                                                searchBatchSize, bgspecseed, recallTNNtopK, recallTNNMinDP,
                                                skipBackgroundScoreCalc, useflankingbins, tophit, plotHistogram);
                }
            }
        }
        ostringstream oss;
        CTimeSummary::getInstance()->print(oss, archive.getNumOfQueriesSearched());

        spdlog::get("A")->info("Time report:\n{}\n", oss.str());
        double timeused = st.secondsElapsed();
        spdlog::get("A")->info("Spectroscape task is finished. \nThe spectra number in archive is {}.\nTotal time elapsed: {:.4f} seconds", archive.size(), timeused);

        return 0;
    }
    catch (const exception &ex) {

        spdlog::get("A")->info("program exit with error: {}", ex.what());
        return -1;
    }

}
