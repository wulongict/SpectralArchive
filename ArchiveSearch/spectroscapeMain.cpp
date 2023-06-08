//
// Created by wulong on 11/14/18.
//

#include "dependentlibs/msbasic/commonfun.h"
#include "dependentlibs/msbasic/CDebugMode.h"
#include "../librarymsms/CMzFileReader.h"
#include "FastCgiInterface.h"
#include "dependentlibs/spectralIndex/CKMeans.h"
#include "FileUtility.h"
#include "CSpectralArchive.h"
#include <thread>
#include "CTimerSummary.h"
#include "Util.h"
// #include <filesystem>
// namespace fs = std::filesystem;
#include <csignal>
#include <atomic>


#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


using namespace std;



void sigint_handler(int signal) {
    
    g_quit_flag.store(true);
    g_quit_key_stroke += 1;
    std::cout << "\nReceived SIGINT signal: " << g_quit_key_stroke << std::endl;
    if(g_quit_key_stroke>1){
        std::cout << "Force to exit Spectroscape." << std::endl;
        exit(0);
    }
    std::cout << "Please wait for a few seconds for Spectroscape save data and exit safely.\n" << std::endl;
    std::cout << "To stop Spectroscape immediately, press ctrl-c again .\n" << std::endl;
}

boost::program_options::variables_map getParam(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description genericOptions("Command line only"), basicOptions("Config file"), hiddenOptions(
            "advancedOptions");
    genericOptions.add_options()
            ("init", po::bool_switch()->default_value(false),
             "check the folder, generate config file, if not exist")
             ("run", po::bool_switch()->default_value(false),
             "run the archive specified by `archivename`")
             ("datasearchpath", po::value<string>()->default_value("./"),
             "search for data files in the given path")
             ("add", po::bool_switch()->default_value(false),
             "add data file to archive by searching for data files specified by option `--datasearchpath`")
             ("yes", po::bool_switch()->default_value(false),
             "ask user to confirm by default. If --yes is specified, directly overwrite existing archive.")

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
             ("gpuidx", po::value<string>()->default_value(
                     "0"),
             "the list of gpu ids for use. by default use first GPU. To use multiple GPUs, e.g. the first three gpus, use \"0;1;2\" as the value.")

             // option to launch the archive.
            ("inputsource", po::value<string>()->default_value("create"),
             "the source of query: 1. socket; 2. msocket; 3. cmd. 4. create.\n"
             "The option `create` will only build a archive and exit; the option`cmd` will search the given datafile, "
             "the option `socket` and `msocket` will keep accepting queries from client side")
            ("numprobe", po::value<int>()->default_value(8),
             "number of buckets to look into, the web API allows user to change this value in each query. default: 8 ")
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
            ("saveAnns", po::value<bool>()->default_value(false),
             "save ANNs to SQL table. default: false ")

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
             ("wwwroot", po::value<string>()->default_value("./"),
             "path of root folder of spectroscape web UI, config this to the SpectralArchiveWeb/arxiv folder. ")
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

    // po::store(po::command_line_parser(argc, argv).options(visible).run(), visiable_vm);
    bool run = vm["run"].as<bool>();
    bool add = vm["add"].as<bool>();

    if (vm.count("config")) {
        string config_file = vm["config"].as<string>();
        if(not File::isExist(config_file)){
            cout << "config file not exist: " << config_file << endl;
            exit(1);
        }
        fs::path config_file_path  = fs::absolute(config_file);
        cout << "The config file is " << config_file_path << endl;
        std::ifstream ifs{config_file.c_str()};
        if (ifs) {
            // po::store will not overwrite existing options.
            po::store(parse_config_file(ifs, config_file_options), vm);
            // po::store(parse_config_file(ifs, config_file_options), visiable_vm);
        }
    } else if (run or add){
        // if run specified, we will try default conf file. 
        string archivie_name = vm["archivename"].as<string>();
        if(archivie_name == ""){
            // vm["archivename"] = boost::lexical_cast<std::string>(fs::current_path().string());
            archivie_name =fs::current_path().string();
            
        }
        cout << "The archive name is " << archivie_name << endl;
        fs::path archivie_path = fs::absolute(archivie_name);
        fs::path config_file = archivie_path / fs::path("conf/spectroscape_auto.conf");
        
        // vm["config"] = boost::lexical_cast<std::string>(config_file);
        fs::path config_file_path  = fs::absolute(config_file);
        cout << "The config file is " << config_file_path << endl;
        if(fs::exists(config_file_path)){
            cout << "The config file is " << config_file_path << endl;
        }else{
            cout << "The config file is not exist" << config_file_path << endl;
            cout << "Please run `<program-name> --init` first to create default config file" << endl;
            exit(1);
        }
        std::ifstream ifs{config_file.c_str()};
        if (ifs) {
            // po::store will not overwrite existing options.
            po::store(parse_config_file(ifs, config_file_options), vm);
            // po::store(parse_config_file(ifs, config_file_options), visiable_vm);
        }
        cout << "vm refreshed" << vm["mzxmlfiles"].as<string>() << endl;

   
    }

    po::notify(vm);
    // po::notify(visiable_vm);
//    cout  <<  vm.at("inputsource").as<string>() << endl;
    if (vm.at("help").as<bool>() or argc == 1) {
        cout << "Spectroscape is a open source tool for building an searchable spectral archive. It can be used in both command line and also as a fastcgi service together with a web interface."
                << endl << endl;

        cout << visible << endl;
        cout << "Help info finish" << endl;
        exit(0);
    }
    displayParam(vm);
    return vm;
}
// SPECTROSCAPE_VERSION should be defined in CMakeLists.txt
#pragma message("The version of Spectroscape is " SPECTROSCAPE_VERSION ".")
#ifndef SPECTROSCAPE_VERSION
#define SPECTROSCAPE_VERSION "0.0.0"
#endif

// new banner for Spectroscape
void displayTitle() {
    spdlog::get("A")->info("\n\n"
                            ".............................................................................\n"
                           R"(.  _____                     _                                              )"".\n"
                           R"(. /  ___|                   | |                                             )"".\n"
                           R"(. \ `--.  _ __    ___   ___ | |_  _ __  ___   ___   ___  __ _  _ __    ___  )"".\n"
                           R"(.  `--. \| '_ \  / _ \ / __|| __|| '__|/ _ \ / __| / __|/ _` || '_ \  / _ \ )"".\n"
                           R"(. /\__/ /| |_) ||  __/| (__ | |_ | |  | (_) |\__ \| (__| (_| || |_) ||  __/ )"".\n"
                           R"(. \____/ | .__/  \___| \___| \__||_|   \___/ |___/ \___|\__,_|| .__/  \___| )"".\n"
                           R"(.        | |                                                  | |           )"".\n"
                           R"(.        |_|                                                  |_|           )"".\n"
                           ".............................................................................\n"
                           ".                                                                           .\n"
                           ".                          Spectroscape v{}                              .\n"
                           ".                                                                           .\n"
                           ".                    Build Date: {} {}                       .\n"
                           ".              Developed by L. WU  in Henry Lam Group @ HKUST               .\n"
                           ".............................................................................\n\n",SPECTROSCAPE_VERSION, __DATE__, __TIME__);
}



void archive_initialization(bool yes_overwrite, string &archivename, const string &datasearchpath);

string &fix_archive_name(string &archivename);

void
collect_files_before_add(const string &datasearchpath, const string &mzXMLList, string &archivename, bool &update_index,
                         string &new_experimental_data, string &new_experimental_datalist, string &new_search_result,
                         string &new_search_result_list);



void writeConfigFile(string filename, string mzxmlfile){
    File::CFile fileobj(filename);
    const fs::path path(fileobj.path);
    if(fs::exists(path)){
        cout << "Path ready"<< endl;
    }else{
        try{
            fs::create_directories(path);
            cout << "Path not exist, created successfully " << endl;
        }catch (const std::exception & ex) {
            std::cerr << "Fail to create path " << path << " - " << ex.what() << endl;
        }

    }

    ofstream  fout(filename, ios::out);
    fout << "mzxmlfiles= "<< mzxmlfile <<"\n"
 "use_gpu = false\n"
 "removeprecursor = true\n"
 "useflankingbins = true\n"
 "\n"
 "inputsource = socket\n"
 "topn = 10\n"
 "port = 8710\n"
 "maxconnection = 10\n"
 "\n"
 "indexstrs = IVF256,PQ16;IVF256,PQ16\n"
 "indexshuffleseeds = customized:1;2\n"
 "\n"
 "#tolerance low mass accuracy: 1; high mass accuracy: 15\n"
 "# 1 ~ 0.03 Th  15 ~ 0.45 Th\n"
 "tolerance = 1\n"
 "\n"
 "# the background spectra list seed\n"
 "bgspecseed = 42\n"
 "saveBackgroundScore = false" << endl;
    fout.close();
    cout << "config file generate " << filename << endl;

}


// run the web appliaction: ./scripts/webinterface_8710.bash


int main(int argc, char *argv[]) {
    // Set up a signal handler for SIGINT
    std::signal(SIGINT, sigint_handler);
    try
    {
        SimpleTimer st(false);
        spdlog::details::os::create_dir("log");
        initlog("log/spectroscape.log", "A");
        spdlog::set_level(spdlog::level::debug);
        displayTitle();
        try
        {

            spdlog::get("A")->info("maximum number of threads on this computer {}", std::thread::hardware_concurrency());
            spdlog::get("A")->info("CMD: {}", argToStr(argc, argv));

            auto vm = getParam(argc, argv);
            bool init = vm.at("init").as<bool>();
            bool yes_overwrite = vm.at("yes").as<bool>();
            string archivename = vm.at("archivename").as<string>();
            
            string datasearchpath = vm.at("datasearchpath").as<string>();
            if (init){
                archivename = fix_archive_name(archivename);
                archive_initialization(yes_overwrite, archivename, datasearchpath);
                return 0;
            }
            bool run = vm.at("run").as<bool>();
            if(run){
                archivename = fix_archive_name(archivename);
            }
            


            string indexfilename = vm.at("indexfile").as<string>();
            string mzXMLList = vm.at("mzxmlfiles").as<string>();
            cout << "The mzXML file is " << mzXMLList << endl;
            
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
            bool add = vm.at("add").as<bool>();
            if(add){
                collect_files_before_add(datasearchpath, mzXMLList, archivename, update_index, new_experimental_data,
                                         new_experimental_datalist,
                                         new_search_result, new_search_result_list);

            }


            string datafile = vm.at("datafile").as<string>();
            string searchfile = vm.at("validation").as<string>();
            bool newImp = vm.at("newimp").as<bool>();
            bool tophit = vm.at("tophit").as<bool>();
            bool search_range = vm.at("search_in_range_of_archive").as<bool>();
            bool saveAnns = vm.at("saveAnns").as<bool>();
            long first = vm.at("first").as<long>();
            long last = vm.at("last").as<long>();

            int port_listening = vm["port"].as<int>();
            int max_connection_allowed = vm["maxconnection"].as<int>();

            bool removeprecursor = vm["removeprecursor"].as<bool>();
            bool useflankingbins = vm["useflankingbins"].as<bool>();

            bool debug_cpu_gpu = vm["debug_cpu_gpu"].as<bool>();
            bool skipBackgroundScoreCalc = vm["skipBackgroundScoreCalc"].as<bool>();
            int initcenoption = vm.at("centroidinit").as<int>();
            string indexstrs = vm.at("indexstrs").as<string>();
            string wwwroot = vm.at("wwwroot").as<string>();
            string indexshuffleseeds = vm.at("indexshuffleseeds").as<string>();
            string gpuidx = vm.at("gpuidx").as<string>();
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
                                     createFileNameBlackList, saveBackgroundScore, verbose, archivename, indexshuffleseeds, gpuidx);
            archive.setRecallTNN(recallTrueNeighbor, recallTNNtopK, recallTNNMinDP, minPeakNumInSpec); // the parameter is set here.
            // July 2 2019:  fixed a bug for GPU index, the setnprobe and toGPU are not interchangable. First to setnprobe();
            archive.setnProbe(numprobe);
            int removeFileNumToShrink = archive.getNumOfFilesToRemoveForShrinkingArchiveTo(shrinkto);

            popDataFileNum = popDataFileNum > removeFileNumToShrink ? popDataFileNum : removeFileNumToShrink;
            if (popDataFileNum > 0)
            {
                archive.remove(popDataFileNum);
            }
            else if (update_index)
            {
                archive.update(new_experimental_data, new_search_result, new_search_result_list, new_experimental_datalist);
                spdlog::get("A")->info("Index updated!");
            }
            else
            {
                if (inputsource == "socket")
                {
                    // nginx
                    spdlog::get("A")->info("Starting index server via socket.. ");
                    shared_ptr<CSocketServerSummary> socketSummary = make_shared<CSocketServerSummary>();
                    CFastCGIServer fastcgiserver(max_connection_allowed, port_listening, archive, topn, 0, socketSummary, wwwroot);
                    fastcgiserver.startFastCGIServer();
                }
                else if (inputsource == "msocket")
                {
                    // nginx
                    spdlog::get("A")->info("Starting index server via socket.. ");
                    vector<shared_ptr<CFastCGIServer>> multiServer;
                    int threadNum = 10;

                    shared_ptr<CSocketServerSummary> socketSummary = make_shared<CSocketServerSummary>();
                    for (int i = 0; i < threadNum; i++)
                    {
                        multiServer.push_back(make_shared<CFastCGIServer>(max_connection_allowed, port_listening, archive, topn, i, socketSummary, wwwroot));
                    }
                    // create threads
                    vector<thread> threads;

                    for (int i = 0; i < threadNum; i++)
                    {
                        threads.push_back(std::thread(std::bind(&CFastCGIServer::startFastCGIServer, multiServer[i])));
                    }

                    for (int i = 0; i < threadNum; i++)
                    {
                        threads[i].join();
                    }
                }
                else if (inputsource == "cmd")
                {
                    // all the other cases goes here,
                    if (search_range)
                    {
                        archive.searchNeighborsWithin(minDpOfNeighborRecordedInSqlDB, first, last, searchBatchSize);
                    }
                    else if (datafile.empty())
                    {

                        cout << "datafile should not be empty" << endl;
                        throw runtime_error("--datafile should not be empty!");
                    }
                    else
                    {

                        double mzTol = 2 * tolerance * 2000.0 / 65535;
                        CMzFileReader mzfile(datafile, false, false, removeprecursor, mzTol, minPeakNum, verbose);
                        mzfile.init(false, false, removeprecursor); // added back.

                        archive.searchMzFileInBatch(mzfile, first, last, searchfile, topn, numPvalueCalculator,
                                                    recallTrueNeighbor,
                                                    searchBatchSize, bgspecseed, recallTNNtopK, recallTNNMinDP,
                                                    skipBackgroundScoreCalc, useflankingbins, tophit, plotHistogram);
                    }
                }else if (inputsource == "create"){
                    // do nothing, just create archive and exit. 
                    spdlog::get("A")->info("spectral archive created.");
                }
            }
            ostringstream oss;
            CTimeSummary::getInstance()->print(oss, archive.getNumOfQueriesSearched());

            spdlog::get("A")->info("Time report:\n{}\n", oss.str());
            double timeused = st.secondsElapsed();
            spdlog::get("A")->info("Spectroscape archive size: {}. Time elapsed: {:.1f} s.", archive.size(), timeused);

            return 0;
        }
        catch (const exception &ex)
        {

            spdlog::get("A")->info("program exit with error: {}", ex.what());

            return -1;
        }
    }
    catch (const exception &ex)
    {
        std::cout << "Error: " << ex.what() << endl;
        return -1;
    }
    catch (...)
    {
        std::cout << "Error: program exit with unknown error " << endl;
    }
}

void
collect_files_before_add(const string &datasearchpath, const string &mzXMLList, string &archivename, bool &update_index,
                         string &new_experimental_data, string &new_experimental_datalist, string &new_search_result,
                         string &new_search_result_list) {
    // checking parameters, if not in current path, change as current path.
    // archivename = fix_archive_name(archivename);

    cout << "adding data " << endl;
    fs::path archive_path(archivename);
    // in this task, we add new data file to the archive.
    if(update_index){
//                    cout << "Error:\t--add and --update cannot be used together. " << endl;
        throw runtime_error("options --add and --update can not be used together. ");
    }
    // collect data from data path.
    vector<string> data_extensions = {".mzml", ".mzxml", ".mgf", ".sptxt", ".scanlist"};
    fs::path fs_datasearchpath = fs::absolute(datasearchpath);
    vector<fs::path> new_experimental_datafiles = get_file_list(fs_datasearchpath, data_extensions);
    vector<string> searchfile_extensions = {".ipro.pep.xml", ".sptxt", ".spectroscape.tsv", "pep.xml", "pepxml"};
    vector<fs::path> new_search_result_files = get_file_list(fs_datasearchpath, searchfile_extensions);
    // write the data file list to a temp file.
    new_experimental_data="";
    new_search_result="";
    new_experimental_datalist = (fs::absolute(archive_path) / fs::path("new_experimental_datafile_list.txt")).string();
    new_search_result_list = (fs::absolute(archive_path) / fs::path("new_search_result_file_list.txt")).string();

    sort(new_experimental_datafiles.begin(), new_experimental_datafiles.end());
    sort(new_search_result_files.begin(), new_search_result_files.end());
    ofstream fout(new_experimental_datalist);
    for(auto & file: new_experimental_datafiles){
        fout << file.string() << endl;
    }
    fout.close();
    fout.open(new_search_result_list);
    for(auto & file: new_search_result_files){
        fout << file.string() << endl;
    }
    fout.close();
    // update the index.
    update_index = true;
    cout << "start to update " << mzXMLList << " archive is " << archivename << endl;
    cout << "new experimental data file list is " << new_experimental_datalist << endl;
    cout << "new search result file list is " << new_search_result_list << endl;
    // output file number
    cout << "new experimental data file number is " << new_experimental_datafiles.size() << endl;
    cout << "new search result file number is " << new_search_result_files.size() << endl;
}

// if archive_name is empty, change it to current folder.
string &fix_archive_name(string &archivename) {// checking parameters, if not in current path, change as current path.
    if(archivename==""){
        archivename = fs::current_path().string();
    }
    return archivename;
}

// three steps to initialization.
// 1. create folder with archive name.
// 2. collecting data files into archive file.
// 3. creating configuration files.
// throw runtime_error when there are error inside.
void archive_initialization(bool yes_overwrite, string &archivename, const string &datasearchpath) {
    // archivename = fix_archive_name(archivename);

    // check if the archive name is already exist.
    if(fs::exists(archivename)){
        cout << "archive path '"<< archivename<<"' already exist" << endl;

    }else{
        bool archive_folder_created = fs::create_directory(archivename);
        if(not archive_folder_created){
            throw runtime_error("No such file or directory: " + archivename);
        }
    }


    archivename = fs::absolute(archivename).string();

    // put raw files into a text file called archive under the archive folder.
    fs::path archive_path(archivename);
    fs::path archive_file = archive_path / "archive";

    if(fs::exists(archive_file) and not fs::is_regular_file(archive_file) ){
        throw runtime_error("can not access file " + archive_file.string() + ", please make sure it is a regular file.");
    }

    if(fs::exists(archive_file) and fs::is_regular_file(archive_file) and not yes_overwrite){
        cout << "archive file already exist, type yes to continue and overwrite the existing file. type no to quit"
        << endl << "(yes/no):" << flush;
        string input;
        cin >> input;
        if (input == "yes"){
            cout << "overwrite the existing archive file" << endl;
        }
        else{
//                        cout << "quit" << endl;
            throw runtime_error("initialization already done. ");
        }

    }

    // look for mzXML and mzML files in the datasearch path.
    fs::path fs_datasearchpath = fs::absolute(datasearchpath);
    auto datafiles = get_file_list(fs_datasearchpath, {".mzxml", ".mzml"});
    // auto datafiles = get_mz_file_list(datasearchpath);
    // sort the datafiles
    sort(datafiles.begin(), datafiles.end());

    // confirmed with yes or yes-overwrite is used.
    int number_of_files = datafiles.size();
    if(number_of_files == 0){
        cout << "Error:\tno mzXML or mzML files found in the search path. " << endl;
        cout <<"\tPlease make sure `--datasearchpath` is set to a folder with mzXML/mzML data files. " << endl;
        cout << "\tCurrent `--datasearchpath` is set to '" << datasearchpath <<"'"<< endl;
        cout << "\tSpectroscape will search in sub-folders. " << endl;
        throw runtime_error("No raw files specified  in path " + datasearchpath);

    }
    ofstream fout(archive_file.string());
    for(auto & file: datafiles){
        fout << file.string() << endl;
    }


    // generate config file if not exist.
    fs::path config_file = archive_path / "conf/spectroscape_auto.conf";
    if(fs::exists(config_file)){
        cout << "config file " << config_file << " already exist"<< endl;
    }else{
        writeConfigFile(config_file.string(), archive_file.string());
    }
    cout << "the archive is initialized successfully. " << endl;


}
