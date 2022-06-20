//
// Created by wulong on 2/26/18.
//


#include <string>
#include "spdlog/spdlog.h"
//#include <spdlog/spdlog.h>

#include "../../../librarymsms/ProteomicsDataTypes.h"
#include "../../../librarymsms/XMLFileParser.h"
#include "../../../librarymsms/Util.h"

#include "ConcretePSMFeatures.h"
#include "CDebugMode.h"
#include "CFlow.h"
#include "FragmentationScore.h"
#include "MGFReader.h"
#include "CPeakPairsImporter.h"
#include "classifier.h"
#include "CFragScore.h"
#include <regex>
#include <mutex>
#include <memory>
#include "BasePSMFeature.h"

#include "../../../External/SpectraST/SpectraSTPeakList.hpp"
//#include "../../../gnuplot-iostream/gnuplot-iostream.h"
#include "../../../External/gnuplot-iostream/gnuplot-iostream.h"
#include "../../../librarymsms/CThreadsPool.h"
#include <random>
#include <stdexcept>

using namespace std;


CFlow::CFlow() { m_name = "Base flow"; }

CFlow::~CFlow() {}

void CFlow::run() {}

void CFlow::attach(CFlow *viewer) {
    following_steps.push_back(viewer);
}

void CFlow::subscribe(CFlow *previous_step) {
    previous_step->attach(this);
}

void CFlow::notify(const string &msg) {
    spdlog::get("A")->info("Send message: {}", msg);
    for (auto x: following_steps) {
        x->update(msg);
    }
}

void CFlow::update(string msg) {
    spdlog::get("A")->info("Get message: {}", msg);
    vector<string> tokens;
    split_string(msg, tokens);
    for (int i = 0; i < tokens.size() - 1 and tokens.size() >= 2; i += 2) {
        update_with_key_value_pair(tokens[i], tokens[i + 1]);
    }
}

string CFlow::getname() const {
    return m_name;
}

FeatureWorkFlow::FeatureWorkFlow(string inputfile, bool ghost, double minInt, string &fragmodelfile,
                                 string fragscoretype, string annotationMethod, bool debug, const string &basename,
                                 string featurelistfile, bool isLowMassAcc, string binaryPath) : CFlow() {
    m_binaryPath = binaryPath;
    m_isLowMassAcc = isLowMassAcc;
    m_name = "PSM feature extraction";
    m_featurelistfile = featurelistfile;
    string featurenameBase = File::CFile(m_featurelistfile).basename;
    m_inputfile = inputfile;
    m_useGhostPeaks = ghost;
    m_minInt = minInt;
    m_fragmodelfile = fragmodelfile;
    m_fragscoretype = fragscoretype;
    m_annotationMethod = annotationMethod;
    File::CFile filename(m_inputfile);
    string tprefix = filename.path + "/target_" + filename.basename;
    string dprefix = filename.path + "/decoy_" + filename.basename;

    m_searchresult = dprefix + ".pep.xml";

    m_debug = debug;
    if (m_annotationMethod == "title") {
        m_feature_outputname = to_string(m_inputfile, "_", featurenameBase, ".feat");
        m_name += " target";
    } else {
        m_feature_outputname = to_string(m_searchresult, "_", featurenameBase, ".feat");
        m_name += " decoy";
    }
}

void FeatureWorkFlow::run() {
    spdlog::get("A")->info("PSM Feature Extraction workflow started!");
    if (not File::isExist(m_fragmodelfile)) {
        cout << "[Info] fragmentation model could not found: \"" << m_fragmodelfile << "\"" << endl;
        exit(-1);
    }
    if (m_annotationMethod == "searchresult" and not File::isExist(m_searchresult)) {
        cout << "[Info] search result file could not be found: \"" << m_searchresult << "\"" << endl;
        exit(-1);
    }

    bool isTruthKnown = true;
    // Step 1. Read MGF file with title annotation...
    MGFReader mgfReader(m_inputfile, m_isLowMassAcc, isTruthKnown);

    bool verbosity = m_debug ? true : false;
    vector<Feature *> features = createFeatureVector(m_fragmodelfile, m_useGhostPeaks, m_minInt, m_fragscoretype,
                                                     m_featurelistfile, m_binaryPath);
    // Attention: when you use the model, you need to know the minIntenFC.
    //  it should be part of the model
    // we should save it into another file, that is the model

    SimpleTimer st("annotation");
    int scan = -1;
    string debug_peptide;
    if (m_debug) {
        cout << "Please type in the scan number you want to check" << endl;
        cin >> scan;
        cout << "Please type in the peptide you want to check" << endl;
        cin >> debug_peptide;
        cout << "debug peptide: " << debug_peptide << " debug scan: " << scan << endl;
    }

    // Annotation of mgfr with different peptide sequences.
    // From title (MGF) or Search-Results (PepXML)
    // Target/Positive: title, MGF
    // Decoy/Negative: Search-Results, PepXML
    string message;

    bool istraining = true;
    if (m_annotationMethod == "title") {
        mgfReader.annotate_with_title(verbosity, m_debug, scan, debug_peptide,
                                      istraining); // problem. annotation with what?
        message = "posfile " + m_feature_outputname;
    } else if (m_annotationMethod == "searchresult" && not m_searchresult.empty()) {
        int rank = CDebugMode::callDebug()->rank;
        mgfReader.annotation_with_search_result(m_searchresult, rank, m_debug, scan, 0, debug_peptide, istraining);
        message = "negafile " + m_feature_outputname;
    } else {
        spdlog::get("A")->info("Invalid annotation method: {}", m_annotationMethod);
        throw CException("invalid annotation method!");
    }
    notify(message);

    st.restart("feature export");
    mgfReader.export_features_with_PSM(features, m_feature_outputname, m_debug, scan, debug_peptide);

    // doing debug!!!
    if (CDebugMode::callDebug()->getMdebug() && scan != -1) {
        // todo: fix bug for invalid scan, e.g. -1
        st.restart("Ploting Spectrum");
        vector<tuple<double, double, string>> onespec;
        mgfReader.getSpectra(scan, onespec);
        Gnuplot gnuplot("gnuplot -persist ");
        gnuplot << "set yrange [0:12000];" << endl;
        gnuplot << "set xrange [0:2000];" << endl;
        gnuplot << "plot '-' using 1:2 with impulses, '-' using 1:2:3 with labels rotate left;" << endl;
        gnuplot.send1d(onespec);
        gnuplot.send1d(onespec);
    }
    cout << "[Info] Releasing resources..." << endl;
    for (auto x : features) {
        delete x;
    }
}

void FeatureWorkFlow::update_with_key_value_pair(string key, string value) {
    if (key == "fragmodel") {
        m_fragmodelfile = value;
    }
    if (key == "decoysearch") {
        m_searchresult = value;
    }
}

void FragModel::extractFeatures(bool verbosity, CPeakPairsImporter &positivePeakPairs,
                                CPeakPairsImporter &negativePeakPairs) {
    bool isTruthKnown = true;
    MGFReader mgfr(m_mgf_file, m_isLowMassAcc, isTruthKnown);
    bool doTrain = true;
    mgfr.annotate_with_title(verbosity, false, 0, "", doTrain);

    positivePeakPairs.getPeakPairSampleFrom(mgfr, m_use_ghost_peak);
    positivePeakPairs.setSourceFileName(mgfr.getMgfFileName());

    m_testdecoy = false;
    if (File::isExist(m_decoysearchfile)) {
        mgfr.annotation_with_search_result(m_decoysearchfile, 0, false, 0, 0, "", true);
        negativePeakPairs.getPeakPairSampleFrom(mgfr, m_use_ghost_peak);
        negativePeakPairs.setSourceFileName(m_decoysearchfile);
        m_testdecoy = true;
    }
}

void FragModel::update_with_key_value_pair(string key, string value) {
    if (key == "decoysearch") {
        m_decoysearchfile = value;
    }
}

FragModel::FragModel(const string &mgf_input, bool verbosity, double minIntenFC, bool overwrite,
                     string &fragPatternModel,
                     bool outputfeature, bool useGhostPeak, bool isLowMassAcc, string binaryPath) : CFlow() {
    m_isLowMassAcc = isLowMassAcc;
    m_name = "Train fragmentation model";
    m_mgf_file = mgf_input;
    m_use_ghost_peak = useGhostPeak;
    m_verbosity = verbosity;
    m_minIntenFC = minIntenFC;
    m_overwritemodel = overwrite;
    m_modelFileName = fragPatternModel;
    m_exportFeatureToFile = outputfeature;
    m_testdecoy = false;
    m_binaryPath = binaryPath;
    File::CFile filename(m_mgf_file);

    string tsearchprefix = filename.path + "/target_" + filename.basename,
            dsearchprefix = filename.path + "/decoy_" + filename.basename;

    m_decoysearchfile = dsearchprefix + ".pep.xml";
}


void FragModel::run() {
    spdlog::get("A")->info("Fragmentation model training workflow!");
    CPeakPairsImporter target_ds, decoy_ds;
    cout << "[Info] Loading fragment features:" << endl;
    extractFeatures(m_verbosity, target_ds, decoy_ds);
    CPeakPairsImporter target_SubData;
    cout << "[Info] refine fragment features: min intensity fold change greater than " << m_minIntenFC << endl;
    target_ds.getSubsetOnMinIntenFoldChange(target_SubData, m_minIntenFC);
    classifier lr(m_modelFileName, m_binaryPath);
    lr.setTestingMode(true);
    lr.setOutputMethod(m_exportFeatureToFile);


    lr.trainModel(target_SubData, m_overwritemodel);

    CFragScore *restarget = lr.predict(target_SubData);
    delete restarget;

    if (m_testdecoy) {
        CPeakPairsImporter decoy_subData;
        decoy_ds.getSubsetOnMinIntenFoldChange(decoy_subData, m_minIntenFC);
        CFragScore *resdecoy = lr.predict(decoy_subData);
        delete resdecoy;
    }

    notify(string("fragmodel ") + m_modelFileName);
}

FragModel::~FragModel() = default;

FragPatternScoreFlow::FragPatternScoreFlow(const string &inputfile, bool verbosity, double minIntenFC, bool overwrite,
                                           string outputmodelfile,
                                           string scoretype, const string &searchresult, bool isLowMassAcc,
                                           string binaryPath)
        : CFlow() {
    m_isLowMassAcc = isLowMassAcc;
    m_name = "Calculate fragmentation score";
    m_inputfile = inputfile;
    m_verbosity = verbosity;
    m_minIntenFC = minIntenFC;
    m_overwritemodel = overwrite;
    m_modelfile = outputmodelfile;
    m_scoretype = scoretype;
    m_searchresult = searchresult;
    m_binaryPath = binaryPath;

    // todo: make another function: check the search result
    // todo: why we need search result in workflow? We can remove it and use the annotation alone.
    if (searchresult.empty()) {
        cout << "[Error] Make sure provide the comet search result file; pep.xml" << endl;
        exit(0);
    }
}

void FragPatternScoreFlow::run() {
    spdlog::get("A")->info("Calculating PSM Frag-Score with fragmentation pattern!");
    CPeakPairsImporter fragmentation_sample;
    bool isTruthKnown = true;
    MGFReader mgfr(m_inputfile, m_isLowMassAcc, isTruthKnown);

    mgfr.annotation_with_search_result(m_searchresult, 0, false, 0, 0, "");
    fragmentation_sample.getPeakPairSampleFrom(mgfr, true);
//    cout << "Line " << __LINE__ << "\t" << __FILE__ << endl;
    fragmentation_sample.setSourceFileName(m_searchresult);
    cout << "[Info] Loading fragment features:" << endl;
    CPeakPairsImporter fragmentation_subsample;
    fragmentation_sample.getSubsetOnMinIntenFoldChange(fragmentation_subsample, m_minIntenFC);
    classifier lr(m_modelfile, m_binaryPath);

    FragmentationScore Y(fragmentation_subsample, m_scoretype, lr);
    string filename = Y.exportScore(m_scoretype);
    cout << "PSM Frag-Score saved to : " << filename << endl;
}

FragPatternScoreFlow::~FragPatternScoreFlow() = default;

CometSearch::CometSearch(const string &cometbinary, const string &database, const string &paramfile,
                         const string &outputname, const string &specfile)
        : CFlow() {
    m_name = "Comet search";
    m_comet_binary = cometbinary;
    m_paramfile = paramfile;
    m_database = database;
    m_outputname = outputname;
    m_specfile = specfile;
}

void CometSearch::run() {
    spdlog::get("A")->info("Comet search workflow started!");
    if (not File::isExist(m_database)) {
        spdlog::get("A")->error("Database {} does not exist!", m_database);
        throw CException("database file does not exist!");
    }
    if (not File::isExist(m_comet_binary)) {
        spdlog::get("A")->error("Comet binaryfile {} does not exist!", m_comet_binary);
        throw CException("Comet binary file does not exist!");
    }
    if (not File::isExist(m_paramfile)) {
        spdlog::get("A")->error("Comet paremeter file {} does not exist!", m_paramfile);

        throw CException("Comet paremter file does ot exist!");
    }
    if (m_outputname.empty()) {
        spdlog::get("A")->error("Comet output name is empty!");

        throw CException("Comet outputname should not be empty!");

    }
    string cmdline =
            m_comet_binary + " -P" + m_paramfile + " -N" + m_outputname + " -D" + m_database + " " + m_specfile;

    spdlog::get("A")->info("Running: {0}", cmdline);
    int retval = system(cmdline.c_str());
    spdlog::get("A")->info("Comet exit with code: {0}", retval);
    if (retval != 0) {
        throw CException("Invalid argument");
    }

    string datapth, specfilename;
    File::parent(m_specfile, datapth, specfilename);
    m_outputname += ".pep.xml";
    if (File::isExist(m_outputname)) {
        cout << "output to " << m_outputname << endl;
    }
    notify(m_database + " " + m_outputname);
}

CometSearch::~CometSearch() = default;

CometSearchTD::CometSearchTD(string cometbinary, string targetdb, string decoydb, string paramfile, string specfile)
        : CFlow() {
    m_name = "Comet target-decoy separated search";
    m_comet_binary = cometbinary;
    m_paramfile = paramfile;
    m_targetdb = targetdb;
    m_decoydb = decoydb;
    m_specfile = specfile;
    File::CFile filename(m_specfile);
    string target_outname = filename.path + "/target_" + filename.basename,
            decoy_outname = filename.path + "/decoy_" + filename.basename;

    m_pTarget = new CometSearch(m_comet_binary, m_targetdb, m_paramfile, target_outname, m_specfile);
    m_pDecoy = new CometSearch(m_comet_binary, m_decoydb, m_paramfile, decoy_outname, m_specfile);
    this->subscribe(m_pTarget);
    this->subscribe(m_pDecoy);
}

CometSearchTD::~CometSearchTD() {
    delete m_pTarget;
    delete m_pDecoy;
}

void CometSearchTD::run() {
    spdlog::get("A")->info("Target decoy separated search workflow started!");

    if (not File::isExist(m_targetdb)) {
        spdlog::get("A")->error("Target-only database \"{}\" does not exist!", m_targetdb);
        throw CException("Target-only database does not exist!");
    }
    if (not File::isExist(m_decoydb)) {
        spdlog::get("A")->error("Decoy-only database \"{}\" does not exist!", m_decoydb);
        throw CException("Decoy-only database does not exist!");
    }

    if (m_pTarget == nullptr or m_pDecoy == nullptr) {
        cout << "[Error] Fail to create target decoy task" << endl;
        throw CException("Fail to initialize target and decoy search tasks!");
    }
    m_pTarget->run();
    m_pDecoy->run();
}

void CometSearchTD::update_with_key_value_pair(string key, string value) {
    if (key == m_targetdb) {
        notify(string("targetsearch ") + value);
    }
    if (key == m_decoydb) {
        notify(string("decoysearch ") + value);
    }
}

// todo: rangerbinary changed
FlowAll::FlowAll(const string &fragscoretype, bool useGhostPeak, bool outputfeature, bool overwrite, double minIntenFC,
                 bool verbosity,
                 string fragPatternModelFilename, const string &inputfile, string cometbinary, string targetdb,
                 string decoydb,
                 string paramfile, bool isTrainingRF, bool isRFProbOut, string writeRFModelTo, string validatorRFmodel,
                 int mtry, int ntree, const string &featurelistfile, bool isLowMassAcc, int trainingSampleSize,
                 int maxdepth,
                 const string &rangerbinary, const string &binaryPath) : CFlow() {
    if (not File::isExist(inputfile)) {
        spdlog::get("A")->error("Input file does not exist!");
        throw CException("Invalid input file!");
    }
    m_name = "All tasks";
    m_inputspec = inputfile; // todo: should be mgf now. in the future, we will support mzXML
    File::CFile filename(m_inputspec);
    string tsearch = filename.path + "/target_" + filename.basename,
            dsearch = filename.path + "/decoy_" + filename.basename;

    // if target search or decoy search result does not exist
    if (not File::isExist(tsearch + ".pep.xml") or not File::isExist(dsearch + ".pep.xml")) {
        cout << "[Info] target-decoy separated search not complete" << endl;
        m_pFlow.push_back(new CometSearchTD(cometbinary, targetdb, decoydb, paramfile, inputfile));
    }

    if (not File::isExist(fragPatternModelFilename, true)) {
//        spdlog::get("A")->info("Fragmentation model \"{}\" does not exist!", fragPatternModelFilename);

        m_pFlow.push_back(new FragModel(inputfile, verbosity, minIntenFC, overwrite, fragPatternModelFilename,
                                        outputfeature, useGhostPeak, isLowMassAcc, binaryPath));
    }
    string featurelistnameBase = File::CFile(featurelistfile).basename;
    string psmfeaturefile = to_string(inputfile, "_",
                                      featurelistnameBase, ".feat");
    string decoypsmfeaturefile = to_string(dsearch + ".pep.xml",
                                           "_", featurelistnameBase, ".feat");

    // two different feat files generated here.
    bool debug = CDebugMode::callDebug()->getMdebug();
    if (not File::isExist(psmfeaturefile, true) or overwrite) {
//        cout << "[Info] PSM feature file \"" << psmfeaturefile << "\" does not exist!" << endl;
        m_pFlow.push_back(new FeatureWorkFlow(inputfile, useGhostPeak, minIntenFC,
                                              fragPatternModelFilename, fragscoretype,
                                              "title", debug, "",
                                              featurelistfile, isLowMassAcc, binaryPath));

    }
    if (not File::isExist(decoypsmfeaturefile, true) or overwrite) {
//        cout << "[Info] PSM feature file \"" << decoypsmfeaturefile << "\" does not exist!" << endl;
        m_pFlow.push_back(new FeatureWorkFlow(inputfile, useGhostPeak, minIntenFC,
                                              fragPatternModelFilename, fragscoretype,
                                              "searchresult", debug, "",
                                              featurelistfile, isLowMassAcc, binaryPath));

    }

    // post-processing using: alglib or mlpack.
    // todo: put the last step into python
    string combinedPNfeature = to_string(inputfile, "_", featurelistnameBase, "trN", trainingSampleSize,
                                         "_PNcombined.feat");
//            inputfile + "_PNcombined.feat";

    if (not File::isExist(combinedPNfeature, true) or overwrite) {
        m_pFlow.push_back(new RangerFormat(psmfeaturefile, decoypsmfeaturefile, combinedPNfeature, trainingSampleSize));

    }
    if (isTrainingRF) {
        RFModelConfig rfconfig(combinedPNfeature, isTrainingRF, writeRFModelTo, true, mtry, ntree, maxdepth,
                               rangerbinary, writeRFModelTo + ".conf");
        rfconfig.write();
        m_pFlow.push_back(new RangerWraper(rfconfig));

    } else {
        RFModelConfig rfconfig(combinedPNfeature, isTrainingRF, validatorRFmodel, true, mtry, ntree, maxdepth,
                               rangerbinary, validatorRFmodel + ".conf");
        rfconfig.write();
        m_pFlow.push_back(new RangerWraper(rfconfig));

    }

    // build notification chain
    for (int i = 1; i < m_pFlow.size(); i++) {
        for (int j = 0; j < i; j++) {
            m_pFlow[i]->subscribe(m_pFlow[j]);
        }
    }
}

void FlowAll::run() {
    spdlog::get("A")->info("Start to run combined workflow!");
    spdlog::get("A")->info("There are {0} tasks in the queue!", m_pFlow.size());
    for (int i = 0; i < m_pFlow.size(); i++) {
        if (m_pFlow[i] != nullptr) {
            spdlog::get("A")->info("Task {}: {}", i + 1, m_pFlow[i]->getname());
        } else {
            spdlog::get("A")->error("NULL pointer detected for {}-th task", i + 1);
            cerr << __func__ << " -> " << __LINE__ << endl;
            throw CException("Pointer to workflow should not be NULL!");
        }
    }
    cout << endl;

    for (int i = 0; i < m_pFlow.size(); i++) {
        if (m_pFlow[i] != nullptr) {
            spdlog::get("A")->info("Running task {}: {}", i + 1, m_pFlow[i]->getname());
            m_pFlow[i]->run();
        } else {
            spdlog::get("A")->error("NULL pointer detected for {}-th task", i + 1);
            cerr << __func__ << " -> " << __LINE__ << endl;
            throw CException("Pointer to workflow should not be NULL!");
        }
    }
}

FlowAll::~FlowAll() {
    for (auto &i : m_pFlow) {
        if (i != nullptr) {
            delete i;
        }
    }
}

void scoreToROC(vector<double> &postiveScores, vector<double> &negativeScores, vector<tuple<double, double>> &roc) {
    int postive_num = postiveScores.size(), negative_num = negativeScores.size();
    sort(postiveScores.begin(), postiveScores.end());
    sort(negativeScores.begin(), negativeScores.end());

    double pcount = 0, ncount = 0;
    roc.emplace_back(0, 0);
    while (pcount < postive_num and ncount < negative_num) {
        if (postiveScores[pcount] < negativeScores[ncount]) {
            pcount++;
        } else if (postiveScores[pcount] > negativeScores[ncount]) {
            ncount++;
            double TPR = pcount / postive_num, FPR = ncount / negative_num;

            roc.emplace_back(FPR, TPR);
        } else {
            pcount++;
            ncount++;
            double TPR = pcount / postive_num, FPR = ncount / negative_num;
            roc.emplace_back(FPR, TPR);
        }
    }
    while (pcount < postive_num) {
        pcount++;
        double TPR = pcount / postive_num, FPR = 1.0;
        roc.emplace_back(FPR, TPR);
    }

    while (ncount < negative_num) {
        ncount++;
        double TPR = 1.0, FPR = ncount / negative_num;
        roc.emplace_back(FPR, TPR);
    }
}

void saveROCtoTXT(const string &roc_data, vector<tuple<double, double>> &roc) {
    ofstream fout(roc_data.c_str(), ios::out);
    for (auto x: roc) {
        fout << get<0>(x) << " " << get<1>(x) << endl;
    }
    fout.close();
}

double getAUCfromROC(vector<tuple<double, double>> &roc) {
    double auc = 0;
    for (int i = 1; i < roc.size(); i++) {
        auc += (get<0>(roc[i]) - get<0>(roc[i - 1])) * (get<1>(roc[i]) + get<1>(roc[i - 1])) / 2.0;
    }
    return auc;
}


// Assumption of this function:
// the first half of the data is positive, the second half is negative
// The assumption can be invalid sometimes...
void getscore(vector<vector<string>> &data, vector<double> &postiveScores, vector<double> &negativeScores) {
    int num_skip_header_lines = 3;
    int sample_num = (int) data.size() - num_skip_header_lines;
    cout << "Sample Number: " << sample_num << endl;
    int postive_num = sample_num / 2, negative_num = sample_num / 2;
    if (sample_num % 2 == 1) {
        cout << "[Warning]: For testing, the sample number should be even, Not odd number???" << endl;
        return;
    }

    for (int i = num_skip_header_lines; i < data.size(); i++) {
        vector<string> &entry = data[i];
        if (i < postive_num + num_skip_header_lines) {
            postiveScores.push_back(atof(entry[1].c_str()));
        } else {
            negativeScores.push_back(atof(entry[1].c_str()));
        }
    }
}

void
plotScoreDistribution(vector<double> &postiveScores, vector<double> &negativeScores,
                      const string &scoreDistributionPlotFile,
                      const string &xlabel) {
//          "set boxwidth 0.01 absolute\n"
    Gnuplot gp;
    gp << "set terminal pngcairo enhanced font \"Arial,18\" size 540,480 " << endl;
    gp << "set output '" << scoreDistributionPlotFile << "'" << endl;
    gp << "set xrange [0:1]\n"
          "# Each bar is half the (visual) width of its x-range.\n"
          "set style fill transparent solid 0.5\n"
          "set xlabel '" << xlabel << "'\n"
                                      "set ylabel 'Frequency'\n"
                                      "set style fill transparent solid 0.5\n"
                                      "bin_width = 0.1;\n"
                                      "\n"
                                      "bin_number(x) = floor(x/bin_width)\n"
                                      "set boxwidth 0.9*bin_width\n"
                                      "rounded(x) = bin_width * ( bin_number(x) + 0.5 )\n"
                                      "\n" << endl;

    gp << "plot '-' using (rounded($1)):(100.0/" << postiveScores.size()
       << ") smooth frequency with boxes title 'positive', '-' using (rounded($1)):(100.0/" << negativeScores.size()
       << ") smooth frequency with boxes title 'negative'" << endl;
    gp.send1d(postiveScores);

    gp.send1d(negativeScores);
}


void plotImportanceBarChart(const string &importance_png_file, vector<tuple<string, double>> &Importance) {
    Gnuplot gp;
    gp << "set terminal png noenhanced font arial 18 size 3072,1024 " << endl;
//        gp << "set samples 10240" << endl;
    gp << "set output '" << importance_png_file << "'" << endl;
    gp << "set bmargin 15" << endl;
    gp << "set rmargin 10" << endl;
    gp << "set lmargin 15" << endl;
    gp << "set xtics right offset 0,0" << endl;

    gp << "set xtics ( ";
    for (int i = 0; i < Importance.size(); i++) {
        gp << " \"" << get<0>(Importance[i]) << "\" " << i;
        if (i < Importance.size() - 1) {
            gp << ", ";
        }
    }
    gp << ")" << endl;
//        gp << "set xtics center offset 0,-1" << endl;
    gp << "set xtics rotate by 45  " << endl;
//        gp << "set ytics font \"Arial,14\"" << endl;
    gp << "set border 1+2 back\n"
          "set tics scale 0"
       << endl;
    gp << "set style fill transparent solid 0.5 " << endl;

//        gp << "unset key" << endl;
    gp << "plot '-' using 1 with boxes lc rgb\"red\" notitle" << endl;
    vector<double> importancevalue;
    for (auto x: Importance) {
        importancevalue.push_back(get<1>(x));
    }
    gp.send1d(importancevalue);
}

void FeatureImportanceProcess(const string &importance_file) {
    vector<vector<string>> data;
    vector<tuple<string, double>> Importance;

    load_tsv_file(importance_file, data);
    while (data[data.size() - 1].empty()) { data.pop_back(); }

    for (auto x: data) {
        Importance.emplace_back(x[0].substr(0, x[0].length() - 1), atof(x[1].c_str()));
    }
    //todo:  tobe fixed in the future.
    // we could not pass the file name to stream.
    string importance_png_file = importance_file + ".png";

    plotImportanceBarChart(importance_png_file, Importance);
    string sorted_importance_png_file = importance_file + "_sorted.png";
    sort(Importance.begin(), Importance.end(),
         [](const tuple<string, double> &x, const tuple<string, double> &y) { return get<1>(x) > get<1>(y); });
    plotImportanceBarChart(sorted_importance_png_file, Importance);
}

void RangerWraper::run() {
    if (m_isTraining) {
        train();
        // the input file should be be hard coded here
        string importance_file = m_rfModel.substr(0, m_rfModel.length() - 7) + ".importance";
        FeatureImportanceProcess(importance_file);
        notify(string("RFModel ") + m_rfModel);
    } else {
        predict();

    }
}

RangerWraper::RangerWraper(RFModelConfig &rfconfig) {
    m_name = "Run ranger (random-forest)";
    m_tsvfile = rfconfig.getMTsvfile();
    m_isTraining = rfconfig.isMIsTraining();
    m_rfModel = rfconfig.getMRfModel();
    m_mtry = rfconfig.getMMtry();
    m_ntree = rfconfig.getMNtree();
    m_maxDepth = rfconfig.getMMaxDepth();


    m_probPrediction = rfconfig.isMProbPrediction();
    m_predictionprefix = m_tsvfile.substr(0, m_tsvfile.length() - 5); // magic value
    m_ranger_binary_path = rfconfig.getMRangerBinaryPath();
    m_threadnum = rfconfig.getMThreadnum();
}

void RangerWraper::predict() const {
    string ranger_commandline =
            m_ranger_binary_path + " --treetype 1 --ntree " + to_string(m_ntree) + "  --verbose --depvarname class  "
                                                                                   " --targetpartitionsize 1 --impmeasure 1 --outprefix " +
            m_predictionprefix +
            " --nthreads " + to_string(m_threadnum) + " --seed 100 --predict " + m_rfModel + " --mtry " +
            to_string(m_mtry);
    if (m_probPrediction) {
        ranger_commandline = ranger_commandline + " --probability ";
    }
    ranger_commandline = ranger_commandline + " --file " + m_tsvfile;

    spdlog::get("A")->info("Running: {0}", ranger_commandline);
    int ret = system(ranger_commandline.c_str());
    spdlog::get("A")->info("Ranger exits with code {0}", ret);

    if (ret != 0) {
        spdlog::get("A")->info("Fail to predict score with random forest model! Program will exit!");
        exit(0);
    }
}

void RangerWraper::train() const {
    string ranger_commandlinex = to_string(m_ranger_binary_path, " ",
                                           " ", "--treetype", 1, "--ntree", m_ntree,
                                           "--maxdepth", m_maxDepth,
                                           "--verbose --depvarname class --targetpartitionsize 1 --write",
                                           "--outprefix", m_rfModel.substr(0, m_rfModel.length() - 7),
                                           "--impmeasure", 1,
                                           "--nthreads", m_threadnum,
                                           "--seed", 100,
                                           "--mtry", m_mtry);

    string ranger_commandline = m_ranger_binary_path + " --treetype 1 --ntree " + to_string(m_ntree)
                                + "  --maxdepth " + to_string(m_maxDepth) +
                                " --verbose --depvarname class --targetpartitionsize 1 --write --outprefix "
                                + m_rfModel.substr(0, m_rfModel.length() - 7)
                                + " --impmeasure 1  --nthreads " + to_string(m_threadnum) + " --seed 100 --mtry " +
                                to_string(m_mtry);

    cout << ranger_commandline << endl << ranger_commandlinex << endl;
    if (m_probPrediction) {
        ranger_commandline = ranger_commandline + " --probability ";
    }
    ranger_commandline = ranger_commandline + " --file " + m_tsvfile;

    spdlog::get("A")->info("Running: {0}", ranger_commandline);
    int ret = system(ranger_commandline.c_str());
    cout << "[Ranger] exit with return code: " << ret << endl;
    if (ret != 0) {
        cout << "Program will exit since it fails to train the random forest model" << endl;
        exit(0);
    }
}

void RangerWraper::update_with_key_value_pair(string key, string value) {
    if (key == "rangerfeature") {
        m_tsvfile = value;
        m_predictionprefix = m_tsvfile.substr(0, m_tsvfile.length() - 5); // remove ext
    }
    if (key == "RFModel") {
        m_rfModel = value;
    }
}


void load_tsv_file(const string &filename, vector<vector<string>> &data) {
    string line;
    ifstream fin(filename.c_str(), ios::in);
    if (fin.is_open()) {
        while (getline(fin, line)) {
            vector<string> tokens;

            trim_space_only(line);
            split_string(line, tokens);
            if (tokens.size() == 1) {
                cout << "skip one token line: " << line << endl;
                // todo: this is unexpected behavior!
                continue;
            }
            data.push_back(tokens);
        }
        fin.close();
    }
    cout << "[Info] " << data.size() << " lines loaded from " << filename << endl;
}

string getfragmodelfile(const string &validatorModel, double m_minIntenFC, bool isTraining) {
    File::CFile filename(validatorModel);
    if (not filename.isFileExist(true) and not isTraining) {
        cout << "RF model file does not exist: " << validatorModel << endl;
        throw CException("Run time Error! File does not exist!");
    }
    string basename = filename.basename;
    int found = basename.find("mtry_");
    if (found != string::npos) {
        basename = basename.substr(0, found);
    }
    cout << "basename change from " << filename.basename << " to " << basename << endl;
    string fragModelFileName = filename.path + "/" + basename + "_" + to_string(m_minIntenFC) + "_frag.model";
    cout << "FragModelFilename is : " << fragModelFileName << endl;
    return fragModelFileName;
}

#include <unordered_set>

vector<Feature *>
createFeatureVector(const string &fragmodel, bool ghost, double minIntFC, const string &fragscoretype,
                    const string &featurelistfile,
                    string binaryPath) {
    vector<Feature *> features, tmp;
    spdlog::get("A")->debug("Create Feature Vector param: "
                            "ghost={}; fragmodelfile={}; scoretype={}, minIntFC={}", ghost, fragmodel, fragscoretype,
                            minIntFC);
    vector<string> featurenames = readlines(featurelistfile);
    unordered_set<string> featurenameset(featurenames.begin(), featurenames.end());

    features.push_back(new Peptide_Length);
    features.push_back(new Precursor_Charge);
    features.push_back(new Num_of_Peaks);

    // fragment ion mass error
    features.push_back(new Mass_Error_Std);
    features.push_back(new Mass_Error_AbsAvg);

    // Matched Peaks -------------------------------
    // Matched Peaks: Counts, b, y and b&y
    // Numbers, not Normalized
    features.push_back(new Num_Of_b_ions);
    features.push_back(new Num_Of_y_ions);
    features.push_back(new Num_Of_Matched_Peaks);

    // Ratio, Normalized
    features.push_back(new Num_Of_b_ions_Prop);
    features.push_back(new Num_Of_y_ions_Prop);
    features.push_back(new Num_Of_Matched_Peaks_Prop);
    features.push_back(new Num_Of_Matched_Peaks_Prop_Calibrated);

    // Match Intensity: ------------------
    // absoluten (NOT relative) values Intensity b, y and b&y
    features.push_back(new Matched_b_Intensity);
    features.push_back(new Matched_y_Intensity);
    features.push_back(new Matched_b_y_Intensity);
    // relative ratio,devided by total intensity
    features.push_back(new Matched_b_Intensity_Prop);
    features.push_back(new Matched_y_Intensity_Prop);
    features.push_back(new Matched_b_y_Intensity_Prop);


    // Amino Acid counts ------------------------------------
    string legalAAs = "ACDEFGHIKLMNPQRSTVWY";
    for (char legalAA : legalAAs) {
        features.push_back(new Num_of_AA(legalAA));
    }
    // Charge status -----------------------------------------
    vector<int> legalChargeStates = {1, 2, 3, 4, 5, 6, 7};
    for (int legalChargeState : legalChargeStates) {
        features.push_back(new isChargeX(legalChargeState));
    }

    // Unassigned Intensity -----------------------------------
    // spectrum side: peak intensity
    // absolute value
    features.push_back(new Unassigned_Intensity_Prop);
    // relative value
    features.push_back(new UnassignedIntensity);

    // total intensity
    features.push_back(new Total_Intensity);                // NOT Promising

    // Explained Cleavages------------------------------------------
    // Explained peptide bonds by b, y or b&y at any charge state
    // Not Normalized ones
    features.push_back(new Explained_Cleavage);             //
    features.push_back(new Explained_b_Cleavage);
    features.push_back(new Explained_y_Cleavage);
    // Ratios Normalized ones
    features.push_back(new Explained_Cleavage_Prop);
    features.push_back(new Explained_b_Cleavage_Prop);
    features.push_back(new Explained_y_Cleavage_Prop);

    // Explained Cleavages, with Only 1+ ions --------------------------
    // Explained peptide bonds by specific ion types: b+, y+ b+&y+
    // Not normalized
    features.push_back(new Explained_cs1_Cleavage);             //
    features.push_back(new Explained_b_cs1_Cleavage);
    features.push_back(new Explained_y_cs1_Cleavage);
    // Ratios Normailized
    features.push_back(new Explained_cs1_Cleavage_Prop);             //
    features.push_back(new Explained_b_cs1_Cleavage_Prop);
    features.push_back(new Explained_y_cs1_Cleavage_Prop);

    // Longest ion series --------------------------------------------------------
    // Longest Ion series calculated with b, y or b&y ion types at any charge state
    // Numbers Not normalized
    features.push_back(new Longest_ion_series);
    features.push_back(new Longest_b_ion_series);
    features.push_back(new Longest_y_ion_series);
    //Ratios normalized
    features.push_back(new Longest_ion_series_Prop);
    features.push_back(new Longest_b_ion_series_Prop);
    features.push_back(new Longest_y_ion_sereis_Prop);

    // Longest ion series, with only 1+ ions ---------------------------------
    // Longest Ion series calculated with b, y or b&y when charge=1
    // Numbers Not normalized
    features.push_back(new Longest_cs1_ion_series);
    features.push_back(new Longest_b_cs1_ion_series);
    features.push_back(new Longest_y_cs1_ion_series);
    // Ratios, normalized
    features.push_back(new Longest_cs1_ion_series_Prop);
    features.push_back(new Longest_b_cs1_ion_series_Prop);
    features.push_back(new Longest_y_cs1_ion_series_Prop);

    // Fragmentation pattern score ------------------------------------------------
    features.push_back(new FragmentationPatternScore(ghost, minIntFC, fragmodel, fragscoretype, binaryPath));

    cout << "Num of Features used : " << features.size() << endl;
    tmp.swap(features);
    for (int i = 0; i < tmp.size(); i++) {
        bool keepFeature = false;
        if (featurenameset.count(tmp[i]->m_feature_name)) {
            keepFeature = true;
            features.push_back(tmp[i]);
        }
        cout << "Feature #" << i + 1 << ": " << tmp[i]->m_feature_name << ": " << std::boolalpha << keepFeature << endl;
        if (not keepFeature) {
            delete tmp[i];
            tmp[i] = nullptr;
        }
    }
    cout << "--------end-----------" << endl;
    return features;
}

// The FDR
// Number of Identified PSM and FDR
// ACS: 3.25 in
// ACS: 600dpi
void plot_FDR_curve(const string &outputfilename, const string &title, vector<tuple<double, double>> &fdr_counts) {
    int ACS_DPI = 600;
    int defaultDPI=72;
    double DPI_fold = ACS_DPI / (double) defaultDPI;
    double ACS_single_img_width = 3.25; // inch
    double width = ACS_single_img_width * ACS_DPI;
    double height = ACS_single_img_width * ACS_DPI;
    // font size measured in points. 1 point is 1/72 inch.
    // image of 3.25 inch height, with 12 point size, corresponds to :
    // font is X, width/72 is the points number. width/600 is the previous points number.
    int newFontSize = int(12*DPI_fold); //
//    cout << "font size is : " << newFontSize << endl;
    Gnuplot gnuplot;
    gnuplot << "dpi = 300 ## dpi (variable)\n"
               "width = 90 ## mm (variable)\n"
               "height = 90 ## mm (variable)\n"
               "\n"
               "in2mm = 25.4 # mm (fixed)\n"
               "pt2mm = 0.3528 # mm (fixed)\n"
               "\n"
               "mm2px = dpi/in2mm\n"
               "ptscale = pt2mm*mm2px\n"
               "round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)\n"
               "wpx = round(width * mm2px)\n"
               "hpx = round(height * mm2px)\n"
               "#set term pngcairo size 720,720 fontscale scale linewidthscale scale\n"
               "\n"
               "set terminal pngcairo size wpx,hpx font 'Vadana,10' fontscale ptscale linewidth ptscale pointscale ptscale "
               << endl;
    //gnuplot << "set terminal pngcairo enhanced font \"Arial "<< newFontSize<<"\" size " << width << "," << height <<  endl;
    gnuplot << "set output '" << outputfilename << "'" << endl;
    gnuplot << "set size square " << endl;
    gnuplot << "set xlabel 'FDR (%)' font ',8'" << endl;
    gnuplot << "set ylabel 'No of Identified PSMs (K) ' font ',8'" << endl;
    gnuplot << "set xrange[0:5]" << endl; // 5%
    gnuplot << "set tics font ',8'" << endl;
    gnuplot << "set key bottom right" << endl;
    gnuplot << "n=1000" << endl;
    gnuplot << "m=100" << endl; // percentages.
    gnuplot << "plot '-' using ($1*m):($2/n) with line lw 2 " <<  " linecolor rgb \"black\" " << " title \"" << title << "\"" << endl;

    gnuplot.send1d(fdr_counts);
}


void
get_feature_from_spec(PSMInfo *psmInfo, CSpectrum *spec, vector<double> *OneSpecFeature, vector<Feature *> *features,
                      int idx, string mzMLsourcefile, const bool fixMz, bool highMassAcc, int hitrank = 0) {

    SpectraSTPeakList *specpkl = new SpectraSTPeakList(1.0, 1);
    string fragType = highMassAcc ? "HCD" : "CID";

    specpkl->setFragType(fragType);
    if (!specpkl) {
        spdlog::get("A")->error("Get nullptr in feature extraction step ");
        throw CException("Invalid pointer");
    }
    specpkl->setParentMz(psmInfo->getParentMz());
    specpkl->setParentCharge(psmInfo->charge);

    for (int pkidx = 0; pkidx < spec->getPeakNum(); pkidx++) {
        double mz, intensity;
        spec->getOnePeak(mz, intensity, pkidx);
        spdlog::get("A")->trace("peak {} : {} \t {}", pkidx, mz, intensity);
        specpkl->insert(mz, intensity, "", "");
    }

    if (psmInfo->searchhits.size() <= hitrank) {
//        cout << "the " << hitrank << "-th peptide is not found: total candiates: " << psmInfo->searchhits.size() << endl;
        delete specpkl;
        return;
    }
    string modified_peptide = psmInfo->searchhits[hitrank]->m_modified_peptide;
    string plain_peptide = regex_replace(modified_peptide, regex("[[0-9.]*]"), "");
    string newPep = regex_replace(modified_peptide, regex("C"), "C[160]");

    spdlog::get("A")->debug("Modified peptide to be used for annotation: {}",
                            psmInfo->searchhits[0]->m_modified_peptide);
    Peptide *ipep = new Peptide(newPep, psmInfo->charge); // peptide done
    if (!ipep) {
        spdlog::get("A")->error("Get nullptr in feature extraction step ");
        throw CException("Invalid pointer");
    }
    specpkl->setPeptidePtr(ipep);

    // If set to true, it will fix the Mz to the theoretical value of the first annotated peak.
    // tobe continued!
    specpkl->annotate(true, fixMz);
    specpkl->normalizeTo(10000);

    PSMFeature *pPSM = new PSMFeature(specpkl, ipep, idx, mzMLsourcefile);
    if (!pPSM) {
        spdlog::get("A")->error("Get nullptr in feature extraction step ");
        throw CException("Invalid pointer");
    }

    for (int i = 0; i < features->size(); i++) {
        Feature *f = features->at(i);
        if (f == nullptr) {
            throw CException(("feature pointer is null for i=" + to_string(i)).c_str());
        }
        (*OneSpecFeature)[i] = f->calculate_feature(pPSM);
    }

    delete specpkl;
    delete ipep;
    delete pPSM;
}


int RangerFormat::exportTrainingFeaturesToCSV(vector<vector<string>> &feature, const string &label, ofstream &fout,
                                              int start_col,
                                              int start_row, const int MAX_SAMPLE_NUM) const {
    int seed = 42;
    int len = feature.size() - 1;
    vector<int> idx(len, 0);
    iota(idx.begin(), idx.end(), start_row);
    gcc5shuffle(idx.begin(), idx.end(), mt19937(seed));
    int sample_counts = 0;
    for (int i = 0; i < idx.size() and i < MAX_SAMPLE_NUM; i++) {
        vector<string> &entry = feature.at(idx[i]);
        for (int j = start_col; j < entry.size(); j++) {
            string &item = entry.at(j);
            fout << item << ",";
        }
        fout << label << endl;
        sample_counts++;
    }

    return sample_counts;
}

void RangerFormat::run() {
    vector<vector<string>> posdata, negdata;
    const string positive_label = "1";
    const string negative_label = "0";
    load_tsv_file(m_positivefile, posdata);
    load_tsv_file(m_negativefile, negdata);
    spdlog::get("A")->info("Training samples loaded! Positive: {}, Negative: {}",
                           posdata.size(), negdata.size());

    if (posdata.size() > m_MIN_SAMPLE_NUM and negdata.size() > m_MIN_SAMPLE_NUM) {
        ofstream fout(m_tsvfile.c_str(), ios::out);
        if (fout.is_open()) {
            const int HEADER_ROW_NUM = 0, COL_NUM = posdata[HEADER_ROW_NUM].size();
            int start_col = 1;  // skip first column, as it is PEPTIDE
            // output header...
            cout << "Exporting header" << endl;
            for (int j = start_col; j < COL_NUM; j++) {
                fout << posdata[HEADER_ROW_NUM][j] << ",";
            }
            fout << "class" << endl;
            // output content
            int start_row = 1;  // skip the first row, as it is headerline!
            cout << "exporting content" << endl;
            int pos_out_num = exportTrainingFeaturesToCSV(posdata, positive_label, fout, start_col, start_row,
                                                          m_MAX_SAMPLE_NUM);
            int neg_out_num = exportTrainingFeaturesToCSV(negdata, negative_label, fout, start_col, start_row,
                                                          m_MAX_SAMPLE_NUM);
            fout.close();
            spdlog::get("A")->info("{} Samples (POS: {}, NEG: {}) exported in ranger's csv format to {}",
                                   pos_out_num + neg_out_num, pos_out_num, neg_out_num, m_tsvfile);
            notify(string("rangerfeature ") + m_tsvfile);
        } else {
            spdlog::get("A")->info("Fail to create/open file in ranger's csv format: {0} ", m_tsvfile);
            throw CException("Fail to create csv file for RF training!");
        }
    } else {
        spdlog::get("A")->info("Minimum number of samples required {}", m_MIN_SAMPLE_NUM);
        throw CException("Too few sample to train RF model!");
    }
}


RangerFormat::RangerFormat(const string &posFeature, const string &negFeature, const string &outPNsampleName,
                           int trainingSampleSize) :
        m_MIN_SAMPLE_NUM(10), m_MAX_SAMPLE_NUM(trainingSampleSize) {
    m_name = "Convert positive/negative feature into csv file";
    m_positivefile = posFeature;
    m_negativefile = negFeature;
    m_tsvfile = outPNsampleName;

}

void RangerFormat::update_with_key_value_pair(string key, string value) {
    if (key == "posfile") {
        m_positivefile = value;
    }
    if (key == "negafile") {
        m_negativefile = value;
    }
}


ValidatePSM::ValidatePSM(string pepxmlfile, const string &validatorModel, string fragmodelscoretype, double minInt,
                         bool ghost,
                         int threadNum, shared_ptr<CTruth> truth, int mtry, int ntree, string featurelistfile,
                         int maxdepth,
                         bool useAlternativeProt, string rangerbinary, string binaryPath) {
    m_binaryPath = binaryPath;
    m_useAlternativeProtein = useAlternativeProt;
    m_name = "Feature extraction";
    fixMz = false;
    m_highMassAcc = true;
    m_featurelistfile = featurelistfile;
    m_validatorModelFile = validatorModel;
    m_pepxmlfile = pepxmlfile;
    m_fragmodelfile = getfragmodelfile(validatorModel, minInt);;
    m_fragmodelscoretype = fragmodelscoretype;
    m_minIntFC = minInt;
    m_ghost = ghost;
    m_basename = File::CFile(validatorModel).basename;

    m_output_feature_file = m_pepxmlfile + m_basename + ".feat";

    RFModelConfig rfconfig(m_output_feature_file, false, validatorModel, true, mtry, ntree,
                           maxdepth, rangerbinary, validatorModel + ".conf");
    rfconfig.write();
    m_other_steps.push_back(new RangerWraper(rfconfig));
    m_other_steps.push_back(
            new resultAnalysis(m_pepxmlfile, truth, File::CFile(validatorModel).basename, m_useAlternativeProtein));
    m_threadNum = threadNum <= 0 or threadNum > getProperThreads() ? getProperThreads() : threadNum;
}

void ValidatePSM::getProcessed_i(DataFile *df, PeptideProphetParser &ppp, vector<int> &index) {
    string debug_peptide;
    int debug_scan = -1;
    if (CDebugMode::callDebug()->getMdebug()) {
        cout << "Please type in debug peptide : " << endl;
        cin >> debug_peptide;
        cout << "Debug peptide is : " << debug_peptide << endl;
        cout << "Please type in debug scan: " << endl;
        cin >> debug_scan;
        cout << "Debug scan num is : " << debug_scan << endl;
        spdlog::get("A")->set_level(spdlog::level::trace);
    }
    int skipped_scan_num = 0;
    for (int i = 0; i < df->getSpectrumNum(); i++) {
        CSpectrum *spec = df->getSpectrum(i);
        if (spec->getMSLevel() != 2) {
            continue;
        }
        // MS2
        int scannum = spec->getScanNum();
        PSMInfo psmInfo;
        bool isfound = ppp.getPSMInfobyScanFileName(df->getSourceFileName(), scannum, psmInfo);
        if (not isfound) {
            skipped_scan_num++;
            //cout << scannum << " not found" << endl;
            continue;
        }

        string peptide = psmInfo.searchhits[0]->m_modified_peptide;
        if (CDebugMode::callDebug()->getMdebug()) {
            if (peptide != debug_peptide) {
                continue;
            }
            if (debug_scan != -1 and scannum != debug_scan) {
                continue;
            }
        }
        index.push_back(i);
    }
    cout << skipped_scan_num << " MS2 spectra skipped by search engine " << endl;
    cout << "Remaining MS2: " << index.size() << endl;
}

void ValidatePSM::run() {
    bool alwaysReCreateFeature = true;
    if (not File::isExist(m_output_feature_file) or alwaysReCreateFeature) {
        vector<Feature *> features = createFeatureVector(m_fragmodelfile,
                                                         m_ghost, m_minIntFC, m_fragmodelscoretype,
                                                         m_featurelistfile, m_binaryPath);
        Peptide::defaultTables();
        vector<vector<double>> featuretable;
        CTable psmtable;
        vector<string> header = {"filename", "id", "scan", "modpep", "charge", "protein", "ppprob", "iprob"};
        psmtable.setHeader(header);
        PeptideProphetParser ppp(m_pepxmlfile);
        for (const auto &eachfile : ppp.m_allSourceFiles) {
            DataFile *df = new DataFile(eachfile, 0, -1);

            vector<int> index;
            getProcessed_i(df, ppp, index);
            vector<vector<double>> onefeaturetable(index.size(), vector<double>(features.size(), 0.0));
            vector<vector<string>> onepsmtable;

            int MinSize = index.size() / m_threadNum + 1;
            vector<int> batches;
            for (int i = 0; i < index.size(); i += MinSize) {
                batches.push_back(i);
            }
            batches.push_back(index.size());
            vector<std::thread *> tasks;
            shared_ptr<Progress> ps(new Progress(index.size(), "Exporting feature"));

            for (int i = 1; i < batches.size(); i++) {  // there might be multiple source! in ppp
                spdlog::get("A")->trace("Start batch {}", i);
                std::thread *t = new std::thread(
                        std::bind(workOnOneBatch, &ppp, eachfile, df, &index, &batches, i, &features, &onefeaturetable,
                                  ps, fixMz, m_highMassAcc));
                if (!t) {
                    cout << "Wrong!----------------" << endl;
                }
                tasks.push_back(t);
            }
            for (auto &task : tasks) {
                if (task)
                    task->join();
            }
            if (CDebugMode::callDebug()->getMdebug()) {
                for (int i = 0; i < onefeaturetable.size(); i++) {
                    cout << "Starting row " << i + 1 << endl;
                    for (int j = 0; j < features.size(); j++) {
                        cout << features[j]->m_feature_name << ": " << onefeaturetable[i][j] << endl;
                    }
                    cout << "Ending row " << i + 1 << endl;
                }
            }
            updatePsmTable(ppp, df, index, psmtable);

            featuretable.insert(featuretable.end(), onefeaturetable.begin(), onefeaturetable.end());
        }

        exportTestingFeature(features, m_output_feature_file, featuretable);

        string corresponding_psm = m_output_feature_file.replace(m_output_feature_file.length() - strlen(".feat"),
                                                                 strlen(".feat"), ".psm");

        psmtable.saveAs(corresponding_psm, true, ',');

        for (auto x: features) {
            delete x;
        }
    }

    for (auto x: m_other_steps) {
        spdlog::get("A")->info("Running {}", x->getname());
        x->run();
    }
}

void
ValidatePSM::updatePsmTable(PeptideProphetParser &ppp, DataFile *df, const vector<int> &index, CTable &psmtable) const {
    int decoy2targetNum = 0, decoynum = 0, targetnum = 0;
    unordered_set<string> decoy2target_peptides;
    bool verbose = false;
    for (int idx : index) {
        CSpectrum *spec = df->getSpectrum(idx);
        PSMInfo psmInfo;
        bool isfound = ppp.getPSMInfobyScanFileName(df->getSourceFileName(), spec->getScanNum(), psmInfo);
        string proteinACNum = psmInfo.searchhits.at(0)->m_protein;
        SearchHit &sh0 = *psmInfo.searchhits[0];

        if (m_useAlternativeProtein and sh0.m_protein != psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(
                m_useAlternativeProtein, 0)) {
            proteinACNum = psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(m_useAlternativeProtein, 0);
            if (verbose) {

                cout << "Protein changed from: " << sh0.m_protein
                     << " to " << proteinACNum
                     << " for peptide " << sh0.m_peptide << endl;
            }
            decoy2targetNum++;
            decoy2target_peptides.insert(sh0.m_peptide);
        }
        if (psmInfo.isDecoy(m_useAlternativeProtein, 0)) {
            decoynum++;
        } else {
            targetnum++;
        }
//        psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(), // use proteinnACNum now
        vector<string> tmp = {
                df->getSourceFileName(),
                to_string(idx),
                to_string(spec->m_scanNum),
                sh0.m_modified_peptide,
                to_string(psmInfo.charge),
                proteinACNum,
                to_string(sh0.m_peptideprophet_prob),
                to_string(sh0.m_iprophet_prob)
        };
        psmtable.addRow(tmp);
    }
    cout << "[Info] " << decoy2targetNum << " PSMs (" << decoy2target_peptides.size()
         << " Peptides) changed from DECOY to TARGET; Decoy: " << decoynum
         << " Target: " << targetnum << " Total:  " << index.size()
         << endl;
}

void ValidatePSM::exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                                       const vector<vector<double>> &featuretable) const {
    ofstream fout(feature_outfile.c_str(), ios_base::out);
    for (auto x: features) {
        fout << x->m_feature_name << ",";
    }
    fout << "class" << endl;

    for (const auto &i : featuretable) {
        for (double j : i) {
            fout << j << ",";
        }
        fout << "-1" << endl;
    }
    fout.close();
    cout << "Feature exported to file " << feature_outfile << endl;
}

void
workOnOneBatch(PeptideProphetParser *ppp, const string &mzMLsourcefile, DataFile *df, vector<int> *index,
               vector<int> *batches,
               int i, vector<Feature *> *features, vector<vector<double>> *featuretable, shared_ptr<Progress> ps,
               bool fixMz, bool highMassAcc) {
    spdlog::get("A")->trace("Running batch {} from {} to {}", i, batches->at(i - 1), batches->at(i));
    string name = File::CFile(df->getSourceFileName()).basename;
    for (int j = batches->at(i - 1); j < batches->at(i); j++) {
        ps->increase();
        int idx = index->at(j);
        CSpectrum *spec = df->getSpectrum(idx);
        if (spec == nullptr) {
            spdlog::get("A")->info("Error: spec is null idx = {}", idx);
            return;
        }
        PSMInfo psmInfo(name, spec->getScanNum());
        bool isfound = ppp->getPsmInfo(psmInfo);
        get_feature_from_spec(&psmInfo, spec, &(*featuretable)[j], features, idx, mzMLsourcefile, fixMz, highMassAcc);
    }
}

void plot_ROC_with_score(CTable &psmtable, int column, const string &scoreName, string outprefix, bool onlytarget,
                         HTMLReporter *reporter) {
    cout << "Plot ROC curve with ground truth, with TARGET_ONLY = " << onlytarget << endl;
    if (onlytarget) {
        outprefix += "_onlyT_";
    }
    vector<double> pos, neg;
    for (int i = 0; i < psmtable.m_row; i++) {
        if (psmtable.getEntry(i, psmtable.getColByHeader("isdecoy")) == "D" and onlytarget) {
            continue;
        }

        double score = atof(psmtable.getEntry(i, column).c_str());
        if (psmtable.getEntry(i, psmtable.getColByHeader("istruth")) == "1") {
            pos.push_back(score);
        } else {
            neg.push_back(score);
        }
    }
    plotScoreDistribution(pos, neg, outprefix + "_scoreDist.png", scoreName);

    CROCPlot rocplot(pos, neg);
    rocplot.saveROCtoTXT(outprefix + "_roc.txt");
    double auc = rocplot.getAUC();
    string outfigure = outprefix + "_roc_plot.png";
    rocplot.plotROCtoPNG(auc, outfigure, scoreName);
    if (reporter != nullptr) {
        reporter->addImage(outfigure, outfigure, "");
        string second_out_fig = outfigure.substr(0, outfigure.length() - 4) + "_full.png";
        reporter->addImage(second_out_fig, second_out_fig, "");
    }
    spdlog::get("A")->info("AUC of score {} is {:.4f}; ", scoreName, auc);
}

shared_ptr<CTruth> CreateTruth(const string &filename, const string &method) {
    cout << "Truthfile " << filename << endl << "method: " << method << endl;
    if (not File::isExist(filename, true)) {
        cout << "File does not exist!!!" << filename << endl;
        return nullptr;
    } else if (method == "scan+pep") {
        return shared_ptr<CTruth>(new CScanPepList(filename));
    } else if (method == "pepset") {
        return shared_ptr<CTruth>(new CTruthPepList(filename));
    } else {
        spdlog::get("A")->info("Invalid Truth Method: {}", method);
        throw CException(R"(Invalid truth method! Try "scan+pep" or "pepset"!)");
    }
}


void CROCPlot::plotROCtoPNG(double auc, const string &roc_png_filename, const string &scorename) {
    vector<tuple<double, double, double>> roc_tuple;
    for (auto &x: m_struct_roc) {
        roc_tuple.emplace_back(x.fpr, x.tpr, x.score);
    }

    string strAUC = to_string_with_precision(auc, 4);
    vector<string> figurefilelist;
    {
        Gnuplot gnuplot;
        gnuplot << "set terminal pngcairo enhanced size 600,600" << endl;
        gnuplot << "set size square" << endl;
        gnuplot << "set title 'ROC curve (partial)'" << endl;
        gnuplot << "set output '" << roc_png_filename << "'" << endl;
        gnuplot << "set xlabel 'false positive rate (FPR)'" << endl;
        gnuplot << "set ylabel 'true positive rate (TPR)'" << endl;
        gnuplot << "set xrange [0:0.02]" << endl;
        gnuplot << "set yrange [0.5:1]" << endl;
        gnuplot << "set key bottom right" << endl;
        gnuplot << "plot '-' using 1:2 with line lw 2 title \"" << scorename << " score AUC=" << strAUC << "\"" << endl;

        gnuplot.send1d(roc_tuple);
        figurefilelist.push_back(roc_png_filename);
    }
    {
        string fullaxis_png = roc_png_filename.substr(0, roc_png_filename.length() - 4) + "_full.png";
        Gnuplot gnuplot;
        gnuplot << "set terminal pngcairo enhanced size 600,600" << endl;
        gnuplot << "set size square" << endl;
        gnuplot << "set title 'ROC curve (full)'" << endl;
        gnuplot << "set output '" << fullaxis_png << "'" << endl;
        gnuplot << "set xlabel 'false positive rate (FPR)'" << endl;
        gnuplot << "set ylabel 'true positive rate (TPR)'" << endl;
        gnuplot << "set key bottom right" << endl;
        gnuplot << "plot '-' using 1:2 with line lw 2 title \" " << scorename << " score AUC=" << strAUC << "\""
                << endl;

        gnuplot.send1d(roc_tuple);
        figurefilelist.push_back(fullaxis_png);
    }
}

CROCPlot::CROCPlot(vector<double> &postiveScores, vector<double> &negativeScores) {
    int positive_num = postiveScores.size(), negative_num = negativeScores.size();
    spdlog::get("A")->info("positive: {}; negative: {}", positive_num, negative_num);

    sort(postiveScores.begin(), postiveScores.end(), [](const double &x, const double &y) { return x > y; });
    sort(negativeScores.begin(), negativeScores.end(), [](const double &x, const double &y) { return x > y; });

    int pcount = 0, ncount = 0;
    double score = 1;

    m_struct_roc.emplace_back(SROC(0, 0, 1, 0, 0));
    while (pcount < positive_num and ncount < negative_num) {
        if (postiveScores[pcount] > negativeScores[ncount]) {
            score = postiveScores[pcount];
            while (pcount < positive_num and postiveScores[pcount] >= score) pcount++;
            double TPR = (double) pcount / positive_num, FPR = (double) ncount / negative_num;
            m_struct_roc.emplace_back(SROC(FPR, TPR, score, ncount, pcount));
        } else if (postiveScores[pcount] < negativeScores[ncount]) {
            score = negativeScores[ncount];
            while (ncount < negative_num and negativeScores[ncount] >= score) ncount++;
            double TPR = (double) pcount / positive_num, FPR = (double) ncount / negative_num;
            m_struct_roc.emplace_back(SROC(FPR, TPR, score, ncount, pcount));
        } else {

            score = negativeScores[ncount];
            while (pcount < positive_num and postiveScores[pcount] >= score) pcount++;
            while (ncount < negative_num and negativeScores[ncount] >= score) ncount++;
            double TPR = (double) pcount / positive_num, FPR = (double) ncount / negative_num;
            m_struct_roc.emplace_back(SROC(FPR, TPR, score, ncount, pcount));
        }
    }
    while (pcount < positive_num) {
        score = postiveScores[pcount];
        pcount++;
        double TPR = (double) pcount / positive_num, FPR = 1.0;
        m_struct_roc.emplace_back(SROC(FPR, TPR, score, ncount, pcount));
    }

    while (ncount < negative_num) {
        score = negativeScores[ncount];
        ncount++;
        double TPR = 1.0, FPR = (double) ncount / negative_num;
        m_struct_roc.emplace_back(SROC(FPR, TPR, score, ncount, pcount));
    }
}

void CROCPlot::saveROCtoTXT(const string &roc_data) {
    ofstream fout(roc_data.c_str(), ios::out);

    for (auto x: m_struct_roc) {
        fout << x << endl;
    }
    fout.close();
}

double CROCPlot::getAUC() {
    // equation
    // \SUM { dx * (y_left+y_right) /2}
    /*
     *                   _________
     *              / |       |
     *  ________  /   |    y_right
     *  y_left  |     |       |
     *  _______ |_____|  _____|______
     *          | dx  |
     *
     * */
    //

    double auc = 0;
    for (int i = 1; i < m_struct_roc.size(); i++) {
        double dx = m_struct_roc.at(i).fpr - m_struct_roc.at(i - 1).fpr;
        double y_left = m_struct_roc.at(i - 1).tpr, y_right = m_struct_roc.at(i).tpr;
        auc += (dx) * (y_left + y_right) / 2.0;
    }

    return auc;
}

CTruth::CTruth() {

}

CTruth::~CTruth() {}

bool CTruth::validate(int scan, string modified_pep) // interface
{
    cout << "This function should not be called" << endl;
    return true;
}

string CTruth::getTruth(int scan) {
    return "N/A";
}

CScanPepList::CScanPepList(const string &scanpepListFile) {
    m_truthfile = scanpepListFile;
    CTable x(scanpepListFile, '\t', false);
    for (int i = 0; i < x.m_row; i++) {
        int scan = atoi(x.getEntry(i, 0).c_str());
        string peptide = x.getEntry(i, 1);
        string plainPepStr = regex_replace(peptide.c_str(), regex("[[0-9.]*]"), "");
        m_scanPepList.emplace_back(scan, plainPepStr);
    }
}

bool CScanPepList::validate(int scan, string modified_pep) {
    string plainPepStr = regex_replace(modified_pep.c_str(), regex("[[0-9.]*]"), "");
    spdlog::get("A")->debug("Validating scan + peptide: {} {} <-- {} ", scan, plainPepStr, modified_pep);
    bool ret = false;
    for (auto &i : m_scanPepList) {
        if (get<0>(i) == scan and get<1>(i) == plainPepStr) {
            ret = true;
        }
    }
    return ret;

}

string CScanPepList::getTruth(int scan) {
    string answer = "N/A";
    for (auto &i : m_scanPepList) {
        if (get<0>(i) == scan) {
            answer = get<1>(i);
        }
    }
    return answer;
}


CTruthPepList::CTruthPepList(const string &pepListFile) {
    m_truthfile = pepListFile;
    ifstream fin(pepListFile.c_str(), ios::in);
    string line;
    if (fin.is_open()) {
        while (getline(fin, line)) {
            // Error: sometimes there is a extra \r at the end of a line
            // The following code remove \r from the end
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            m_pepList.push_back(regex_replace(line, regex("L"), "I"));
        }
    }
    fin.close();
    spdlog::get("A")->info("{} peptides loaded from {}", m_pepList.size(), pepListFile);
//    for(int i = 0; i < m_pepList.size(); i ++){
//        cout << i << "[\"" <<m_pepList[i]<<"\"]"  << endl;
//    }
}

// 1. allow partial match
// consider the synthetic peptide data set, the peptide can be incorrect. e.g. one AA missing at n terminus?
// This can be true peptide in the sample, so we make them the correct peptide.
// 2. I == L
bool CTruthPepList::validate(int scan, string modified_pep) {
    bool ret = false;
    bool allow_partial_match = true;
    string plainPepStr = regex_replace(modified_pep.c_str(), regex("[[0-9.]*]"), "");
    plainPepStr = regex_replace(plainPepStr, regex("L"), "I"); // L to I

    if (m_pepList.end() != std::find(m_pepList.begin(), m_pepList.end(), plainPepStr)) {
        ret = true;// ret;
//        cout << modified_pep << "---> " << plainPepStr << " Status: " << ret << endl;
    }
    if (not ret and allow_partial_match) { // false, try partial match
        for (int i = 0; i < m_pepList.size(); i++) {
            int found = m_pepList[i].find(plainPepStr);
            if (found == string::npos) {
                // not found
                continue;
            } else {
                ret = true;
                cout << "Partial match " << plainPepStr << " in " << i << "-th peptide: " << m_pepList[i] << endl;
                break;
            }
        }
    }

    // remove the first nterm aa, try find again
    if (not ret and allow_partial_match) { // false, try partial match
        string shortplainPepStr = plainPepStr.substr(1);
        for (int i = 0; i < m_pepList.size(); i++) {
            int found = m_pepList[i].find(shortplainPepStr);
            if (found == string::npos) {
                // not found
                continue;
            } else {
                ret = true;
                cout << "Partial match " << shortplainPepStr << " from " << plainPepStr << " found in " << i
                     << "-th peptide: " << m_pepList[i] << endl;
                break;
            }
        }
    }

    return ret;
}

HTMLReporter::HTMLReporter(string filename, const string &title) {
    m_outfilename = filename;
    m_fout.open(m_outfilename.c_str(), ios::out);
    m_fout << "<!DOCTYPE html>\n"
              "<html lang=\"en\">\n"
              "<head>\n"
              "    <meta charset=\"UTF-8\">\n"
              "    <title>" << title << "</title>\n"
                                        "</head>\n"
                                        "<body>\n"
           << endl;
}

HTMLReporter::~HTMLReporter() {
    m_fout << "</body>\n"
              "</html>" << endl;
    m_fout.close();
    spdlog::get("A")->info("HTML report file {} generated!", m_outfilename);
}

void HTMLReporter::addImage(const string &filename, const string &title, const string &description) {
    string outstr = "<figure>\n";
    outstr += "  <img src=\"" + filename + "\" >\n";
    outstr += "  <figcaption>" + title + " </figcaption>\n";
    outstr += "</figure>\n";
    m_fout << outstr << endl;
}

resultAnalysis::resultAnalysis(string pepxmlfile, shared_ptr<CTruth> truth, string basename, bool useAlternativeProt) {
    m_pepxmlfile = pepxmlfile;
    m_name = "Result analysis";
    m_truth = truth;
    m_basename = basename;
    m_useAlternativeProt = useAlternativeProt;
}

struct SFdrScore {
    double m_fdr;
    int m_corr;
    int m_incorr;
    double m_score;

    SFdrScore(double fdr, int corr, int incorr, double score) {
        m_fdr = fdr;
        m_corr = corr;
        m_incorr = incorr;
        m_score = score;
    }
};

class CFdrNumCorr {
   //vector<tuple<double, double>> m_fdr_counts;
    vector<SFdrScore> m_FdrScore;

public:
    CFdrNumCorr() {
        ;
    }

    void plot(const string &title, const string &filename) {
        SFdrScore afs(0.01, 0, 0, 1.0);
        auto y = upper_bound(m_FdrScore.begin(), m_FdrScore.end(), afs,
                             [](const SFdrScore &u, const SFdrScore &v) -> bool { return u.m_fdr < v.m_fdr; });
        spdlog::get("A")->info("{}: FDR: {:.04f}, #PSM: {}, #PSM_Incorr: {}, Probability cut-off: {:.04f}", title,
                               y->m_fdr, y->m_corr, y->m_incorr, y->m_score);
        vector<tuple<double, double>> tmpobj;
        for (auto each: m_FdrScore) {
            tmpobj.emplace_back(each.m_fdr, each.m_corr);
        }
        plot_FDR_curve(filename, title, tmpobj);
    }

    ~CFdrNumCorr() = default;

    void saveAs(const string &outputfile) {
        ofstream fout(outputfile.c_str(), ios::out);
        for (auto x: m_FdrScore) {
            fout << x.m_fdr  << "\t" << x.m_corr << "\t" << x.m_incorr << "\t" << x.m_score<< endl;
        }
        fout.close();
        cout << "Data exported to file : " << outputfile << endl;
    }

    void export_peptideprophetinfo(PeptideProphetParser &ppp, bool useAlternativeProt = true) {
        vector<double> tProbs, dProbs;
        extract_pep_prob(tProbs, dProbs, ppp, get_peptidePropeht_prob, useAlternativeProt);

        get_FDR_CorrectNum(tProbs, dProbs);
    }

    void get_FDR_CorrectNum(vector<double> &tProbs, vector<double> &dProbs) {
        // sort to be descendent order
        sort(tProbs.begin(), tProbs.end(), [](const double &x, const double &y) -> bool { return x > y; });
        sort(dProbs.begin(), dProbs.end(), [](const double &x, const double &y) -> bool { return x > y; });

        int tcount = 0, dcount = 0;
        // fixed initialization value
        m_FdrScore.emplace_back(SFdrScore(0.0, 0, 0, 1.0));
        //m_fdr_counts.emplace_back(0.0, 0.0);
        while (tcount < tProbs.size() and dcount < dProbs.size()) {
            double score = 1.0;
            if (tProbs[tcount] > dProbs[dcount]) {
                score = tProbs[tcount];
                tcount++;
            } else if (tProbs[tcount] < dProbs[dcount]) {
                score = dProbs[dcount];
                dcount++;
            } else{
                score = tProbs[tcount];
                dcount ++;
                tcount ++;
            }

            //  FDR = D/T x 100%, shall we try pi_0 ?
            //  NO,we are using concatenated search here!
            double FDR = 1.0;
            if(tcount>0){
                FDR = dcount * 1.0 / tcount;
            }
            //m_fdr_counts.emplace_back(FDR, tcount * 1.0);
            if(score < m_FdrScore.back().m_score or fabs(FDR - m_FdrScore.back().m_fdr) > 0.001) {
                // score decreasing or FDR diffs with 0.1% (FDR can increase or decrease)
                m_FdrScore.emplace_back(SFdrScore(FDR, tcount, dcount, score));

            }
        }
        // from FDR to q-value

        for(auto it= m_FdrScore.rbegin(); it+1!=m_FdrScore.rend(); it ++ ){
            // from end to begin
            if(it->m_fdr < (it+1)->m_fdr){
                // if current fdr is smaller than the previous one. reset to current one.
                (it+1)->m_fdr= it->m_fdr;
            }
        }

    }

    void export_iprophetinfo(PeptideProphetParser &ppp, bool useAlternativeProt = true) {
        vector<double> tProbs, dProbs;
        extract_pep_prob(tProbs, dProbs, ppp, get_iprophet_prob, useAlternativeProt);
        get_FDR_CorrectNum(tProbs, dProbs);
    }

};

// this analysis should be called: peptide prophet analysis
// it will not work with other format.
void resultAnalysis::run() {
    cout << "\n--------------\nResult analysis\n"
            "Input: interact-.ipro.pep.xml file (Sequence Database Search(TDC)+xinteract(-dDECOY -OAPd)), RF model prediction score;\n"
            "Search engine can be Comet/MSFragger;\n"
            "Output: FDR\n"
            "\n"
            "Step 1: import rf scores\n"
            "Step 2: collect decoy rf scores and target rf scores\n"
            "Step 3: plot FDR -- CorrectNum figure; and AUC" << endl;
    string nameprefix = m_pepxmlfile + m_basename;
    HTMLReporter hr(nameprefix + ".html");

    string m_psmfile = nameprefix + ".psm";
    string m_predictionfile = nameprefix + ".prediction";

    // output to table of everything
    CTable psmtable(m_psmfile, ',', true, 0),
            predictiontable(m_predictionfile, ' ', false, 3);

    predictiontable.setHeader({"rf+", "rf-"});
    psmtable.Join(predictiontable);

    vector<double> tProb, dProb;
    psmtable.appendHeader({"isdecoy", "istruth", "truthpep"});

    // Todo: improved later.
    map<string, int> col2index = {{"protein", 5},
                                  {"rf+",     8},
                                  {"iProb",   7},
                                  {"pProb",   6},
                                  {"peptide", 3},
                                  {"scan",    2}};

    for (int i = 0; i < psmtable.m_row; i++) {
        string protein = psmtable.getEntry(i, col2index["protein"]);// to check todo:X
        double rfscore = atof(psmtable.getEntry(i, col2index["rf+"]).c_str());// key

        if (string::npos == protein.find("DECOY")) {
            spdlog::get("A")->debug("not a decoy protein {}", protein);
            tProb.push_back(rfscore);
            psmtable.appendEntry(i, "T");
        } else {
            dProb.push_back(rfscore);
            psmtable.appendEntry(i, "D");
        }
    }

    if (m_truth != nullptr) {
        spdlog::get("A")->info("using ground truth: {}", m_truth->getTruthFileName());
        for (int i = 0; i < psmtable.m_row; i++) {
            int scan = atoi(psmtable.getEntry(i, col2index["scan"]).c_str());

            if (m_truth->validate(scan, psmtable.getEntry(i, col2index["peptide"]))) {
                psmtable.appendEntry(i, "1");
            } else {
                psmtable.appendEntry(i, "0");
            }
            psmtable.appendEntry(i, m_truth->getTruth(scan));
        }
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("rf+"), "Random Forest", nameprefix + "_rf", false, &hr);
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("iprob"), "iProphet", nameprefix + "_ip", false, &hr);
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("ppprob"), "PeptideProphet", nameprefix + "_pp", false,
                            &hr);
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("rf+"), "Random Forest", nameprefix + "_rf", true, &hr);
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("iprob"), "iProphet", nameprefix + "_ip", true, &hr);
        plot_ROC_with_score(psmtable, psmtable.getColByHeader("ppprob"), "PeptideProphet", nameprefix + "_pp", true,
                            &hr);

    } else {
        spdlog::get("A")->info("no ground truth");
    }


    CFdrNumCorr fdrRF, fdrPep, fdrIpro;
    fdrRF.get_FDR_CorrectNum(tProb, dProb);
    fdrRF.saveAs(nameprefix + "_FDR_CorrectNum_myProb.txt");
    PeptideProphetParser ppp(m_pepxmlfile);
    fdrPep.export_peptideprophetinfo(ppp, m_useAlternativeProt);
    fdrPep.saveAs(nameprefix + "_FDR_CorrectNum_pepprob.txt");
    fdrIpro.export_iprophetinfo(ppp, m_useAlternativeProt);
    fdrIpro.saveAs(nameprefix + "_FDR_CorrectNum_iprob.txt");


    map<string, string> filenames;
    filenames["RF"] = nameprefix + "_FDR_CorrectNum_myProb.png";
    filenames["iProb"] = nameprefix + "_FDR_CorrectNum_iprob.png";
    filenames["pProb"] = nameprefix + "_FDR_CorrectNum_pepprob.png";

    fdrRF.plot("RF score", filenames["RF"]);
    fdrIpro.plot("iProb", filenames["iProb"]);
    fdrPep.plot("pProb", filenames["pProb"]);

    hr.addImage(filenames["RF"], "FDR-Corrrect ID Plot for Random Forest score", "");
    hr.addImage(filenames["iProb"], "FDR-Corrrect ID Plot for iProphet score", "");
    hr.addImage(filenames["pProb"], "FDR-Corrrect ID Plot for PeptideProphet score", "");

    string pepxmlfeat = nameprefix + ".feat";
    CTable featuretable(pepxmlfeat, ',', true, 0);
    psmtable.Join(featuretable);
    // true FDR vs estimated FDR
    psmtable.saveAs(m_psmfile + ".all", true, ',');
}

//----------extract features ---
ExtractFeatures::ExtractFeatures(string pepxmlfile, const string &validatorModel, string fragmodelscoretype,
                                 double minInt,
                                 bool ghost, int threadNum, string featurelistfile, bool useAlternativeProt,
                                 string binaryPath) : m_delimiter_psm_file(',') {

    m_binaryPath = binaryPath;
    m_useAlternativeProtein = useAlternativeProt;
    m_name = "Feature extraction";
    fixMz = false;
    m_highMassAcc = true;
    m_featurelistfile = featurelistfile;
    m_validatorModelFile = validatorModel;
    m_pepxmlfile = pepxmlfile;
    m_fragmodelfile = getfragmodelfile(validatorModel, minInt);;
    m_fragmodelscoretype = fragmodelscoretype;
    m_minIntFC = minInt;
    m_ghost = ghost;
    m_basename = File::CFile(validatorModel).basename;

    m_output_feature_file = m_pepxmlfile + m_basename + ".feat";

    m_threadNum = threadNum <= 0 or threadNum > getProperThreads() ? getProperThreads() : threadNum;
}

void ExtractFeatures::getIndexPairs(DataFile *df, PeptideProphetParser &ppp, vector<pair<int, int>> &indexPair) {

//    for(int i = 0; i < ppp.getPSMNum(); i ++){
//        PSMInfo psminfo;
//        ppp.getPSMInfobyindex(i,psminfo);
//        psminfo.start_scan
//    }
}

void ExtractFeatures::getProcessed_i(DataFile *df, PeptideProphetParser &ppp, vector<int> &index) {
    string debug_peptide;
    int debug_scan = -1;
    if (CDebugMode::callDebug()->getMdebug()) {
        cout << "Please type in debug peptide : " << endl;
        cin >> debug_peptide;
        cout << "Debug peptide is : " << debug_peptide << endl;
        cout << "Please type in debug scan: " << endl;
        cin >> debug_scan;
        cout << "Debug scan num is : " << debug_scan << endl;
        spdlog::get("A")->set_level(spdlog::level::trace);
    }
    int skipped_scan_num = 0;
    for (int i = 0; i < df->getSpectrumNum(); i++) {
        CSpectrum *spec = df->getSpectrum(i);
        if (spec->getMSLevel() != 2) {
            continue;
        }
        // MS2
        int scannum = spec->getScanNum();
        PSMInfo psmInfo;
        bool isfound = ppp.getPSMInfobyScanFileName(df->getSourceFileName(), scannum, psmInfo);
        if (not isfound) {
            skipped_scan_num++;
            //cout << scannum << " not found" << endl;
            continue;
        }

        string peptide = psmInfo.searchhits[0]->m_modified_peptide;
        if (CDebugMode::callDebug()->getMdebug()) {
            if (peptide != debug_peptide) {
                continue;
            }
            if (debug_scan != -1 and scannum != debug_scan) {
                continue;
            }
        }
        index.push_back(i);
    }
    cout << skipped_scan_num << " MS2 spectra skipped by search engine " << endl;
    cout << "Remaining MS2: " << index.size() << endl;
}

void ExtractFeatures::run() {
    string filebasename = m_output_feature_file;
    int found = m_output_feature_file.find_last_of('.');
    if (found != string::npos) {
        filebasename = m_output_feature_file.substr(0, found);
    }
    cout << "filebasename: " << filebasename << endl;
    string corresponding_psm = filebasename + ".psm";
    string corresponding_psmtsv = filebasename + ".psm.tsv";
    cout << "Export to two files; " << corresponding_psm << "\t" << corresponding_psmtsv << endl;

    bool alwaysReCreateFeature = true;
    if (not File::isExist(m_output_feature_file) or alwaysReCreateFeature) {
        // todo: you do not know whether it is Xinteract exported ... PepXML or not...
        vector<Feature *> features = createFeatureVector(m_fragmodelfile,
                                                         m_ghost, m_minIntFC, m_fragmodelscoretype,
                                                         m_featurelistfile, m_binaryPath);
        Peptide::defaultTables();
        vector<vector<double>> featuretable;
        CTable psmtable;
        vector<string> header = {"filename", "id", "scan", "modpep", "charge", "protein", "ppprob", "iprob"};
        psmtable.setHeader(header);
        PeptideProphetParser ppp(m_pepxmlfile);
        for (const auto &eachfile : ppp.m_allSourceFiles) {
            DataFile *df = new DataFile(eachfile, 0, -1);

            vector<int> index;
            getProcessed_i(df, ppp, index);
            vector<vector<double>> onefeaturetable(index.size(), vector<double>(features.size(), 0.0));
            vector<vector<string>> onepsmtable;

            int MinSize = index.size() / m_threadNum + 1;
            vector<int> batches;
            for (int i = 0; i < index.size(); i += MinSize) {
                batches.push_back(i);
            }
            batches.push_back(index.size());
            vector<std::thread *> tasks;
            shared_ptr<Progress> ps(new Progress(index.size(), "Exporting feature"));

            for (int i = 1; i < batches.size(); i++) {  // there might be multiple source! in ppp
                spdlog::get("A")->trace("Start batch {}", i);
                std::thread *t = new std::thread(
                        std::bind(workOnOneBatch, &ppp, eachfile, df, &index, &batches, i, &features, &onefeaturetable,
                                  ps, fixMz, m_highMassAcc));
                if (!t) {
                    cout << "Wrong!----------------" << endl;
                }
                tasks.push_back(t);
            }
            for (auto &task : tasks) {
                if (task)
                    task->join();
            }
            if (CDebugMode::callDebug()->getMdebug()) {
                printFeatureTable(features, onefeaturetable);
            }
            updatePsmTable(ppp, df, index, psmtable);

            featuretable.insert(featuretable.end(), onefeaturetable.begin(), onefeaturetable.end());
        }
        exportTestingFeature(features, m_output_feature_file, featuretable);
        psmtable.saveAs(corresponding_psm, true, ',');
        psmtable.saveAs(corresponding_psmtsv, true, '\t');

        for (auto x: features) {
            delete x;
        }
    }
}

void
ExtractFeatures::printFeatureTable(const vector<Feature *> &features,
                                   const vector<vector<double>> &onefeaturetable) const {
    for (int i = 0; i < onefeaturetable.size(); i++) {
        cout << "Starting row " << i + 1 << endl;
        for (int j = 0; j < features.size(); j++) {
            cout << features[j]->m_feature_name << ": " << onefeaturetable[i][j] << endl;
        }
        cout << "Ending row " << i + 1 << endl;
    }
}

void
ExtractFeatures::updatePsmTable(PeptideProphetParser &ppp, DataFile *df, const vector<int> &index,
                                CTable &psmtable) const {
    int mixtureSpectraNum = 0;
    int decoy2targetNum = 0, decoynum = 0, targetnum = 0;
    unordered_set<string> decoy2target_peptides;
    bool verbose = true;
    for (int idx : index) {
        CSpectrum *spec = df->getSpectrum(idx);
        PSMInfo psmInfo;
        vector<PSMInfo> allPsmInfo;
        bool isfound = ppp.getPSMInfobyScanFileName(df->getSourceFileName(), spec->getScanNum(), psmInfo);
        bool isfoundAll = ppp.getAllPSMsWithScanFileName(df->getSourceFileName(), spec->getScanNum(), allPsmInfo);
        if (isfoundAll and allPsmInfo.size() > 1) {
            cout << "Warning: more than one PSM found for this scan (Mixture spectra)" << endl;
            for (auto &eachpsm: allPsmInfo) {
                eachpsm.print();
            }
            mixtureSpectraNum++;
        }
        string proteinACNum = psmInfo.searchhits.at(0)->m_protein;
        SearchHit &sh0 = *psmInfo.searchhits[0];

        if (m_useAlternativeProtein and sh0.m_protein != psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(
                m_useAlternativeProtein, 0)) {
            proteinACNum = psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(m_useAlternativeProtein, 0);
            if (verbose) {
                cout << "Protein changed from: " << sh0.m_protein
                     << " to " << proteinACNum
                     << " for peptide " << sh0.m_peptide << endl;
            }
            decoy2targetNum++;
            decoy2target_peptides.insert(sh0.m_peptide);
        }
        if (psmInfo.isDecoy(m_useAlternativeProtein, 0)) {
            decoynum++;
        } else {
            targetnum++;
        }
//        psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(), // use proteinnACNum now
        vector<string> tmp = {
                df->getSourceFileName(),
                to_string(idx),
                to_string(spec->m_scanNum),
                sh0.m_modified_peptide,
                to_string(psmInfo.charge),
                proteinACNum,
                to_string(sh0.m_peptideprophet_prob),
                to_string(sh0.m_iprophet_prob)
        };
        psmtable.addRow(tmp);
    }
    cout << "[Info] " << decoy2targetNum << " PSMs (" << decoy2target_peptides.size()
         << " Peptides) changed from DECOY to TARGET; Decoy: " << decoynum
         << " Target: " << targetnum << " Total:  " << index.size()
         << endl;
    cout << "[Info] " << mixtureSpectraNum << " mixture spectra detected" << endl;
}

void ExtractFeatures::exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                                           const vector<vector<double>> &featuretable) const {
    ofstream fout(feature_outfile.c_str(), ios_base::out);
    for (auto x: features) {
        fout << x->m_feature_name << ",";
    }
    fout << "class" << endl;

    for (const auto &i : featuretable) {
        for (double j : i) {
            fout << j << ",";
        }
        fout << "-1" << endl;
    }
    fout.close();
    cout << "Feature exported to file " << feature_outfile << endl;
}


ExtractFeaturesFromPepXML::ExtractFeaturesFromPepXML(string pepxmlfile, const string &validatorModel,
                                                     string fragmodelscoretype, double minInt,
                                                     bool ghost, int threadNum, string featurelistfile,
                                                     bool useAlternativeProt,
                                                     string binaryPath, int hitrank) {
    m_hitrank = hitrank;
    m_binaryPath = binaryPath;
    m_useAlternativeProtein = useAlternativeProt;
    m_name = "extract features from pepXML file";
    fixMz = false;
    m_highMassAcc = true;
    m_featurelistfile = featurelistfile;
    m_validatorModelFile = validatorModel;
    m_pepxmlfile = pepxmlfile;
    m_fragmodelfile = getfragmodelfile(validatorModel, minInt);;
    m_fragmodelscoretype = fragmodelscoretype;
    m_minIntFC = minInt;
    m_ghost = ghost;
    m_basename = File::CFile(validatorModel).basename;

    m_output_feature_file = m_pepxmlfile + m_basename + ".feat";

    m_threadNum = threadNum <= 0 or threadNum > getProperThreads() ? getProperThreads() : threadNum;

}

class FeatureExtractionThread : public ICThreadTask {
    ICPepXMLParser *cpx;
    int i;
    DataFile *df;
    string eachfile;
    bool m_highMassAcc;
    bool fixMz;
    vector<Feature *> *features;
    vector<vector<double>> *onefeaturetable;
    int m_hitrank;
public:
    FeatureExtractionThread(ICPepXMLParser *cpx_, int psmId, DataFile *df_, string filename, bool highMassAcc,
                            bool fixMZ, vector<Feature *> *features_, vector<vector<double>> *OneSpecFeature_,
                            int hitrank) {
        cpx = cpx_;
        i = psmId;
        df = df_;
        eachfile = filename;
        m_highMassAcc = highMassAcc;
        fixMz = fixMZ;
        features = features_;
        onefeaturetable = OneSpecFeature_;
        m_hitrank = hitrank;
    }

    void run() override {
        PSMInfo psminfo;
        cpx->getPSMInfobyindex(i, psminfo);
        int scanNum = psminfo.start_scan;
        // search for spec with scanNum
        int idx = df->getIdxByScan(scanNum);
        CSpectrum *spec = df->getSpectrum(idx);

        get_feature_from_spec(&psminfo, spec, &((*onefeaturetable)[i]), features,
                              idx, eachfile, fixMz, m_highMassAcc, m_hitrank);

    }
};

void ExtractFeaturesFromPepXML::run() {
    int found = m_output_feature_file.find_last_of('.');
    string filebasename;
    if (found == string::npos) {
        // not find
        filebasename = m_output_feature_file;
    } else {
        filebasename = m_output_feature_file.substr(0, found);
    }
    cout << "filebasename: " << filebasename << endl;
    string corresponding_psm = filebasename + ".psm";
    string corresponding_psmtsv = filebasename + ".psm.tsv";
    cout << "Export to two files; " << corresponding_psm << "\t" << corresponding_psmtsv << endl;

    bool alwaysReCreateFeature = true;
    if (File::isExist(m_output_feature_file) and not alwaysReCreateFeature) return;


    vector<Feature *> features = createFeatureVector(m_fragmodelfile,
                                                     m_ghost, m_minIntFC, m_fragmodelscoretype,
                                                     m_featurelistfile, m_binaryPath);
    Peptide::defaultTables();
    CTable psmtable;
    vector<string> header = {"filename", "id", "scan", "modpep", "charge", "protein", "ppprob", "iprob"};
    psmtable.setHeader(header);

    ICPepXMLParser *cpx = createPepXML(m_pepxmlfile);

    vector<vector<double>> featuretable(cpx->getPsmNum(), vector<double>(features.size(), 0.0));

    vector<string> allSourceFiles = cpx->getAllSoruceFiles();
    for (const auto &eachfile : allSourceFiles) {
        DataFile *df = new DataFile(eachfile, 0, -1);
        auto range = cpx->getIndexRange(eachfile);
        int sampleNum = range.second - range.first;

        vector<int> psmIdx(range.second - range.first, 0);
        std::iota(psmIdx.begin(), psmIdx.end(), range.first);

        vector<ICThreadTask *> tasks;
        vector<FeatureExtractionThread> vtask;
        // create task, then run them.
        for (int i = range.first; i < range.second; i++) {
            vtask.emplace_back(cpx, i, df, eachfile, m_highMassAcc, fixMz, &features,
                               &featuretable, m_hitrank);
        }
        for (auto &i : vtask) {
            tasks.push_back(&i);
        }
        CTaskPool::issueTasks(tasks, false, true, m_threadNum, eachfile);
        updatePsmTable(cpx, df, psmIdx, psmtable);
        delete df;
    }
    exportTestingFeature(features, m_output_feature_file, featuretable);//
    psmtable.saveAs(corresponding_psm, true, ',');
    psmtable.saveAs(corresponding_psmtsv, true, '\t');

    for (auto x: features) {
        delete x;
    }
    cout << "release cpx" << endl;
    delete cpx;

}

void ExtractFeaturesFromPepXML::exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                                                     const vector<vector<double>> &featuretable) const {
    ofstream fout(feature_outfile.c_str(), ios_base::out);
    for (auto x: features) {
        fout << x->m_feature_name << ",";
    }
    fout << "class" << endl;

    for (const auto &row_i : featuretable) {
        for (double feature_j : row_i) {
            fout << feature_j << ",";
        }
        fout << "-1" << endl;
    }
    fout.close();
    cout << "Feature exported to file " << feature_outfile << endl;
}

string i_th(int i) {
    string ret = "th";
    switch (i) {
        case 1:
            ret = "st";
            break;
        case 2:
            ret = "nd";
            break;
        case 3:
            ret = "rd";
            break;
        default:
            ret = "th";
    }
    return ret;
}

void ExtractFeaturesFromPepXML::updatePsmTable(ICPepXMLParser *pepxml, DataFile *df, const vector<int> &psmIdx,
                                               CTable &psmtable) const {
    int mixtureSpectraNum = 0;
    int decoy2targetNum = 0, decoynum = 0, targetnum = 0;
    unordered_set<string> decoy2target_peptides;
    bool verbose = true;
    string name = File::CFile(df->getSourceFileName()).basename;
    //pair<int, int> range;
    int num_search_hit_not_found = 0;
    for (int j : psmIdx) {
        PSMInfo psmInfo;
        pepxml->getPSMInfobyindex(j, psmInfo);
        int scan = psmInfo.start_scan;
        int idx = df->getIdxByScan(scan);
        CSpectrum *spec = df->getSpectrum(idx);


        vector<string> tmp = {
                df->getSourceFileName(),
                to_string(idx),
                to_string(spec->m_scanNum),
                "?",
                to_string(psmInfo.charge),
                "?",
                "?",
                "?"
        };

        if (psmInfo.searchhits.size() <= m_hitrank) {
            num_search_hit_not_found++;

        } else {
            string proteinACNum = psmInfo.searchhits.at(m_hitrank)->m_protein;
            SearchHit &sh0 = *psmInfo.searchhits[m_hitrank];

            if (m_useAlternativeProtein and sh0.m_protein != psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(
                    m_useAlternativeProtein, m_hitrank)) {
                proteinACNum = psmInfo.getProtein_UseAlterProteinIfItsNotDecoy(m_useAlternativeProtein, m_hitrank);
                if (verbose) {
                    cout << "Protein changed from: " << sh0.m_protein
                         << " to " << proteinACNum
                         << " for peptide " << sh0.m_peptide << endl;
                }
                decoy2targetNum++;
                decoy2target_peptides.insert(sh0.m_peptide);
            }
            if (psmInfo.isDecoy(m_useAlternativeProtein, m_hitrank)) {
                decoynum++;
            } else {
                targetnum++;
            }

            tmp[3] = sh0.m_modified_peptide;
            tmp[5] = proteinACNum.substr(0, proteinACNum.find_first_of(' '));
            tmp[6] = to_string(sh0.m_peptideprophet_prob);
            tmp[7] = to_string(sh0.m_iprophet_prob);

        }
        psmtable.addRow(tmp);
    }
    cout << "[Info] " << decoy2targetNum << " PSMs (" << decoy2target_peptides.size()
         << " Peptides) changed from DECOY to TARGET; Decoy: " << decoynum
         << " Target: " << targetnum << " Total:  " << psmIdx.size()
         << endl;
    cout << "[Info] " << mixtureSpectraNum << " mixture spectra detected" << endl;
    cout << "[Info] " << num_search_hit_not_found << " MS2 spectra do not have " << m_hitrank << i_th(m_hitrank)
         << " hit (0-hit is top hit)" << endl;
}

RFModelConfig::RFModelConfig(const string &configFileName) {
    m_exportName = configFileName;
    string line;
    ifstream fin(configFileName, ios::in);
    map<string, string> key_values;
    while (std::getline(fin, line)) {
        cout << line << endl;
        vector<string> tokens;
        split_string(line, tokens, '=');
        trim_space_only(tokens[0]);
        trim_space_only(tokens[1]);
        key_values[tokens[0]] = tokens[1];
    }
    m_ranger_binary_path = key_values.at("ranger_binary_path");
    m_tsvfile = key_values.at("feature_tsv_file");
    m_mtry = stringTo<int>(key_values.at("mtry"));
    m_ntree = stringTo<int>(key_values.at("ntree"));
    m_isTraining = stringTo<bool>(key_values.at("is_training"));
    m_maxDepth = stringTo<int>(key_values.at("maxDepth"));
    m_threadnum = stringTo<int>(key_values.at("threadnum"));
    m_rfModel = key_values.at("rf_model_file");
    m_probPrediction = stringTo<bool>(key_values.at("is_prob_pred"));

}

void RFModelConfig::write() {
    ofstream fout(m_exportName, ios::out);
    fout << "ranger_binary_path = " << m_ranger_binary_path << endl;
    fout << "mtry = " << m_mtry << endl;
    fout << "ntree = " << m_ntree << endl;
    fout << "feature_tsv_file = " << m_tsvfile << endl;
    fout << "threadnum = " << m_threadnum << endl;
    fout << "is_training = " << m_isTraining << endl;
    fout << "maxDepth = " << m_maxDepth << endl;
    fout << "is_prob_pred = " << m_probPrediction << endl;
    fout << "rf_model_file = " << m_rfModel << endl;
}

RFModelConfig::RFModelConfig(string featuretsv, bool isTraining, string RFmodelfile, bool probPrediction, int mtry,
                             int ntree, int maxdepth, string rangerbinary, string exportName) {
    //m_name = "Run ranger (random-forest)";
    m_tsvfile = featuretsv;
    m_isTraining = isTraining;
    m_rfModel = RFmodelfile;
    m_mtry = mtry;
    m_ntree = ntree;
    m_maxDepth = maxdepth;
    m_probPrediction = probPrediction;

    m_ranger_binary_path = rangerbinary;
    m_threadnum = getProperThreads();
    m_exportName = exportName;
}

RFModelConfig::~RFModelConfig() = default;

void RFModelConfig::print() {
    cout << "m_tsvfile = " << m_tsvfile << endl;
    cout << "m_isTraining = " << m_isTraining << endl;
    cout << "m_rfModel = " << m_rfModel << endl;
    cout << "m_mtry = " << m_mtry << endl;
    cout << "m_ntree = " << m_ntree << endl;
    cout << "m_maxDepth = " << m_maxDepth << endl;
    cout << "m_probPrediction = " << m_probPrediction << endl;

    cout << "m_ranger_binary_path = " << m_ranger_binary_path << endl;
    cout << "m_threadnum = " << m_threadnum << endl;
    cout << "m_exportName = " << m_exportName << endl;

}

const string &RFModelConfig::getMTsvfile() const {
    return m_tsvfile;
}

bool RFModelConfig::isMIsTraining() const {
    return m_isTraining;
}

const string &RFModelConfig::getMRfModel() const {
    return m_rfModel;
}

bool RFModelConfig::isMProbPrediction() const {
    return m_probPrediction;
}

int RFModelConfig::getMMtry() const {
    return m_mtry;
}

int RFModelConfig::getMNtree() const {
    return m_ntree;
}

int RFModelConfig::getMMaxDepth() const {
    return m_maxDepth;
}

const string &RFModelConfig::getMRangerBinaryPath() const {
    return m_ranger_binary_path;
}

int RFModelConfig::getMThreadnum() const {
    return m_threadnum;
}

const string &RFModelConfig::getMExportName() const {
    return m_exportName;
}


CException::CException(const char *msg) : std::runtime_error(msg) {}
