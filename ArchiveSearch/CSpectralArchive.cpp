//
// Created by wulong on 8/6/18.
//

#include <set>
#include <string>
#include <fstream>
#include <spdlog/spdlog.h>
#include "../librarymsms/XMLFileParser.h"
#include "../librarymsms/CPSMAnnotation.h"
#include "../librarymsms/pathManager.h"
#include "../librarymsms/CArchiveSearchReply.h"
#include "../librarymsms/DatabaseManager.h"
#include "../librarymsms/CAnnotationDB.h"
#include "../librarymsms/CMzFileReader.h"
#include "../librarymsms/ProteomicsDataTypes.h"
#include "../librarymsms/ConcreteQueries.h"
#include "../librarymsms/CFragmentIndex.h"
#include "../librarymsms/Util.h"

#include "dependentlibs/cudaproj/dpcuda.h"
#include "dependentlibs/pvalueEstimation/EMAlgorithm.h"
#include "../gnuplot-iostream/gnuplot-iostream.h"
#include "dependentlibs/spectralIndex/CKMeans.h"
#include "CSpectralArchive.h"
#include "dependentlibs/spectralIndex/CMultiIndices.h"
#include "dependentlibs/pvalueEstimation/pvalueCalc.h"
#include "dependentlibs/pvalueEstimation/CTailEstimation.h"
#include "dependentlibs/randomSpecParser.h" // getPeakList
#include "dependentlibs/MzRTMap.h"
#include "CTimerSummary.h"
using namespace std;

bool add_new_data_file_to_mzXMLList(string mzXMLList, string newdatafile) {
    bool gotnewdata = true;
    CTable mzxmls(mzXMLList, '\t', false, 0);
    for (int i = 0; i < mzxmls.m_row; i++) {
        string currentfilename = mzxmls.getEntry(i, 0);
        if (currentfilename == newdatafile) {
            gotnewdata = false;
            cout << "file exist already! " << currentfilename << " file id: " << i << endl;
            break;
        }
    }
    if (gotnewdata) {
        mzxmls.addRow({newdatafile});
        mzxmls.saveAs(mzXMLList, false, '\t');// mzXMLList updated
        cout << "new file added: " << newdatafile << endl;
    }
    return gotnewdata;
}

class CGtOutput {
    struct PSMSummary {
        string name;
        int significant_num;
        const int CHARGE_LEN;
        const int PEP_LEN;
        vector<int> chargeNum;//[CHARGE_LEN];
        vector<int> peplenNum;//[PEP_LEN];
        int otherpeplen;
        int othercharge;


        PSMSummary(string name):CHARGE_LEN(7), PEP_LEN(100),name(name) {
            significant_num = 0;
            otherpeplen = 0;
            othercharge = 0;
            chargeNum.assign(CHARGE_LEN,0);
            peplenNum.assign(PEP_LEN,0);

        }
        void printSummary(){
            cout << name << " ==========Summary=========" << endl
            << "significant: " << significant_num << endl
            << "charge: other " << othercharge << endl;
            for (int i = 0; i < CHARGE_LEN; i++) cout << "charge_" << i + 1 << ": " << chargeNum[i] << endl;
            cout << "peplen: other " << otherpeplen << endl;
            for (int i = 0; i < PEP_LEN; i++) {
                if (peplenNum[i] > 0) {
                    cout << "peplen_" << i + 1 << ": " << peplenNum[i] << endl;
                }
            }
            cout << name << " =========Summary=======done==" << endl;
        }
        // tsv
        void printSummary_tsv(){
            cout << name << " ==========Summary=========" << endl
            << "significant: " << significant_num << endl
            << "charge: other " << othercharge << endl;
            cout << "Charge";
            for(int i = 0; i < CHARGE_LEN; i ++){
                cout << "\t" << i+1 ;
            }
            cout << endl;
            cout << "#Spec";
            for(int i = 0; i < CHARGE_LEN; i ++){
                cout << "\t" << chargeNum[i] ;
            }
            cout << endl;

//            for (int i = 0; i < CHARGE_LEN; i++) cout << "charge_" << i + 1 << ": " << chargeNum[i] << endl;
            cout << "peplen: other (peptide lenght  0 for unidentified spectra) " << otherpeplen << endl;

            int col = 10;
            for(int j = 0; j < PEP_LEN; j += col){
                cout << "PepLen";
                for(int i = j; i < PEP_LEN and i < j + col; i ++){
                    cout << "\t" << i +1 ;
                }
                cout << endl;
                cout << "#Spec";
                for(int i = j; i < PEP_LEN and i < j+col; i ++){
                    cout << "\t" << peplenNum[i] ;
                }
                cout << endl;
            }
//            for (int i = 0; i < PEP_LEN; i++) {
//                if (peplenNum[i] > 0) {
//                    cout << "peplen_" << i + 1 << ": " << peplenNum[i] << endl;
//                }
//            }
            cout << name << " =========Summary=======done==" << endl;
        }
        ~PSMSummary() {
//            printSummary();
            printSummary_tsv();

        }

        void set(int charge, int peplen) {
            if (charge >= 1 and charge <= CHARGE_LEN) chargeNum[charge - 1]++;
            else othercharge++;
            if (peplen >= 1 and peplen <= PEP_LEN) peplenNum[peplen - 1]++;
            else otherpeplen++;
        }

    };

    vector<SPsmAnnotation> m_gtInDB;
    string m_outfile;
    int m_size;
public:
    CGtOutput(string name, vector<long> &queryidx, CAnnotationDB *annodb, string outputfilename, ICMzFile *mzfile) {
        cout << "generating summary of current index  " << queryidx.size() << endl;
        m_size = queryidx.size();
        m_outfile = outputfilename;
        m_gtInDB.assign(m_size, SPsmAnnotation());
//        cout << "retrival of ground truth " << endl;
        for (int i = 0; i < m_size; i++) {
            annodb->retrieveGtinfo(queryidx[i], m_gtInDB[i]);
        }
//        cout << "retrieval of peak num" << endl;
        for (int i = 0; i < m_size; i++) {
            m_gtInDB[i].peaknum = mzfile->getPeakNum(queryidx[i]);
        }
//        cout << "start summary" << endl;

        PSMSummary summary("background peptides");
        ofstream fout(m_outfile.c_str(), ios::out);
        for (int i = 0; i < m_size; i++) {
            SPsmAnnotation &gt = m_gtInDB[i];
            summary.set(gt.charge, gt.peptideseq == "UNKNOWN" ? 0 : gt.peptideseq.size());
            summary.significant_num += gt.significance;
            fout << gt << endl;
        }
//        cout << "done summary" << endl;
    }

    ~CGtOutput() {

    }
};


// init by the name, then add the list of mzXMLs one by one.
CSpectralArchive::CSpectralArchive(string mzXMLList, string pepxml, string indexfile, bool removeprecursor,
                                   bool useflankingbins, int tol, int minPeakNum, bool myOwnIndex, CPQParam option,
                                   string indexstrings,
                                   bool usegpu, bool rebuildsqldb, int seedpvalue, const int topPeakNum,
                                   bool createfilenameBlackList, bool saveBackgroundScore, bool verbose,
                                   string archivename, string indexshuffleseeds) : PeakNumPerSpec(
        topPeakNum) {
            m_savebackgroundscore = saveBackgroundScore;
    // checking the parameters.
    if(mzXMLList.empty()){
        cerr << "The archive file name is not provided. use --mzxmlfiles or -m to specify the file contains a list of ms/ms data files. (or use the config file to provide the archive filename.)" << endl;
        throw runtime_error("empty archive name. ");
    }
    m_verbose = verbose;
    cout << "Start to create archive object" << endl;
    cout << "three steps" << endl 
    << "1. create sql database " << endl 
    << "2. create mz file " << endl 
    << "3. create index " << endl;
    m_AnnotationDB = make_shared<CAnnotationDB>(createfilenameBlackList);
    m_usegpu = usegpu;
    m_mzXMLListFileName = mzXMLList;
    if (archivename.empty()){
        // empty archive name, then use the mzXMListFileName.
        archivename = m_mzXMLListFileName;
    }
    m_tol = tol;
    m_removeprecursor = removeprecursor;
    m_useflankingbins = useflankingbins;
    m_pepxmlFileName = pepxml;
    
    // fixed index path issue. using the same path as archive settings.
    if(indexfile.empty()){
        m_indexFileName = File::CFile(m_mzXMLListFileName).path + "/";
    } else{
        m_indexFileName = indexfile;
    }
    
    m_minPeakNum = minPeakNum;
    m_dim = 4096;
    // -------------------------------------------parameters ready ----------------------------




    cout << "[info] 1/3 Connecting to SQL Database " << endl;
    m_AnnotationDB->createDatabase(rebuildsqldb, m_mzXMLListFileName + ".sqlite3db", m_verbose);

    cout << "[info] 2/3 Initialize parameters of peak list file " << endl;
    m_csa = CMzFileReader::makeShared(m_mzXMLListFileName, false, true, m_removeprecursor, intTol2double(m_tol),
                                       m_minPeakNum, verbose);

    cout << "[info] 3/3 Initialize index" << endl;
    // pass seeds in.
    createIndices(myOwnIndex, make_shared<CPQParam>(option), indexstrings, indexshuffleseeds);


    // loading index  create if not exist.
    // load index, create if not exist.
    if(m_indices->isExist()){
        m_indices->loadIndices();
    }    else{
        m_indices->trainOnFileList(m_mzXMLListFileName);
    }
    m_indices->display();
    // index ready. either trained or loaded.
    // database connection ready.
    // mz file to be loaded.



    // --------------------Now add some data files.
    this->update("","","",m_mzXMLListFileName);

    cout << "Index created" << endl;
    m_AnnotationDB->populateGtfilesTable(m_pepxmlFileName);  // GTFILES not used!
    if (m_usegpu)m_indices->toGpu(); // should we delay this? todo:
    // to start the service.
    getScorerFactoryPtr();
    int threadnum=1;
    m_pc = make_shared<CPValueCalculator>(m_pScorer->getSpecNum(), threadnum, seedpvalue, m_tol, saveBackgroundScore,verbose);
    m_pc->buildFragIndex(m_pScorer, verbose);

    // generate summary of archive.
    vector<long> queryidxlist = m_pc->getIndexList();
    CGtOutput gtoutput("background spectra summary", queryidxlist, m_AnnotationDB.get(), m_mzXMLListFileName + "_gtout.tsv", m_csa.get());

    updateIndex(m_verbose);
    m_AnnotationDB->createBlackListWithCE(false);
    m_AnnotationDB->fixAllChargeState();
}

// function not used!
void CSpectralArchive::createMzFileObject() {
    m_csa = CMzFileReader::makeShared(m_mzXMLListFileName, false, true, m_removeprecursor, intTol2double(m_tol),
                                       m_minPeakNum, false);
    getScorerFactoryPtr();
}

// init parameters for each index. filename, seed.etc..
// create path multiIndices under the same path as m_indexFileName.
// get number of indices to be created, using indexstr, split by ';'
// create multiIndices object. each corresponding to a file in multiIndices folder. 
void CSpectralArchive::createIndices(bool myOwnIndex, shared_ptr<CPQParam> option, string &indexstrings, string indexshuffleseeds) {
    int cnts = count(indexstrings.begin(), indexstrings.end(), ';') + 1;
    if (cnts > 6 and cnts < 1) {
        cout << "Error: invalid cnts " << cnts << endl;
        throw runtime_error("Invalid parameters!");
    }

    File::CFile fileObj(m_indexFileName);
    string path = fileObj.path + "/multiIndices";
    createDirectory(path);
    double index_tol = intTol2double(m_tol);
    // pass seeds from here.
    m_indices = make_shared<CMultiIndices>(indexstrings, indexshuffleseeds,path, fileObj.basename, true, index_tol, myOwnIndex, option,
                                           PeakNumPerSpec, m_removeprecursor, m_useflankingbins, m_dim);


}

CSpectralArchive::~CSpectralArchive() {}

void CSpectralArchive::update(string new_experimental_data, string new_search_result, string new_search_result_list,
                              string new_experimental_datalist) {
    updateIndex(m_verbose);
    bool newFilesAdded = false;
    addRawData(new_experimental_data, newFilesAdded);
    addListOfRawData(new_experimental_datalist, newFilesAdded);
    // after updated raw and list of raw, update the list of mzxmllist file. 
    if(newFilesAdded){
        vector<string> datafiles = m_AnnotationDB->getListOfSpecFiles();
        File::saveas(datafiles, m_mzXMLListFileName, true);
        cout << m_mzXMLListFileName << " updated, number of raw files " << datafiles.size() << endl;

        string scanFilename = m_mzXMLListFileName + ".scan";
        bool found = File::isExist(scanFilename, true);
        cout << "Removing scan file: " << scanFilename << endl;
        int rc = std::remove(scanFilename.c_str());
        if(found and rc) { 
            perror("remove of .scan file fails"); 
            throw runtime_error("\nERROR: Fail to remove "+ scanFilename +" file. Program will exit. \nPlease try to manually remove the file, and rerun. ");
        }
    }
    
    addSearchResult(new_search_result);
    addListOfSearchResults(new_search_result_list);
    m_indices->write();  // as there is only one file, we could write here.

    // mz file loaded.
    m_csa->refresh_mz_object_from_disk();

    cout << "[Info] size of archive: " <<this->size() << endl;
    
}

// update list of spectra
void CSpectralArchive::updateListOfSpectra(string spectraList) {
    // if there are too many spectra, try in batch mode.
    // checking existense of spectra
    // add to raw data; mz and index

    // add annotation. SQL-DB

}

void CSpectralArchive::addListOfRawData(const string &new_experimental_datalist, bool &newFileAdded) {
    if (new_experimental_datalist != "") {
        vector<string> files = readlines(new_experimental_datalist);
        int L = files.size();
        for(int i = 0; i < L; i ++){
            if(files[i].empty() or m_AnnotationDB->isSpecFileExist(files[i])) continue;
            files.push_back(files[i]);
        }
        files.erase(files.begin(), files.begin()+L);
        cout << files.size() << " out of " << L << " files will be added" << endl;
        if ( files.empty()){
            return; 
        }
        else{
            SimpleTimer st;
            for (int i = 0; i < files.size(); i++) {
                cout << "Processing " << i + 1 << " / " << files.size() << " file " << files[i] << endl;
                addRawData(files[i], newFileAdded);
                double time_elapsed = st.secondsElapsed();
                double speed = time_elapsed/(i+1);
                double time_to_be_used = (files.size()-1-i)*speed;

                spdlog::get("A")->info("time used {:.1f}s, {:.1f}s per file, ETA: {:.1f}s, -> {} hours {} minutes", time_elapsed, speed, time_to_be_used, int(time_to_be_used)/3600, (int(time_to_be_used)%3600)/60 );
            }
        }
       
    }
}

string CSpectralArchive::searchPeptide(string peptide) {
    cout << "Start searching for peptide " << peptide << endl;
    auto dbentry = m_AnnotationDB->searchPeptide(peptide);
//    shared_ptr<CDBEntry> dbentry = nullptr;
//    m_AnnotationDB->searchPeptide(peptide, dbentry);
    return m_AnnotationDB->toJsonNodes(*dbentry);
//
//    m_AnnotationDB->searchPeptide(peptide, results);
//    cout << "End of searching for peptide " << peptide << endl;
}

string
CSpectralArchive::searchFileName(string filename, string startscan, string endscan) {
    cout << "start searching for filename " << filename << endl;
    shared_ptr<CDBEntry> dbentry = nullptr;
    m_AnnotationDB->searchGTWithFileName(filename, startscan, endscan, dbentry);
    cout << "end of searching for filename" << endl;
    return m_AnnotationDB->toJsonNodes(*dbentry);
}

CMzSpec CSpectralArchive::getMzSpec(long queryindex) {
    return m_pScorer->getMzSpec(queryindex);
}

long CSpectralArchive::size() {
    long specNumDB = m_AnnotationDB->getTotalSpecNum();
    long specNum_index = m_indices->total();
    long specNum_mz = m_csa->getSpecNum();
    if(specNumDB != specNum_index or specNum_index !=specNum_mz){
        cout << "Inconsistent archive: \nSQL-DB\tIndex\tMzFile\n" << specNumDB << "\t" << specNum_index << "\t" << specNum_mz << endl;
        throw runtime_error("Inconsistent Archive size");
//        exit(-1);
    }
    return m_indices->total();
}

void CSpectralArchive::setnProbe(int nprobe) {
    m_indices->setNprobe(nprobe, false);
    agtsummary.setNProbe(nprobe);
    // if(nprobe != agtsummary.m_nprobe){
    //     agtsummary.m_nprobe = nprobe;
    //     agtsummary.correctTNNnum=0;
    //     agtsummary.totalTNNnum = 0;
    // }
    
}

// this is relies on the mzxml and pepxml relationship.
// todo: change the algorithm to allow interact ...
void CSpectralArchive::addSearchResult(string pepxmlfile) {
    if (pepxmlfile.empty()) {
        // cout << "No search result file provides! The following files will not be updated: \n"
        //         "pepxml table, and gt table. "
        //      << endl;
        return;
    }
    SimpleTimer st("update groundtruth table");
    doAddSearchResult(pepxmlfile);
}



void CSpectralArchive::addRawData(string mzXMLfile, bool &newFileAdded) {
    if (mzXMLfile.empty()) {
        return;
    }

    if (m_AnnotationDB->isSpecFileExist(mzXMLfile)) {
        cout << "File already exist in sqlite3 database, will not add again: " << mzXMLfile << endl;
    } else {
        newFileAdded = true;
        DataFile df(mzXMLfile, 0, -1);
        m_AnnotationDB->appendNewDataFile(df);
        cout << "DB updated" << endl;
        m_csa->append(df, m_removeprecursor, nullptr);
        cout << "DB MZ updated" << endl;
        m_indices->append(df);
        cout << "DB MZ INDEX updated" << endl;
        // m_indices->write();  // as there is only one file, we could write here.
        // the mzXML file will not be updated.
        // appendFileName(mzXMLfile);
    }
}

// this function should not be used.
void CSpectralArchive::appendFileName(const string &mzXMLfile) {
    // append the new mzXMLfile into the exist file list.
    vector<string> datafiles = readlines(m_mzXMLListFileName);

    auto it = find(datafiles.begin(), datafiles.end(), mzXMLfile);
    if (it == datafiles.end()) {
        cout << "data file not find in " << m_mzXMLListFileName << endl;
        datafiles.push_back(mzXMLfile);
        File::saveas(datafiles, m_mzXMLListFileName, true);
    }
}

void CSpectralArchive::updateIndex(bool verbosity) {
    m_AnnotationDB->createTable("SPECFILES", false,m_verbose);
    syncIndicesWithSpecFileTable();
    if (not m_AnnotationDB->isTableExist("GROUNDTRUTH")) {
        m_AnnotationDB->createTable("GROUNDTRUTH", false,m_verbose);
        populateGroundTruthTable();
    }

    m_AnnotationDB->addColumnCEIfNotExist(verbosity);
    m_AnnotationDB->addColumnAlterPepIfNotExist(verbosity);
    m_AnnotationDB->addColumnRFScoreIfNotExist(verbosity);
    m_AnnotationDB->addColumnNeighborIfNotExist(verbosity);
}

void CSpectralArchive::addListOfSearchResults(string pepxmlfilelist) {
    if (pepxmlfilelist.empty()) {
        // cout << "No search result list provided! The following files files will not be udpated:\n "
        //         "pepxml table, and gt table." << endl;
        return;
    }
    CTable pepxmllist(pepxmlfilelist, '\t', false, 0);
    cout << "reading pepxml file from " << pepxmlfilelist << endl;

    // Progress ps(pepxmllist.m_row, "Updating ground truth");
    for (int i = 0; i < pepxmllist.m_row; i++) {
        // ps.increase();
        string currentfile = pepxmllist.getEntry(i, 0);
        cout << "\nUpdating GROUNDTRUTH table: " << i + 1 << " / " << pepxmllist.m_row << "  " << endl;
        doAddSearchResult(currentfile);
    }
}


// todo: to be removed addpepxmlfile doing nothing
void CSpectralArchive::doAddSearchResult(string gtfile) {
    addpepxmlfile(m_pepxmlFileName, gtfile);
    shared_ptr<CGtUpdater> updator = m_AnnotationDB->buildGtUpdater(gtfile);
}

int CSpectralArchive::getDim() {
    return m_dim;
}

void CSpectralArchive::populateGroundTruthTable() {
    cout << "under construction" << endl;
}

void CSpectralArchive::syncIndicesWithSpecFileTable() {
    if (not m_AnnotationDB->isTableExist("SPECFILES")) {
        spdlog::get("A")->error("table SPECFILES does not exist!");
    }
    vector<string> results;
    long totalspecnum = m_AnnotationDB->getTotalSpecNum();
    if (m_indices->total() < totalspecnum) {
        spdlog::get("A")->error("index size {} is different from spectra number in spec file table {}",
                                m_indices->total(), totalspecnum);

        long index_total = m_indices->total();
        int startFileId = m_AnnotationDB->getFileIDWithEndIdx(index_total);

        int numFile = m_AnnotationDB->getTotalFileNum();
        for (int i = startFileId; i < numFile; i++) {
            string filename = m_AnnotationDB->getSpecFileNameByID(i);
            DataFile df(filename, 0, -1);
            spdlog::get("A")->info("updating {} / {}\tFile name: {}", i+1, numFile, filename);
            m_indices->append(df);
            if(i%5==0 or i == numFile-1) m_indices->write();
        }
//        if (m_indices) m_indices->write();
    }
}

// Note: query id can be -1
// When first query id is -1, we search with the spectrum in query pointer!
void CSpectralArchive::searchQuery(long query_index, string &jsonstring, int topN, int calcEdge, int nprobe,
                                   vector<uint16_t> &query, bool visualize, double minTNNDP, int indexNum, int TNNtopK) {
    setnProbe(nprobe);
    if(agtsummary.m_recallOfTrueNeighbor){
        if(minTNNDP>=0.0 and minTNNDP<=1.0){
            agtsummary.setRecallTNNMinDP(minTNNDP);
        }
        if(indexNum>=1 and indexNum<=m_indices->getNum()){
            agtsummary.setIndexNum(indexNum);
        }
        if(TNNtopK>=1){
            agtsummary.setRecallTNNTopK(TNNtopK);
        }
    }
    SimpleTimer st("query");
    CArxivSearchResult annResults;
    vector<long> queryIds = {query_index};
    vector<int> dpscores; // background score
    ICQuery *q = createICQuery(&query, &queryIds, m_useflankingbins, getDim(), *m_pScorer);
    bool verbose = false;
    searchICQuery(topN, *q, m_tol, verbose, annResults,indexNum); // index FAISS, use topN=100 to filter

    st.restart("mass diff analysis");
    bool do_massDiff_analysis = false;
    if (do_massDiff_analysis) {
        massdiff(*q, annResults);
    }
    st.restart("pairwise distance");
    vector<CAnnSpectra *> annVec;
    // and also among neighbors topN * topN
    // MZ file peaklist, 50 number per spec
    do_pairwise_distance(calcEdge, query_index, annVec, annResults); // real distance between query and topN

    
    if(m_savebackgroundscore){
// save the pairwise scores of close neighbors.
    string dpscore_neighbor_file = "dpscore_neighbor_query_" + to_string(query_index) + ".txt";
     ofstream fout (dpscore_neighbor_file.c_str(), ios::out);
    fout << "From\tTo\tDist\tDP" << endl;
    for (auto each : annVec){
        for(int i = 0; i < each->m_anns.size(); i ++)  {
            fout << each->getQueryIdx()
            << "\t" << each->m_anns[i].idx
            << "\t" << each->m_anns[i].dist
            << "\t" << each->m_anns[i].dotprod  << endl;
        }
     }
    }
   

    st.restart("pvlaue");
    bool do_pvalue_on_all = false;  // all: all the spectra in archive. NOT 10,000 background spectra
    if (do_pvalue_on_all) {
        cout << "[Warning] Never called" << endl;
        do_pvalue_onall(query_index, annVec);
    }
    bool plotfigure = false;  //20200806 generate figure
    SLinearRegressionModel modelLR;

    if (plotfigure)cout << "[Info] Figure will be generated " << endl;
    bool do_partial_score_calucation = true;
    if (do_partial_score_calucation) {
        //cout << "[Info] calcualting p-values on 10,000 spectra " << endl;
        do_pvalue_onfraction(annVec, *q, query_index, plotfigure, &dpscores, &modelLR);
    }

    bool do_neighbor_pvalue_calculation = false;
    if (do_neighbor_pvalue_calculation) {
        CAnnSpectra *R1 = annResults.get(0);
        do_neighbor_pvlaue_on_all(query_index, annVec, R1);
    }

    st.restart("visulization");
    if (visualize) {
        cout << "Running visualization step" << endl;
        visualization_of_score_distribution(query_index, annVec, calcEdge == 1, topN);
    }
    cout << "Start generating result in jsonstring" << endl;
    st.restart("response");
    // SQL DB used to annotated
    exportResponse(query_index, annVec, jsonstring, &dpscores, &modelLR);
    cout << "jsonstring size: " << jsonstring.size() / 1024 << "KB" << endl;
    delete q;
}


void CSpectralArchive::calcPvalue(int query_index, int tol, string outputbasename, bool normalize) {
    cout << "[Warning] This function will be removed soon; use CTailEstimation instead!" << endl;
    outputbasename = to_string(outputbasename, "_", "tol", tol, "normscore", normalize);

    vector<int> all_scores = m_pScorer->distributionAll(tol, query_index, normalize);

    exportTable(all_scores, outputbasename + ".allscore.txt", '\t');

    visualizeRawScores(outputbasename, all_scores, true);
    visualizeRawScores(outputbasename, all_scores, false);

    visualizeCDF(outputbasename, all_scores, true, true);
    visualizeCDF(outputbasename, all_scores, false, true);
}


void CSpectralArchive::visualizeCDF(string outputbasename, vector<int> &all_scores, bool noZero, bool logCDF) {
    double tailfraction = 0.1;
    vector<double> rx(all_scores.size(), 0), ry(all_scores.size(), 0);
    cout << "rx ry size " << rx.size() << "\t" << ry.size() << endl;
    int shift = noZero ? 1 : 0;
    double total = 0;
    for (size_t i = all_scores.size() - 1; i >= shift & i != size_t(-1); i--) {
        rx[i] = i;
        total += all_scores[i];
        if (i == all_scores.size() - 1) {
            ry[i] = all_scores[i];
        } else {
            ry[i] = all_scores[i] + ry[i + 1];
        }
    }
    cout << "total  " << total << endl;
    if (total > EPSILON) {
        for (int i = 0; i < all_scores.size(); i++) {
            ry[i] /= total;
        }
    }
    cout << "working on normalized result " << endl;
    CVisual::SImageSize imsize;
    CVisual::gnuplotWrapper info;
    string pngfilename = outputbasename + to_string("", "_", "nozero", noZero) + "_cdf_scores.png";
    info.set_filename(pngfilename)
            .set_height(imsize.m_hpx)
            .set_width(imsize.m_wpx);
    CVisual::gnuplot_curve_topng(rx, ry, "lines", info);
    exportTable(ry, outputbasename + "_cdf_allscores.txt", '\t');
    if (logCDF) {
        const int MAX_TOP50_COS = 42925;
        vector<double> rxPositive, logry;
        vector<vector<double>> rx_logry;
        for (int i = all_scores.size() - 1; i >= 0; i--) {
            if (ry[i] > 0) {
                if (rxPositive.size() == 0) {
                    rxPositive.push_back(rx[i] / MAX_TOP50_COS);
                    logry.push_back(log10(ry[i]));
                    rx_logry.push_back(vector<double>{rx[i] / MAX_TOP50_COS, log10(ry[i])});
                } else if (ry[i] - ry[i + 1] > EPSILON) {
                    rxPositive.push_back(rx[i] / MAX_TOP50_COS);
                    logry.push_back(log10(ry[i]));
                    rx_logry.push_back(vector<double>{rx[i] / MAX_TOP50_COS, log10(ry[i])});
                }
            }
            if (ry[i] > tailfraction) {
                break;
            }
        }
        cout << "working on log cdf result " << endl;
        string pngfilename = outputbasename + to_string("", "_", "nozero", noZero) + "_logpValue_scores.png";
        CVisual::SImageSize imsize;
        CVisual::gnuplotWrapper info;
        info.set_filename(pngfilename)
                .set_height(imsize.m_hpx)
                .set_width(imsize.m_wpx)
                .set_xlabel("normalized dot product")
                .set_ylabel("log_{10}(1-CDF(x))");

        CVisual::gnuplot_curve_topng(rxPositive, logry, "points pt 7 ", info);
        vector<double> coefs;
        double r2;
        bool success = linear_regression_on_logCDF(rx_logry, coefs, r2);
        if (success) {  // plot the image
            {
                Gnuplot gnuplot;
                // plot
                gnuplot << "set term pngcairo enhanced size 8in, 6in" << endl;
                gnuplot << "set output '" << outputbasename + "_rx_logcdf_allscores_adaptive_lr_model.png" << "'"
                        << endl;

                ostringstream osstitle;
                osstitle << setprecision(4) << "log(p)=" << coefs[0] << "x+" << coefs[1] << " R^2=" << r2;
                gnuplot << "f(x)=" << coefs[0] << "*x+" << coefs[1] << endl;

                gnuplot << "set xlabel 'dot product'" << endl;
                gnuplot << "set ylabel 'log_{10}(p-value)'" << endl;

                cout << __LINE__ << " ylabel to be fixed!!" << endl;
                gnuplot << "plot " << gnuplot.file1d(rx_logry)
                        << "  u 1:2 w points lc rgb'green' title 'p-value in log scale' ,"
                           " f(x) w line lc rgb'red' linetype '-' title '" << osstitle.str() << "'" << endl;

            }
            // plot another figure to text form
            {
                Gnuplot gnuplot;
                gnuplot << "set term dumb size 100, 40" << endl;

                ostringstream osstitle;
                osstitle << setprecision(4) << "log(p)=" << coefs[0] << "x+" << coefs[1] << " R^2=" << r2;
                gnuplot << "f(x)=" << coefs[0] << "*x+" << coefs[1] << endl;

                gnuplot << "set xlabel 'dot product'" << endl;
                gnuplot << "set ylabel 'log_{10}(p-value)'" << endl;
                gnuplot << "plot " << gnuplot.file1d(rx_logry)
                        << "  u 1:2 w points lc rgb'green' title 'p-value in log scale' ,"
                           " f(x) w line lc rgb'red' linetype '-' title '" << osstitle.str() << "'" << endl;
            }
        }
        exportTable(rx_logry, outputbasename + "_rx_logcdf_allscores.txt", '\t');
    }
}


void CSpectralArchive::visualizeRawScores(const string &outputbasename, vector<int> &all_scores, bool noZeros) const {
    int shift = noZeros ? 1 : 0;
    double total_score = accumulate(all_scores.begin() + shift, all_scores.end(), 0.0);
    vector<double> y(all_scores.begin() + shift, all_scores.end());
    for (int i = 0; i < y.size(); i++) y[i] /= total_score;
    vector<double> x(all_scores.size() - shift, 0);
    iota(x.begin(), x.end(), 1);
    CVisual::gnuplotWrapper info;
    CVisual::SImageSize imsize;
    string pngFilename = outputbasename + to_string("", "_", "nozero", noZeros) + "_scores.png";
    info.set_height(imsize.m_hpx)
            .set_width(imsize.m_wpx)
            .set_filename(pngFilename);
    CVisual::gnuplot_curve_topng(x, y, "lines", info);
}

struct IdxScore {
    long idx;
    int Score;

    IdxScore(long _idx, int _score) {
        idx = _idx;
        Score = _score;
    }

};

void
CSpectralArchive::getkTrueNearestNeighbor(ICQuery &query, vector<vector<long>> &topKTrueNN, bool isLowMassAcc, int topK,
                                          int dp_UInt) {
    // get True Nearest Neighbors from the spectral archive. 
    if (m_csa->getSpecNum() == 0) {
        cout << "[Error] mz file is empty: " << endl;
    }
    topKTrueNN.assign(query.getSpecNum(), vector<long>());
    Progress ps(query.getSpecNum(), "accurate dotproduct to all");
    for (int i = 0; i < query.getSpecNum(); i++) {
        ps.increase();
        if(query.getMzSpec(i).getPeakNum()<agtsummary.m_minPeakNumInSpec) continue;
        vector<long> allidx(m_csa->getSpecNum(), 0);
        vector<int> allscores(m_csa->getSpecNum(), 0);
        iota(allidx.begin(), allidx.end(), 0);
        if (allidx.empty()) {
            cout << "Empty index list " << endl;
        }
        m_csa->scorePartiallyWithVecForm(50, m_tol, 32, true, allidx, query.getSpecBy(i), allscores);
        // copy only positive numbers:
        vector<IdxScore> res;
        res.reserve(1000);
        for (long k = 0; k < allidx.size(); k++) {
            if (allscores[k] > dp_UInt) {
                res.emplace_back(allidx[k], allscores[k]);
            }
        }
        cout << endl << "number of candidates with dp score larger than " << dp_UInt << ": " << res.size() << endl;
        string tmp = to_string("candidates dp > ", " ", dp_UInt, res.size(), "\n");
        agtsummary.writeToLogFile(tmp);
        stable_sort(res.begin(), res.end(), [](const IdxScore &x, const IdxScore &y) { return x.Score > y.Score; });
        if (res.size() > topK) res.erase(res.begin() + topK, res.end());
//        res.resize(topK);
//        long idx = *std::min_element(allidx.begin(), allidx.end(),
//                                     [&](const long &x, const long &y) { return allscores[x] > allscores[y]; });
        cout << endl << "QueryID\tRank\tIdx\tScore\tDP" << endl;
        tmp = "QueryID\tRank\tIdx\tScore\tDP\n";
        agtsummary.writeToLogFile(tmp);
        for (int k = 0; k < res.size(); k++) {
            topKTrueNN[i].push_back(res[k].idx);
            double dp = res[k].Score * 1.0 / 42925;
            tmp = to_string("", "\t", query.getQueryIndex(i), k, res[k].idx, res[k].Score, dp, "\n");
            agtsummary.writeToLogFile(tmp);
            cout << query.getQueryIndex(i) << "\t" << k << "\t" << res[k].idx << "\t" << res[k].Score << "\t" << dp << endl;
        }

        //m_csa->calcDotProduct(50, m_tol, query.getSpecBy(i), 32, allidx, allscores);
//        vector<float> distscore(m_csa->getSpecNum(), 1.4);
//
//        for (long j = 0; j < m_csa->getSpecNum(); j++) {
//            float q = m_csa->getSquaredNorm(j);
//            float p = m_csa->getSquaredNorm(query.getSpecBy(i));
//            if (q > EPSILON and p > EPSILON) {
//                distscore[j] = sqrt(2.0) * sqrt(1.0 - allscores[j] / sqrt(p * q));
//            } else {
//                cout << "norm of q and p are too small " << endl;
//            }
//        }
//        long idx = *std::min_element(allidx.begin(), allidx.end(),
//                                     [&](const long &x, const long &y) { return distscore[x] < distscore[y]; });
//        topIdx[i] = idx;

    }


}

void CSpectralArchive::getAccurateTopNeighbor(ICQuery &query, vector<long> &topIdx, bool isLowMassAcc) {
    topIdx.assign(query.getSpecNum(), 0);
    Progress ps(query.getSpecNum(), "accurate dotproduct to all");
    if (m_csa->getSpecNum() == 0) {
        cout << "[Error] mz file is empty: " << endl;
    }
//    int bin_tol = isLowMassAcc? 15: 1;
    for (int i = 0; i < query.getSpecNum(); i++) {

        ps.increase();
        vector<long> allidx(m_csa->getSpecNum(), 0);
        vector<int> allscores(m_csa->getSpecNum(), 0);
        iota(allidx.begin(), allidx.end(), 0);
        if (allidx.empty()) {
            cout << "Empty index list " << endl;
        }
        m_csa->calcDotProduct(50, m_tol, query.getSpecBy(i), 32, allidx, allscores);
        vector<float> distscore(m_csa->getSpecNum(), 1.4);

        for (long j = 0; j < m_csa->getSpecNum(); j++) {
            float q = m_csa->getSquaredNorm(j);
            float p = m_csa->getSquaredNorm(query.getSpecBy(i));
            if (q > EPSILON and p > EPSILON) {
                distscore[j] = sqrt(2.0) * sqrt(1.0 - allscores[j] / sqrt(p * q));
            } else {
                cout << "norm of q and p are too small " << endl;
            }
        }
        long idx = *std::min_element(allidx.begin(), allidx.end(),
                                     [&](const long &x, const long &y) { return distscore[x] < distscore[y]; });
        topIdx[i] = idx;
        cout << i << "\t" << idx << endl;
    }
}

void CSpectralArchive::searchICQuery(ICQuery &query, int tolerance, bool verbose, CArxivSearchResult &archiveRes, int indexNum) {
    if (m_indices == nullptr) {
        cout << "Error: Empty index!" << endl;
        throw string("Invalid index in !") + string(__FUNCTION__) + ":" + to_string(__LINE__);
    }

    int ret_num = 1024;


    int indexChoice = -1;
    vector<vector<long>> results;
    int querySize = query.size();
    cout << "changing index num to " << indexNum << endl;
    m_indices->getAnns(query, ret_num, results,indexNum);
//    cout << querySize << "\t" << results.size() << endl;

    if (agtsummary.m_recallOfTrueNeighbor) {
        string tmp = to_string("nprobe","\t",agtsummary.m_nprobe,"\n");
        agtsummary.writeToLogFile(tmp);
        vector<vector<long>> topKTrueNN;
        const int MAX_SCORE = 42925;
        double mindp = agtsummary.m_recallTNNminDP;
        int minDpScoreInt = mindp * MAX_SCORE;
        int topK = agtsummary.m_recallTNNtopK;
        bool isLowMassAcc = false;
        getkTrueNearestNeighbor(query, topKTrueNN, isLowMassAcc, topK, minDpScoreInt);
        if (topKTrueNN.size() == results.size()) {// WHY!!!
            if(verbose) cout << "[OK] same size for top k true neighbors, and top k approx neighbors " << endl;
        }

        for (int j = 0; j < querySize; j++) {
            if(query.getMzSpec(j).getPeakNum()<agtsummary.m_minPeakNumInSpec){
                continue;
            }
            vector<long> collectIdx;
            m_indices->collectANNs(indexChoice, ret_num, results, j, collectIdx, verbose);
            long queryIndex = query.getQueryIndex(j);
            string tmp = "rank\tTruthID\tFound\tReportID\tRecall@rank\tquery\t" + to_string(queryIndex) + "\n";
            cout << "rank\tTruthID\tFound\tReportID\tRecall@rank\tquery\t" << queryIndex << endl;
            map<long, int> z;  
            // mapping id to its last occurence 
            // therefore we got list of unique keys. id.
            for (int k = 0; k < collectIdx.size(); k++) {
                z[collectIdx[k]] = k;
            }
            // stats
            int numOfTNN = topKTrueNN[j].size();
            int numberOfTNNsFound = 0;
            long queryindex = query.getQueryIndex(j);
            for (int k = 0; k < topKTrueNN[j].size(); k++) {
                int idx = -1;
                // bool isfound = (z.count(topKTrueNN[j][k]) == 0) ? false : (idx = z[topKTrueNN[j][k]], true);
                bool isfound = false;
                if (z.count(topKTrueNN[j][k]) == 0){
                    // k-th TNN not found.
                    isfound = false;
                }else{
                    // k-th TNN found.
                    idx = z[topKTrueNN[j][k]];
                    isfound = true;
                    numberOfTNNsFound += 1;
                }


                tmp = to_string(tmp, "\t", k, topKTrueNN[j][k], boolalpha, isfound,(isfound ? collectIdx[idx] : idx), numberOfTNNsFound*1.0/(k+1), "\n");
                cout << k << "\t" << topKTrueNN[j][k] << boolalpha << "\t" << isfound << "\t"
                     << (isfound ? collectIdx[idx] : idx)<< "\t" << numberOfTNNsFound*1.0/(k+1)<< endl;

                
                agtsummary.increase(isfound, queryindex);
            }
            agtsummary.writeToLogFile(tmp);
            tmp = to_string("","\t","QueryIndex",queryIndex, 
            "#TNN", numOfTNN, "#TNNFound",  numberOfTNNsFound, 
            "oneQueryRecall", numberOfTNNsFound*1.0/numOfTNN, "nprobe", agtsummary.m_nprobe, 
            "indexNum", m_indices->getNumOfindexInUse(),"mindp", agtsummary.m_recallTNNminDP, "topK", agtsummary.m_recallTNNtopK, "\n");
            agtsummary.writeToLogFile(tmp);
            cout << tmp << flush;
        }
        agtsummary.print();
    }


    vector<vector<long>> allRetIdx(querySize, vector<long>());
//    verbose =true;
    for (int i = 0; i < querySize; i++) {
        m_indices->collectANNs(indexChoice, ret_num, results, i, allRetIdx[i], verbose);
        filterWithMinPeakNum(verbose, allRetIdx[i]);
        m_AnnotationDB->filterwithblacklist(verbose, allRetIdx[i]);
    }
//    search.restart("dp score");
    vector<vector<float>> accDist(querySize, vector<float>());
    vector<vector<int>> dpscores(querySize,vector<int>());
    int threadnum = getProperThreads() / 5 * 2 + 1; // 66% percent, but now zero by add 1!
    m_pScorer->dpscore(tolerance, allRetIdx, threadnum, accDist, query,dpscores);
    archiveRes.init(query, verbose, allRetIdx, accDist,dpscores);
}

// Find n neighbors via index
// calculate real distances.
// keep top N neighbors
void
CSpectralArchive::searchICQuery(int topN, ICQuery &query, int tolerance, bool verbose, CArxivSearchResult &archiveRes, int indexNum) {
    searchICQuery(query, tolerance, verbose, archiveRes,indexNum);
    archiveRes.keepTopN(topN, verbose);
}

void CSpectralArchive::filterWithMinPeakNum(bool verbose, vector<long> &retIdx) {
    vector<long> retIdx_tmp;
    int counts = 0;
    for (int k = 0; k < retIdx.size(); k++) {
        if (m_pScorer->getPeakNum(retIdx[k]) > m_minPeakNum) {
            retIdx_tmp.push_back(retIdx[k]);
        } else {
            counts++;
        }
    }
    retIdx.swap(retIdx_tmp);
    if (verbose)cout << "Number of spectra removed (peak num < " << m_minPeakNum << "): " << counts << endl;
}

void CSpectralArchive::exportResponse(int query_index, vector<CAnnSpectra *> &annSpecList, string &jsonstring,
                                      vector<int> *ptrDPs, SLinearRegressionModel *ptrLRModel) {
    CArchiveSearchReply asr(annSpecList, query_index);
    asr.tojsonstring(jsonstring, *m_AnnotationDB);
//    cout << "debug: Got Jsonstring: pointers: DPsPtr " << ptrDPs << " lrModelPtr " << ptrLRModel << endl;
    if (not ptrDPs->empty()) {
        // if not empty: let's reply with DP scores as string
        ostringstream oss;
        oss << "\"backgroundscores\":[";
        for (int i = 0; i < ptrDPs->size(); i++) {
            if (i > 0) {
                oss << ", ";
            }
            oss << ptrDPs->at(i);
        }
        oss << "]";
        jsonstring = jsonstring.substr(0, jsonstring.length() - 3) + "," + oss.str() + "\n}\n\n";
//        cout << "debug: dps done!" << endl;
    }

    if (ptrLRModel != nullptr) {
        jsonstring = jsonstring.substr(0, jsonstring.length() - 3) + "," + ptrLRModel->toJsonString() + "\n}\n";
    }
//    cout << "debug: json done" << endl;
}

void CSpectralArchive::do_pvalue_onall(long query_index, vector<CAnnSpectra *> &vqrSimple) {
    cout << "[Warning] Should not be called any more! Not efficient! calculating pvalue with all samples..." << endl;
    vector<int> scores = m_pScorer->distributionAll(1, query_index, true);
    shared_ptr<ICHistdata> histscore = make_shared<CHistInt>(scores);// 20200806
    CTailEstimation cte(histscore, true);
    bool good = cte.buildLinearRegressionModel(true, to_string(query_index), true, false);
    if (not good) return;
    for (int k = 0; k < vqrSimple.size(); k++) {
        for (int m = 0; m < vqrSimple[k]->m_anns.size(); m++) {
            double x = vqrSimple[k]->m_anns[m].dist;
            double dp = 1 - x * x / 2;
            vqrSimple[k]->m_anns[m].pvalueAll = cte.getlogpvalue(dp);
        }
    }
}

void CSpectralArchive::do_pvalue_onfraction(vector<CAnnSpectra *> &vqrSimple, ICQuery &q, long query_index,
                                            bool plot_figure, vector<int> *ptrDPs, SLinearRegressionModel *ptrLR) {
    SimpleTimer st("p-value shuffle");
    cout << "[Info] number of ANNs: " << vqrSimple.size() << endl;
    if (vqrSimple.empty() or vqrSimple[0]->m_anns.empty()) {
        cout << "Empty result " << endl;
        return;
    }
    CArxivSearchResult queryResults;
    for (auto each: vqrSimple) queryResults.push_back(each);
    vector<long> idxlist = {query_index};

    bool verbose_Pvalue_calculation = false; // takes time when set to true
    bool saveBackgroundScore = false; //
    m_pc->run(idxlist, q, queryResults, plot_figure, verbose_Pvalue_calculation, ptrDPs, ptrLR, nullptr, nullptr);

    // release objects
    for (int i = 0; i < queryResults.size(); i++) {
        queryResults.set(i, nullptr, false);
    }
}

void
CSpectralArchive::do_neighbor_pvlaue_on_all(long query_index, vector<CAnnSpectra *> &annSpecList, CAnnSpectra *R1) {
    // We tried to use the best approximation to calculate the model; but it is not very accurate
    // run the following command to visulize the difference;
    //
    // gnuplot -e 'set xrange [0:1]; set yrange [-20:0]; p -18.88*x + 3.068 title "truth", -17.5*x+1.415 title "estimation"'
    // gnuplot -e 'set xrange [0:1]; set yrange [-20:0]; p -19.59*x+3.979 title "truth", -20.78*x + 4.058 title "estimation"'

    cout << "pvalue of neighbor" << endl;
    long neighbor_query = R1->m_anns[0].idx == query_index ? R1->m_anns[1].idx : R1->m_anns[0].idx;
    cout << "--------neighbor: " << neighbor_query << endl;

    vector<int> scores_of_neighbor = m_pScorer->distributionAll(1, neighbor_query, true);
    shared_ptr<ICHistdata> histscore = make_shared<CHistInt>(scores_of_neighbor);//20200806
    CTailEstimation cte_of_neighbor(histscore, true);
    bool good = cte_of_neighbor.buildLinearRegressionModel(true, to_string(query_index), true, false);
    if (not good) return;
    for (int k = 0; k < annSpecList.size(); k++) {
        for (int m = 0; m < annSpecList[k]->m_anns.size(); m++) {
            double x = annSpecList[k]->m_anns[m].dist;
            double dp = 1 - x * x / 2;
            annSpecList[k]->m_anns[m].pvalueAll = cte_of_neighbor.getlogpvalue(dp); // to fix
            annSpecList[k]->m_anns[m].neighborPvalueAll = true;
        }
    }
    cout << "neighbor query done" << endl;
}

void CSpectralArchive::visualization_of_score_distribution(long query_index, vector<CAnnSpectra *> &annSpecList,
                                                           bool tworound, int topN) {
    cout << "[Warning] This function should not be called any more" << endl;
    vector<double> scores;
    for (auto &each: annSpecList) {
        for (auto &ann: each->m_anns) {
            scores.push_back(ann.dist);
        }
    }

    vector<Distribution *> d;
    d.push_back(static_cast<Distribution *> (new Gaussian(-1, 1.001)));// todo: new/delete!
    EM em(d, scores);
    em.run(50);
    em.D[0]->visualize();
    em.visualizeSample();

    int peaknum = m_pScorer->getMzSpec(query_index).getPeakNum();

    SPsmAnnotation gtinfo;
    m_AnnotationDB->retrieveGtinfo(query_index, gtinfo);

    string outputbasename = to_string("", "_", query_index, "peaknum", peaknum, "rnd", tworound, "topn", topN,
                                      "mu_sigma_", em.getNames(), "p", gtinfo.precursormass, "chg", gtinfo.charge);

    double minval = 0, maxval = 1.4142;
    CVisual::gnuplotWrapper info;
    CVisual::SImageSize imsize;
    info.set_filename(outputbasename
                      + to_string(m_tol) + "_query_kNNScores.png")
            .set_minmax(minval, maxval)
            .set_width(imsize.m_wpx)
            .set_height(imsize.m_hpx);
    CVisual::gnuplot_histogram_topng(scores, 100, minval,
                                     maxval, info);
    bool normalize = true;
    calcPvalue(query_index, 15, outputbasename, normalize);
    calcPvalue(query_index, 1, outputbasename, normalize);
}

//// todo: NEVER write size()-1
// Calculate edge/distance from source to targets
// input: a list of n ids
// output: pair-wise connections; a object with n(n-1)/2 (that is, n-choose-2) edges.
void CSpectralArchive::do_pairwise_distance(int calcEdge, long query_index, vector<CAnnSpectra *> &vecAnns,
                                            CArxivSearchResult &annResults) {
    vector<long> annIds = annResults.concatenateIndex(0);
    CAnnSpectra *queryAnns = annResults.get(0);
    if (calcEdge == 1) {
        cout << "[Info] Size of ANNs: " << annIds.size() << endl;
        for (int i = 0; i + 1 < annIds.size(); i++) {
            // calculate distance between source node (sourceId) to target nodes (targetIds).
            long sourceId = annIds.at(i);
            vector<long> targetIds(annIds.begin() + i + 1, annIds.end());
            vector<float> distance;
            vector<int> dpscore;
            m_pScorer->dist(sourceId, targetIds, m_tol, distance,dpscore);
//            cout << "size of dpscore and size of distance score " << dpscore.size() << " " << distance.size() << endl;
//            cout << "First dp and dist: " << dpscore.at(0) << "\t" << distance.at(0) << endl;

            CAnnSpectra *p = new CAnnSpectra(targetIds, distance, sourceId,dpscore);
            if (p == nullptr) { cout << "[Error] Fail to create object: " << __FUNCTION__ << endl; }
            vecAnns.push_back(p);
        }

        if (not queryAnns->isExist(query_index)) {
            vecAnns.push_back(queryAnns);
            annResults.set(0, nullptr, false); // moved!
        }
    } else {
        vecAnns.push_back(queryAnns);
        annResults.set(0, nullptr, false); // moved!
    }
    cout << "[Info] Size of edge clusters (tree-ish object): " << vecAnns.size() << endl;
}

void CSpectralArchive::searchMzFileInBatch(CMzFileReader &querySpectra, long first, long last, string validationfile,
                                           int topHitsNum, int numOfpValues, bool recalltrueneighbor, int batchSize,
                                           int bgspec_seed, int recallTNNtopK, double recallTNNminDP,
                                           bool skipBackgroundScoreCalc, bool useFlankingBins) {
    // agtsummary.setRecallTNN(recalltrueneighbor, recallTNNtopK, recallTNNminDP, recallTNNminPeakNumInSpec);
    agtsummary.setvalidationfile(validationfile);
    string JsonFile = querySpectra.getMzFilename() + to_string("_","_",first,last) +"_out.json";
    string allJsonStrs ="{\n";
    int initSize = allJsonStrs.size();
    SimpleTimer st("Search Mz file");
    if (last < 0 or last > querySpectra.getSpecNum()) {
        last = querySpectra.getSpecNum();
        spdlog::get("A")->info("Warning: resetting range [first, last] -> [{}, {}], because of the total number of querySpectra is {}", first, last, last);
    }

    if (first < 0 or first > querySpectra.getSpecNum() or first >= last) {
        spdlog::get("A")->info("Warning: the range [first, last] -> [{}, {}] is invlaid. ", first, last);
        first = 0;
        spdlog::get("A")->info("Warning: reset the range [first, last] as [{}, {}]. ", first, last);
    }
//    bool verbose = false;
    if(m_verbose)cout << "range to be processed: " << first << "\t" << last << endl;
    int batchsize = batchSize;

    int threadnum = 10;
    long ntotalnum = m_pScorer->getSpecNum();

    shared_ptr<CPvalueMultiCalculator> pvalueMulticalc=nullptr;
    if(not skipBackgroundScoreCalc){
        pvalueMulticalc = make_shared<CPvalueMultiCalculator>(numOfpValues, ntotalnum, threadnum, m_pScorer, bgspec_seed, m_tol);
    }
    vector<CAnnSpectra *> vqrs(last - first, nullptr);
    if(m_verbose) cout << "the data file is : " << querySpectra.getListFile() << endl;
    m_AnnotationDB->createBlackListWithFileName(querySpectra.getListFile(), m_verbose);

    Progress ps( last-first, "Search");
    for (long i = first; i < last; i += batchsize) {
        long start = i, end = i + batchsize < last ? i + batchsize : last;
        long queryNum = end -start;
        if(m_verbose)cout << "processing : " << start << "\t" << end << endl;
        vector<long> idxlist(queryNum, 0);
        iota(idxlist.begin(), idxlist.end(), i);
        CMzQuery query(querySpectra, m_dim, idxlist, useFlankingBins);
        if(m_verbose)cout << "Size of query: " << query.size() << endl;
        CArxivSearchResult searchRes;
        searchICQuery(topHitsNum, query, m_tol, m_verbose, searchRes,-1);
        if(m_verbose)cout << "Searching in archive done"<< endl;
        vector<vector<int>> ptrDPs(numOfpValues, vector<int>());
        vector<vector<vector<int>>> ptrDP_all_pv_query(numOfpValues, vector<vector<int>>(idxlist.size(), vector<int>()));
        vector<SLinearRegressionModel> slRMs(numOfpValues, SLinearRegressionModel());
        vector<vector<SLinearRegressionModel>> slRMs_pv_query(numOfpValues, vector<SLinearRegressionModel>(idxlist.size(), SLinearRegressionModel()));
        if(pvalueMulticalc){
            if(m_verbose)cout << "start pvalue calculation"<< endl;
            pvalueMulticalc->run(idxlist, query, searchRes, false, m_verbose, &ptrDPs, &slRMs, &ptrDP_all_pv_query, &slRMs_pv_query);
            if(m_verbose)cout << "finish pvalue calcualtion"<< endl;
        }


        for(int j = 0; j < queryNum; j ++){
            ps.increase();
            vector<CAnnSpectra *> localvqr{searchRes.get(j)};
            string jsonString ;
            exportResponse(idxlist[j],localvqr,jsonString, &ptrDP_all_pv_query[0][j], &slRMs_pv_query[0][j]);
            if(allJsonStrs.size() > initSize) {allJsonStrs += ",\n";}
            allJsonStrs += "\""+to_string(idxlist[j]) + "\" : " + jsonString ;
        }
        searchRes.moveTo(vqrs, start - first);

    }
    allJsonStrs += "\n}\n";
    ofstream fout(JsonFile.c_str(), ios::out);
    if(fout.is_open()){
        fout << allJsonStrs << endl;
        fout.close();
        if(m_verbose) cout << "Json file saved as " << JsonFile << endl;
    }


    spdlog::get("A")->info("Retrieving annotation from database...");
    CMzFileSearchResults mzSR(*m_AnnotationDB, topHitsNum);
    mzSR.retrieveAnnotation(vqrs);

    spdlog::get("A")->info("Exporting result into tsv file");
    string outfile = to_string(querySpectra.getListFile(), "_", first, last, ".psmtsv");
    mzSR.exportTsv(outfile);

    if (not validationfile.empty() and File::isExist(validationfile)) {
        mzSR.validation(validationfile, querySpectra);
    }
}

void CSpectralArchive::massdiff(ICQuery &q, CArxivSearchResult &annRes) {
    vector<vector<int>> prescore;
    m_csa->get_prescoreMatrix(prescore);
    for (int m = 0; m < annRes.get(0)->m_anns.size(); m++) {
        SAnnSpectrum &each = annRes.get(0)->m_anns[m];
        SPsmAnnotation gtinfo;
        double massA = 0;
        if (q.getQueryIndex(0) != -1) {
            m_AnnotationDB->retrieveGtinfo(q.getQueryIndex(0), gtinfo);
            massA = gtinfo.precursorNeutralMass();
        }
        m_AnnotationDB->retrieveGtinfo(each.idx, gtinfo);
        double massB = gtinfo.precursorNeutralMass();
        vector<double> diff;
        cout << q.getQueryIndex(0) << " -- " << each.idx << "----------delta mass of precursor: " << massA - massB
             << endl;
        uint16_t *x = q.getPtrUint16(0);
        uint16_t *y = m_pScorer->getSpecBy(each.idx);

        bool print_matched_pks = false;
        bool debug = false;
        long score = m_csa->calculate_dot_product_with_hist(m_tol, 50, prescore, diff, x, y, print_matched_pks);
    }
}

void CSpectralArchive::addRemark(long query_index, string &remarks) {
    m_AnnotationDB->updateSpecRemarks(query_index, remarks);
}
// todo
void CSpectralArchive::getRemark(long query_index, string &remarks) {
    vector<vector<string>>  remark_vec;
    m_AnnotationDB->searchRemarks(query_index, remark_vec);
    if(not remark_vec.empty()){
        remarks = remark_vec[0][2];
    }

}


void CSpectralArchive::getScorerFactoryPtr() {
//#ifndef __CUDA__
//    if(m_usegpu)
//    {
//        m_usegpu = false;
//        cout << "[Warning] GPU not available, will keep using CPU" << endl;
//    }
//#endif
    if (m_usegpu) {
#ifdef __CUDA__
        m_scorerFactory = make_shared<CMzCUDAFactory>();
#else
        m_usegpu = false;
        cout << "GPU is not available. Switch to use CPU" << endl;
        m_scorerFactory = make_shared<CMzFileReaderFactory>();
#endif
    } else {
        m_scorerFactory = make_shared<CMzFileReaderFactory>();
    }
    m_scorerFactory->setfilename(m_mzXMLListFileName + ".mz");
    m_scorerFactory->setMzFileObject(m_csa);
    m_pScorer = m_scorerFactory->create();

}

string CSpectralArchive::getPlatform() {
    return m_usegpu ? "GPU" : "CPU";
}

void CSpectralArchive::addpepxmlfile(string pepxmllist, string new_gt_file) {
    if (not File::isExist(pepxmllist,true)) {
//        cout << "pepxmllist file not found: " << pepxmllist << endl;
//        cout << "[Error] pepxml file does not exist! Doing Nothing!!! Function to be removed " << __LINE__ << __FILE__ << __FUNCTION__ << endl;
        return;
    }
    CTable pepxml(pepxmllist, '\t', false, 0);

    int rownum = pepxml.getRowByKey(new_gt_file, 0);

    if (rownum != -1) {
        cout << "ground truth file already exist in line: " << rownum << endl;
        cout << "Line " << rownum << ": " << pepxml.getEntry(rownum, 0) << endl;
    } else {
        cout << "ground truth file list will be updated!" << endl;
        pepxml.addRow({new_gt_file});
        pepxml.saveAs(pepxmllist, false, '\t');
    }
}


string convertChromatogramToSVG(shared_ptr<Chromatogram> chr, string filename, string output_title) {
    vector<vector<double>> x;
    chr->toVector(x, 200);
    cout << "x.size(): " << x.size() << endl;

    string outfile = "tmp.svg";
    {
        Gnuplot gnuplot;
        gnuplot << "set terminal svg" << endl;
        gnuplot << "set output 'tmp.svg'" << endl;
        gnuplot << "unset key" << endl;
        if (output_title != "") {
            gnuplot << "set title '" << output_title << "' noenhanced" << endl;
        } else {
            gnuplot << "set title '" << File::CFile(filename).filename << " scan=" << chr->m_ms2Scan << "' noenhanced"
                    << endl;

        }
        gnuplot << "set xlabel ' RT(sec) ' " << endl;
        gnuplot << "set ylabel ' Intensity ' " << endl;

        gnuplot << "set arrow from " << int(chr->m_ms2_rt_in_second) << ", graph(0,0) to "
                << int(chr->m_ms2_rt_in_second) << ", graph(1,1) nohead dt '.'" << endl;
        gnuplot << "plot '-' with lines lt 1 pt 1" << endl;
        gnuplot.send1d(x);
    }

    return readlinesfromfile(outfile);

}

void CSpectralArchive::getRawSpec(long queryindex, vector<double> &mz, vector<double> &intensity) {
    SPsmAnnotation gtinfo;
    m_AnnotationDB->retrieveGtinfo(queryindex, gtinfo);
    int scan_num = gtinfo.ms2_scan;
    string filename = m_AnnotationDB->getSpecFileNameByID(gtinfo.fileid);
    cout << "got filename and scan num  " << filename << " at " << scan_num << endl;
    // if file does not exist, append the prefix path of m_mzXMLFilename
    if(not File::isExist(filename, true)){
        // try add prefix.
        filename = File::CFile(m_mzXMLListFileName).path + "/" + filename;
        cout << "filename is fixed: " << filename <<endl;
    }
    
    // if file does not exist, append the prefix path of m_mzXMLFilename
    bool status = getPeakList(filename, scan_num, mz, intensity);
    cout << "got the real peaks, with status " << status << endl;
}

void CSpectralArchive::getChromatograms(long queryindex, bool getChromatogram, bool updateChromSVG, double mzTolerance,
                                        int charge, map<int, shared_ptr<Chromatogram>> &ChromatogramOfIsotopes,
                                        string &mzrtmap) {
    SPsmAnnotation gtinfo;
    m_AnnotationDB->retrieveGtinfo(queryindex, gtinfo);
    int scan_num = gtinfo.ms2_scan;
    string filename = m_AnnotationDB->getSpecFileNameByID(gtinfo.fileid);
    cout << "got filename and scan num  " << filename << " at " << scan_num << endl;

    if (getChromatogram and File::CFile(filename).ext != "sptxt" and File::CFile(filename).ext != "mgf") {
        DataFile df(filename, 0, -1);
        double delta_mass = 1.0032;
        for (auto it = ChromatogramOfIsotopes.begin(); it != ChromatogramOfIsotopes.end(); it++) {
            int key = it->first;
            double mzOffset = key * delta_mass / charge;
            it->second = df.extractChromatogram(mzTolerance, scan_num, mzOffset, charge);
            shared_ptr<Chromatogram> chr = it->second;
            if (updateChromSVG) {
                string output_title = to_string(File::CFile(filename).filename, " ", " scan =", chr->m_ms2Scan,
                                                "Mz =", std::fixed, std::setprecision(3), chr->m_precursorMz,
                                                "Th queryID =", queryindex);
                string svgStr = convertChromatogramToSVG(chr, filename, output_title);
                cout << "svgStr: " << svgStr << endl;

            }
        }

        CMzRtMap m1(df);
        double precursorMz = ChromatogramOfIsotopes[0]->m_precursorMz, rt = ChromatogramOfIsotopes[0]->m_ms2_rt_in_second;
        CMzRtMap::projectionParam p{precursorMz - 1.1, precursorMz + 1.6, rt - 60, rt + 60, 100, 100};
        m1.projection(p);
        mzrtmap = m1.toJsonStr();
    }
}

// add neighbor info to db.
void CSpectralArchive::searchNeighborsWithin(double min_dp, long start, long end) {
    cout << "Searching with range " << start << " " << end << endl;
    int batchSize = 100;
    if (end == -1 or end > size()) end = size();
    if(start < 0 or start>end) start = 0;
    long taskNum = end-start;
    Progress ps(taskNum, "Search for Neighbors");

    for (long i = start; i < end; i += batchSize) {
        CArxivSearchResult archiveRes;
        int len = i + batchSize < end ? batchSize : end - i;
        vector<long> idxlist(len, 0);
        iota(idxlist.begin(), idxlist.end(), i);

        cout << "processing: " << idxlist.front()<< " " << idxlist.back() << endl;
//        exit(-1);

        ICQuery *query = createICQuery(nullptr, &idxlist, true, getDim(), *m_pScorer);
        searchICQuery(*query, m_tol, false, archiveRes,-1);
        archiveRes.keepWithMinDP(min_dp);
        // save the results to database
        ps.increase(len);
//        cout << "add db" << endl;
        m_AnnotationDB->addNeighborInfo(archiveRes);
        cout << "neighbor info added to database " << endl;
    }
    cout << "Done with queries from archive range " << start << "\t" << end << endl;

}



struct SPeptideSet {
    set<string> missed;
    set<string> insig;

    void printSummary() {
        cout << "Peptide (annotated by another search engine) NOT find correct answer in top N Hits : " << missed.size()
             << endl;
        for (auto each: missed) cout << each << endl;
        cout << endl;

        cout << "Peptide (annotatted by another search engine) FOUND but lower rank, not significant!! : "
             << insig.size() << endl;
        for (auto each: insig) cout << each << endl;
    }
};

struct SVennStats {
    int gtSigNum;
    int outOfDomainNum;
    int archiveSigOnly;
    int mix_InDomain_DiffSig;
    int outOfDomain_SigNum;
    int intxn_num;
    int intxn_SigNum;
    int intxn_Rank1Num;
    int intxn_Rank1SigNum;

    SVennStats() {
        gtSigNum = 0;
        intxn_num = 0;
        intxn_SigNum = 0;
        outOfDomainNum = 0;
        archiveSigOnly = 0;
        outOfDomain_SigNum = 0;
        intxn_Rank1Num = 0;
        intxn_Rank1SigNum = 0;
        mix_InDomain_DiffSig = 0;
    }

    void updateVenn(bool found, bool significantGT, bool sigArchive, bool inSpace) {
        CANSIConsole a;
        if (not found or not significantGT) {
            if (sigArchive) {
                cout << a.getColorStr("archiveSigOnly", CANSIConsole::GREEN) << endl;
                archiveSigOnly++;
            }
        } else {
            gtSigNum++;
            if (not inSpace) {
                cout << "NOT in SPACE of archive search" << endl;
                outOfDomainNum++;
                if (sigArchive) {
                    cout << a.getColorStr("Mixture", CANSIConsole::GREEN) << endl;
                    outOfDomain_SigNum++;
                }
            }
        }
    }

    void printSummary() {
        spdlog::get("A")->info("Summary of the venn validation step");
        spdlog::get("A")->info("#PSM significant in groundtruth {} ", gtSigNum);
        spdlog::get("A")->info("#PSM significant but not in archive space {} ", outOfDomainNum);
        spdlog::get("A")->info("#PSM significant and in space: {}", gtSigNum - outOfDomainNum);

        double recall = intxn_SigNum * 1.0 / (gtSigNum - outOfDomainNum);
        spdlog::get("A")->info("Intersection(without FDR filter): {}; Recall: {} / ({}-{}) = {}", intxn_num,
                               intxn_SigNum, gtSigNum, outOfDomainNum, recall);
        double rank1recall =
                intxn_Rank1SigNum * 1.0 / (gtSigNum - outOfDomainNum);
        spdlog::get("A")->info("Intersection(without FDR filter): {}; Recall: {} / ({}-{}) = {}", intxn_num,
                               intxn_Rank1SigNum, gtSigNum, outOfDomainNum,
                               rank1recall);

        cout << "---summary table---" << endl;
        string tmp = to_string("table\t", "\t", "Rank1", "Rank2-N", "Sum\n");
        tmp = to_string(tmp, "\t", "FDR<1%", intxn_Rank1SigNum,
                        intxn_SigNum - intxn_Rank1SigNum,
                        intxn_SigNum, "\n");
        tmp = to_string(tmp, "\t", "FDR>1%", intxn_Rank1Num - intxn_Rank1SigNum,
                        (intxn_num - intxn_SigNum) - (intxn_Rank1Num - intxn_Rank1SigNum), intxn_num - intxn_SigNum,
                        "\n");
        tmp = to_string(tmp, "\t", "Sum", intxn_Rank1Num, intxn_num - intxn_Rank1Num, intxn_num, "\n");
//        cout << tmp << endl;
        spdlog::get("A")->info("summary table: \n {}", tmp);

        double up_ratio = archiveSigOnly * 1.0 / intxn_SigNum;// intxn_num;
        spdlog::get("A")->info("Spectral Archive search only psm num: {} / {} = {}", archiveSigOnly,
                               intxn_SigNum, up_ratio);

        double up_ratio_of_gt = (gtSigNum - outOfDomainNum - intxn_SigNum) / 1.0 / intxn_SigNum;
        spdlog::get("A")->info("Groundtruth only psm num: {} / {} = {}",
                               gtSigNum - outOfDomainNum - intxn_SigNum,
                               intxn_SigNum, up_ratio_of_gt);

        spdlog::get("A")->info("Mixture: not in space but significant: {} in space but different: {}",
                               outOfDomain_SigNum, mix_InDomain_DiffSig);
    }

};
// new lines add

class ICExporter {
public:
    virtual ~ICExporter() {}

    virtual void processing(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo,
                            float dist, float evalue, float m_threshold_pvalue, PSMInfo &searchpsminfo,
                            string &prevstring) = 0;

    virtual void
    processingZeroRanking(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo, float dist,
                          float evalue, float m_threshold_pvalue, PSMInfo &searchpsminfo, string &prevstring) = 0;

    virtual void outputstr(string &str) = 0;
};

class ColorExporter : public ICExporter {
    SVennStats *m_venn;
    SPeptideSet *m_pepSet;
    CANSIConsole a;
public:
    ColorExporter(SVennStats *venn, SPeptideSet *pepSet) {
        m_venn = venn;
        m_pepSet = pepSet;
    }

    void outputstr(string &str) override {
        cout << str << endl;
    }

    void processingZeroRanking(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo, float dist,
                               float evalue, float m_threshold_pvalue, PSMInfo &searchpsminfo,
                               string &prevstring) override {
        // in this function, we still need to fix the following function
        // tophit.getProtein_UseAlterProteinIfItsNotDecoy()
        string tmp = to_string("NOT_FOUND(CONSISTENT TOPHIT): ", "\t", 0, scaninfo, "truth",
                               tophit.m_modified_peptide,
                               tophit.getProtein_UseAlterProteinIfItsNotDecoy(), std::defaultfloat, tophit.m_expect,
                               std::defaultfloat,
                               std::fixed, searchpsminfo.precursorNeutMass, tophit.m_peptideprophet_prob,
                               tophit.m_iprophet_prob);
        cout << a.getColorStr(tmp, CANSIConsole::RED) << endl;
        cout
                << "----------------------------------------------ARCHIVE FAIL-------------------------------------------------"
                << endl;
        m_pepSet->missed.insert(tophit.m_modified_peptide);
    }

    void
    processing(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo, float dist, float evalue,
               float m_threshold_pvalue, PSMInfo &searchpsminfo, string &prevstring) override {
        // part L : low rank ones but find consistent
        m_venn->intxn_num++;
        string tmp = "";
        if (ranking == 1) {
            m_venn->intxn_Rank1Num++;
            tmp = to_string("CORRECT(TOP)/CONSISTENT:\t", "\t", ranking, std::fixed, setprecision(4),
                            scaninfo, scaninfo.getPrecursorNeutralMass(), gtinfo->getModifiedPeptide(false), "DIST",
                            dist);
        } else {
            tmp = to_string("CORRECT(LOWRANK)/CONSISTENT:\t", "\t", ranking, std::fixed, setprecision(4),
                            scaninfo, scaninfo.getPrecursorNeutralMass(), gtinfo->getModifiedPeptide(false), "DIST",
                            dist);
        }
        tmp = a.getColorStr(tmp, CANSIConsole::BLUE);

        if (evalue < m_threshold_pvalue) {
            // part M:  significant
            tmp = tmp + "\t" + a.getColorStr(to_string("", "\t", "EVALUE(Sig:GOOD)", std::scientific, evalue, "DIST",
                                                       dist),
                                             CANSIConsole::GREEN);
            m_venn->intxn_SigNum++;
            if (ranking == 1) {
                m_venn->intxn_Rank1SigNum++;
            }
        } else {
            // intersection but not significant under BH
            tmp = tmp + "\t" + a.getColorStr(to_string("", "\t", "EVALUE(Insig:BAD)", std::scientific, evalue , "DIST",
                                                       dist),
                                             CANSIConsole::RED);
            m_pepSet->insig.insert(tophit.m_modified_peptide);
        }
        // print some information
        tmp = to_string(tmp + "\t", "\t", tophit.m_expect, std::defaultfloat,
                        tophit.m_modified_peptide, std::fixed,
                        searchpsminfo.precursorNeutMass, tophit.m_peptideprophet_prob,
                        tophit.m_iprophet_prob, "\n");
        cout << tmp << endl;
    }

};

class FileExporter : public ICExporter {
    ofstream m_fout;
public:
    FileExporter(string outputfile) {
        m_fout.open(outputfile.c_str(), ios::out);

        string outstring = to_string("", "\t", "id", "scannum", "precMz", "charge",
                                     "preNeuMass", "sigArchive", "inGT", "sigGT", "inArchSpace",
                                     "rankingArch", "Pep[PTM]", "dist", "evalue", "neuMassGT", "pepProb", "iProb");
        m_fout << outstring << endl;
    }

    void
    processing(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo, float dist, float evalue,
               float m_threshold_pvalue, PSMInfo &searchpsminfo, string &prevstring) override {
        cout << "file exporter processing" << endl;
        cout << "modified peptide " << gtinfo->getModifiedPeptide(true) << endl;
        cout << "gtinfo " << gtinfo << "\t  is nullptr: " << (gtinfo == nullptr) << endl;
        searchpsminfo.print();
        searchpsminfo.printSearchHit(0);
        cout << "Start to string " << prevstring << endl;
        const SearchHit &sh0=*searchpsminfo.searchhits.at(0);
        prevstring = to_string(prevstring + "\t", "\t", ranking, gtinfo->getModifiedPeptide(false),
                               dist, evalue, searchpsminfo.precursorNeutMass,
                               sh0.m_modified_peptide,
                               sh0.m_peptideprophet_prob,
                               sh0.m_iprophet_prob);
        cout << "start to write to file " << prevstring << endl;
        m_fout << prevstring << endl;
        cout << "done file processing " << __FILE__ << __LINE__ << endl;
    }

    void processingZeroRanking(int ranking, SearchHit &tophit, MzSpecInfo &scaninfo, SPsmAnnotation *gtinfo, float dist,
                               float evalue, float m_threshold_pvalue, PSMInfo &searchpsminfo,
                               string &prevstring) override {
//        const SearchHit &sh0=searchpsminfo.searchhits.at(0);
//        prevstring = to_string(prevstring + "\t", "\t", ranking, gtinfo->getModifiedPeptide(false), dist,
//                               evalue, searchpsminfo.precursorNeutMass,
//                               sh0.m_modified_peptide,
//                               sh0.m_peptideprophet_prob,
//                               sh0.m_iprophet_prob);
//        m_fout << prevstring << endl;
    }

    void outputstr(string &str) override {
        m_fout << str << endl;
    }

};

// todo: for each spectrum in mzfile, compare the identity in search file
void CMzFileSearchResults::validation(string searchfile, CMzFileReader &mzfile) {
    string outputfile = mzfile.getMzFilename() + "qc.tsv";
    spdlog::get("A")->info("Start validation");
    CHistogram histOfRankings(m_topHitsNum + 1);
    PeptideProphetParser groundtruth(searchfile);
    SVennStats venn;
    SPeptideSet pepSet;
    ICExporter *filexpt = new FileExporter(outputfile);
    ICExporter *expt = new ColorExporter(&venn, &pepSet);

    for (int i = 0; i < res.size(); i++) {
        string outstring = "";
        cout << "\n------BEGIN ----Archive Query result ------BEGIN----" << i << endl;
        SArchiveMatch &arcMatch = res[i];

        long queryidx = arcMatch.m_annNodes->getQueryIdx();
        string filename;
        MzSpecInfo &scaninfo = mzfile.getSpecMetaInfo(queryidx, filename);

        PSMInfo searchpsminfo;
        vector<PSMInfo> psms;
        bool found = groundtruth.getPSMInfobyScanFileName(filename, scaninfo.scan, searchpsminfo);
        bool foundAll = groundtruth.getAllPSMsWithScanFileName(filename, scaninfo.scan, psms);
        bool significantGT = false;
        bool inSpace = false;
        if (found) {
            significantGT = groundtruth.isPSMSignificant(searchpsminfo) and not searchpsminfo.isDecoy(true, 0);
            if (psms.size() > 1) {
                cout << "total number found: " << psms.size() << endl;
                for (auto &eachpsminfo: psms) {
                    cout << "Multiple: " << groundtruth.isPSMSignificant(eachpsminfo) << "\t" << eachpsminfo.isDecoy(
                            true, 0)
                         << endl;
                    eachpsminfo.print();
                }
            }
            inSpace = m_AnnotationDB.inSearchSpace(searchpsminfo, mzfile.getListFile());
        }
        bool sigArchive = isSignificantID(arcMatch, scaninfo);
        venn.updateVenn(found, significantGT, sigArchive, inSpace);

        outstring = to_string(outstring, "\t", i, scaninfo.scan, std::fixed, setprecision(3),
                              scaninfo.precursormz, scaninfo.charge, scaninfo.getPrecursorNeutralMass(), sigArchive,
                              found, significantGT, inSpace);
        if (not found or not significantGT or not inSpace) {
            filexpt->outputstr(outstring);
//            fout << outstring << endl;
            continue;
        }

        float dist, evalue;
        SPsmAnnotation *gtinfo = nullptr;
        cout << "Start validation" << endl;
        int ranking = arcMatch.validation(searchpsminfo, dist, evalue, gtinfo);

        outstring = to_string(outstring + "\t", "\t", ranking, gtinfo->getModifiedPeptide(false), dist,
                              evalue, searchpsminfo.precursorNeutMass,
                              searchpsminfo.searchhits[0]->m_modified_peptide,
                              searchpsminfo.searchhits[0]->m_peptideprophet_prob,
                              searchpsminfo.searchhits[0]->m_iprophet_prob);
        filexpt->outputstr(outstring);

        cout << "Ranking(Zero is NOT CONSISTENT Result found by Archive) " << ranking << "\t" << scaninfo << "\t"
             << scaninfo.getPrecursorNeutralMass()
             << "\t" << gtinfo->getModifiedPeptide(false) << endl;
        if (ranking == 0) {
            cout
                    << "### Error: check this case because we did not find correct anwser even with six index and 256 buckets... "
                    << endl;
            cout << "Scan info: " << scaninfo << endl << "gt info " << endl;
            searchpsminfo.searchhits[0]->print();
            if (sigArchive) {
                venn.mix_InDomain_DiffSig++;
            }
        }
        SearchHit &tophit = *searchpsminfo.searchhits[0];
        histOfRankings.add_data(ranking);

        string prev = "";
        if (ranking > 0) {
            // part L : low rank ones but find consistent
            expt->processing(ranking, tophit, scaninfo,
                             gtinfo, dist, evalue, m_pvalue_th,
                             searchpsminfo, prev);
        } else {
            expt->processingZeroRanking(ranking, tophit, scaninfo, gtinfo, dist,
                                        evalue, m_pvalue_th, searchpsminfo, prev);
        }
    }
    venn.printSummary();
    histOfRankings.display();
    pepSet.printSummary();
}

bool CMzFileSearchResults::isSignificantID(const SArchiveMatch &arcMatch, const MzSpecInfo &scaninfo) const {
    bool sigArchive = false;
    int sigHitsIdx = 0;
    CANSIConsole a;
    for (int j = 0; j < arcMatch.size(); j++) {
        cout << "Checking result by archive: " << j + 1 << "/" << arcMatch.size() << endl;
        const SPsmAnnotation &annot = arcMatch.getGt(j);
        if (arcMatch.isSsmSignificant(j, m_pvalue_th) and arcMatch.isGtSignificant(j)) {
            sigArchive = true;
            sigHitsIdx++;
            string tmp = to_string("Match SPEC_WITH_GOOD(sig)_Annotation: ", "\t", "Hits Ranking", sigHitsIdx,
                                   "EVALUE(GOOD)",
                                   arcMatch.getPvalue(j), "DIST",
                                   arcMatch.getDist(j), "DP",
                                   arcMatch.getDP(j),annot.getModifiedPeptide(false),
                                   annot.iProb, annot.pProb, annot.idx, annot.protein);

            cout << defaultfloat << scaninfo << "\t" << a.getColorStr(tmp, CANSIConsole::YELLOW) << endl;
        } else {
            string tmp = to_string("Matched BAD(insig) Annotation: ", "\t", "General Ranking", j + 1, "EVALUE",
                                   arcMatch.getPvalue(j), "DIST",
                                   arcMatch.getDist(j), "DP",
                                   arcMatch.getDP(j),annot.getModifiedPeptide(false),
                                   annot.iProb, annot.pProb, annot.idx, annot.protein);

            cout << defaultfloat << scaninfo << "\t" << a.getColorStr(tmp, CANSIConsole::RED) << endl;
        }
    }
    return sigArchive;
}

void CMzFileSearchResults::toTsv(const string &outfile, bool withheader, bool wrapLines) {
    SimpleTimer st("export search result to file");
    ofstream fout(outfile, ios::out);
    if (withheader) {
        fout << to_string("", "\t", "QueryId", "HitRank", "HitId", "HitL2dist", "p-value", "STD(p-value)", "Score",
                          "GtIdx", "GtPeptide[PTM]", "GtProtein", "GtiProb", "GtpProb", "GtIsSignificant",
                          "parentMass", "charge", "GtPeptide[PTM]") << endl;
    }
    for (int i = 0; i < res.size(); i++) {
        const SArchiveMatch &ssm = res.at(i);
        int hits_num = ssm.size();
        CAnnSpectra &qr = *(ssm.m_annNodes);
        if (hits_num == 0) continue;

        string line = to_string("", "\t", qr.getQueryIdx());
        for (int j = 0; j < hits_num; j++) {
            SAnnSpectrum &nbs = qr.m_anns[j];
            const SPsmAnnotation &gtinfo = ssm.getGt(j);
            if (not ssm.isGtSignificant(j)) continue;

            line = to_string(line, "\t", to_string("Hit_", "", j + 1), nbs.idx, nbs.dist, nbs.pvaluePartial,
                             nbs.getStdOfPvalue(), gtinfo.score, gtinfo.idx, gtinfo.getModifiedPeptide(false),
                             gtinfo.protein != "" ? gtinfo.protein : "N/A", gtinfo.iProb, gtinfo.pProb,
                             gtinfo.significance,
                             gtinfo.precursorNeutralMass(), gtinfo.charge, gtinfo.getModifiedPeptide(false));
            if (wrapLines) {
                fout << line << endl;
                line = to_string("", "\t", qr.getQueryIdx());
            }
        }
        if (not wrapLines) {
            fout << line << endl;
        }
    }
}

void SArchiveMatch::retrieveAnnotation(CAnnotationDB &annoDB) {
    if (m_annNodes->empty()) { return; }
    SPsmAnnotation gtinfo;

    string line = to_string("", "\t", m_annNodes->getQueryIdx());
    for (int j = 0; j < m_annNodes->size(); j++) {
        SAnnSpectrum &nbs = m_annNodes->m_anns[j];
        annoDB.retrieveGtinfo(nbs.idx, gtinfo);
        m_gtInDB.push_back(gtinfo);
    }
}

void CMzFileSearchResults::retrieveAnnotation(vector<CAnnSpectra *> &annOfQueries) {
    spdlog::get("A")->info("Retrieve annotation for each neighbor...");
    res.assign(annOfQueries.size(), SArchiveMatch());
    Progress ps(res.size(), "searching annotation in database...");
    transform(res.begin(), res.end(), annOfQueries.begin(), res.begin(), [&](SArchiveMatch &x, CAnnSpectra *y) {
        ps.increase();
        x.m_annNodes = y;
        x.retrieveAnnotation(m_AnnotationDB);
        return x;
    });

}

void CMzFileSearchResults::calcFDRThresholdPvalue(const string &outfile, bool fdr_use_tophit_pvalues_only) {
    vector<double> p_values;
    createPvalueList(p_values, fdr_use_tophit_pvalues_only);
    CVisual::gnuplot_histogram(p_values, 5, 0.0, 1.0);
    CVisual::gnuplotWrapper info;
    CVisual::SImageSize imsize;
    info.set_filename(outfile + "pvalue.png")
            .set_xlabel("pvalue")
            .set_ylabel("frequency")
            .set_minmax(0, 1)
            .set_height(imsize.m_hpx)
            .set_width(imsize.m_wpx);

    CVisual::gnuplot_histogram_topng(p_values, 20, 0, 1, info);
    File::saveas(p_values, outfile + "_pvalue.txt", true);

    bool use_qvalue = true;
    m_pvalue_th = statistic::BHfdrThreshold(p_values, 0.01, use_qvalue);
    spdlog::get("A")->info("pvalue threshold for FDR < 1% is : {}", m_pvalue_th);
}

void CMzFileSearchResults::createPvalueList(vector<double> &p_values, bool fdr_use_tophit_pvalues_only) {
    SimpleTimer st("export search result to file");
    for (int i = 0; i < res.size(); i++) {
        const SArchiveMatch &arcMatch = res.at(i);
        int hits_num = arcMatch.size();
        if (hits_num == 0) continue;

        for (int j = 0; j < hits_num; j++) {
            if (arcMatch.isGtSignificant(j)) {
                p_values.push_back(arcMatch.getPvalue(j));
                if (fdr_use_tophit_pvalues_only) break;
            }
        }
    }
}

void CMzFileSearchResults::exportTsv(string outfilename) {
    bool outputWithHeader = true;
    bool wrapLines = false;
    bool fdr_use_tophit_pvalues_only = false;
    calcFDRThresholdPvalue(outfilename, fdr_use_tophit_pvalues_only);
    toTsv(outfilename, outputWithHeader, wrapLines);
}

CMzFileSearchResults::CMzFileSearchResults(CAnnotationDB &annodb, int topHitsNum) : m_AnnotationDB(annodb) {
    m_topHitsNum = topHitsNum;
    m_pvalue_th = 0;
}


double intTol2double(int tol) {
    return tol * 2000.0 / 65535 * 2;
}

// ----------------
// to be used in the future
//
class CArchiveMatchValidator {
    float dist;
    float evalue;
    SPsmAnnotation *matched_gtinfo_ptr;
    string gtPeptide;
    int ranking;
    vector<string> promisingAlternativeSolutionPeptides;
    bool found;
    const PSMInfo &psminfo;
public:
    CArchiveMatchValidator(const PSMInfo &searchpsminfo) : psminfo(searchpsminfo) {
        gtPeptide = searchpsminfo.searchhits.at(0)->m_peptide;
        ranking = 0;
    }

    void getResult(float &d, float &e, SPsmAnnotation *&g, int &r) {
        d = dist;
        e = evalue;
        g = matched_gtinfo_ptr;
        r = ranking;
    }

    void validation(SArchiveMatch *arcMatch) {
        if (arcMatch->m_annNodes->empty()) {
            cout << "... no neighbor  0" << endl;
            //ranking = 0;
        } else {
            cout << "Validation " << endl;
            found = arcMatch->validationAllNeighbors(dist, evalue, matched_gtinfo_ptr, gtPeptide,
                                                     promisingAlternativeSolutionPeptides,
                                                     ranking);
            if (not found) {
                cout << "Correct peptide by another search engine not found: " << gtPeptide << "\tTwo probabilities\t"
                     << psminfo.searchhits[0]->m_peptideprophet_prob << "\t"
                     << psminfo.searchhits[0]->m_iprophet_prob << endl;
                cout << "Start" << endl;
                for (int i = 0; i < promisingAlternativeSolutionPeptides.size(); i++) {
                    cout << "i = " << i << "\t" << promisingAlternativeSolutionPeptides[i] << endl;
                }
                cout << "end" << endl;
                ranking = 0; // Zero means not OK
            }
        }

        cout << "End validation " << __FILE__ << __LINE__ << __FUNCTION__ << endl;
    }


};

//------------------
// to be replace
int SArchiveMatch::validation(const PSMInfo &searchpsminfo, float &dist, float &evalue,
                              SPsmAnnotation *&matched_gtinfo_ptr) {
    int rankingInSignificantlyAnnotatedNeighbors = 0;
    const SearchHit &sh0 = *searchpsminfo.searchhits.at(0);
    string gtPeptide = sh0.m_peptide;

    if (not m_annNodes->empty()) {
        bool found = false;
        vector<string> promisingAlternativeSolutionPeptides;
        for (long j = 0; j < m_gtInDB.size(); j++) {
            // for each of our result in archive search
            SPsmAnnotation &gtinfo = m_gtInDB[j];
            if (gtinfo.significance) {
                // if it is significantly labeled node
                rankingInSignificantlyAnnotatedNeighbors++;
                if (comparePeptide_I_equls_L(gtPeptide, gtinfo.peptideseq)) {
                    // if it is significantly labeled node,  and the same as the gt sequence ,we find it, remember it and stop
                    found = true;
                    dist = m_annNodes->m_anns[j].dist;
                    evalue = m_annNodes->m_anns[j].pvaluePartial;
                    matched_gtinfo_ptr = &gtinfo;
                    break;
                } else {
                    // if it is significant, but not the same as the gt sequence, we may find something like a mixture spectra,print it
                    // print it to screen.
                    ostringstream oss;
                    gtinfo.toOstringStreamNoId(oss);
                    string tmp = to_string("Significantly by Archive But DIFFERENT from another Search engine", "\t",
                                           "\nGROUNDTRUTH", gtPeptide, "OurResult", oss.str(),
                                           gtinfo.getModifiedPeptide(false),
                                           m_annNodes->m_anns[j].dist,
                                           "pvalue:",
                                           m_annNodes->m_anns[j].pvaluePartial,
                                           "mass",
                                           gtinfo.precursorNeutralMass()
                    );
                    promisingAlternativeSolutionPeptides.push_back(tmp);
                    // mixture...
                    if (matched_gtinfo_ptr == nullptr) {
                        dist = m_annNodes->m_anns[j].dist;
                        evalue = m_annNodes->m_anns[j].pvaluePartial;
                        matched_gtinfo_ptr = &gtinfo;
                    }
                }
            } else {
                // the result in archive is not significant
                cout << "GTinfo annotation insignificant" << endl;
                if (comparePeptide_I_equls_L(gtPeptide, gtinfo.peptideseq)) {
                    // it is not significantly labeled, but the sequence is the same as gt sequence,
                    // we should make better annotation, maybe this is significant by another search engine, other than msfragger..
                    cout << "WOW------------inSignificant annotation was correct!! improve the annotation!!" << endl;
                    ostringstream oss;
                    gtinfo.toOstringStreamNoId(oss);
                    string tmp = to_string(
                            "Insignificantly annotated Candidates in spectral archive But CONSIST TO another Search engine",
                            "\t", "\nGROUNDTRUTH", gtPeptide, "OurResult", oss.str(),
                            gtinfo.getModifiedPeptide(false),
                            m_annNodes->m_anns[j].dist,
                            "pvalue:",
                            m_annNodes->m_anns[j].pvaluePartial,
                            "mass",
                            gtinfo.precursorNeutralMass()
                    );
                    promisingAlternativeSolutionPeptides.push_back(tmp);
                }
            }
        }
        // end of for loop
        if (not found) {
            // the significantly labeled ones are all wrong...
            // print them.. the ranking is 0
            cout << "Correct peptide by another search engine not found: " << gtPeptide << "\tTwo probabilities\t"
                 << sh0.m_peptideprophet_prob << "\t"
                 << sh0.m_iprophet_prob << endl;
            for_each(promisingAlternativeSolutionPeptides.begin(), promisingAlternativeSolutionPeptides.end(),
                     [](const string &x) { cout << x << endl; });
            rankingInSignificantlyAnnotatedNeighbors = 0; // Zero means not OK
        }
    }
    return rankingInSignificantlyAnnotatedNeighbors;
}

double SArchiveMatch::getPvalue(int i) const {
    return m_annNodes->m_anns.at(i).pvaluePartial;
}

double SArchiveMatch::getDist(int i) const {
    return m_annNodes->m_anns.at(i).dist;
}

double SArchiveMatch::getDP(int i) const {
    const int MAX_TOP50_COS = 42925;
    return m_annNodes->m_anns.at(i).dotprod * 1.0 / MAX_TOP50_COS;
}

bool SArchiveMatch::isSsmSignificant(int i, double pvaluethreshold) const {
    return getPvalue(i) < pvaluethreshold;
}

bool SArchiveMatch::isGtSignificant(int i) const {
    return getGt(i).significance == 1;
}

const SPsmAnnotation &SArchiveMatch::getGt(int i) const {
    return m_gtInDB.at(i);
}

int SArchiveMatch::size() const {
    if (m_annNodes->size() != m_gtInDB.size()) {
        throw logic_error("size not consistent");
    }
    return m_gtInDB.size();
}

long SArchiveMatch::getQueryIdx() const {
    return m_annNodes->getQueryIdx();
}


bool
SArchiveMatch::validationAllNeighbors(float &dist, float &evalue, SPsmAnnotation *&matched_gtinfo_ptr, string gtPeptide,
                                      vector<string> &promisingAlternativeSolutionPeptides, int &ranking) {
    bool found = false;
    for (long j = 0; j < m_gtInDB.size(); j++) {
        SPsmAnnotation &gtinfo = m_gtInDB[j];
        cout << "Processing j=" << j << endl;
        if (isGtSignificant(j)) {
            cout << "IF " << endl;
            ranking++;
            if (comparePeptide_I_equls_L(gtPeptide, gtinfo.peptideseq)) {
                cout << "EQUAL " << endl;
                // if it is significantly labeled node,  and the same as the gt sequence ,we find it, remember it and stop
                found = true;
                dist = m_annNodes->m_anns[j].dist;
                evalue = m_annNodes->m_anns[j].pvaluePartial;
                matched_gtinfo_ptr = &gtinfo;
                break;
            } else {
                cout << "NOT EQUAL " << endl;
                // if it is significant, but not the same as the gt sequence, we may find something like a mixture spectra,print it
                // print it to screen.
                ostringstream oss;
                gtinfo.toOstringStreamNoId(oss);
                string tmp = to_string("Significantly by Archive But DIFFERENT from another Search engine", "\t",
                                       "\nGROUNDTRUTH", gtPeptide, "OurResult", oss.str(),
                                       gtinfo.getModifiedPeptide(false),
                                       m_annNodes->m_anns[j].dist,
                                       "pvalue:",
                                       m_annNodes->m_anns[j].pvaluePartial,
                                       "mass",
                                       gtinfo.precursorNeutralMass()
                );
                promisingAlternativeSolutionPeptides.push_back(tmp);

                if (matched_gtinfo_ptr == nullptr) {
                    dist = m_annNodes->m_anns[j].dist;
                    evalue = m_annNodes->m_anns[j].pvaluePartial;
                    matched_gtinfo_ptr = &gtinfo;
                }
            }
        } else {
            cout << "ELSE" << endl;
            // the result in archive is not significant
            cout << "GTinfo annotation insignificant " << __FILE__ << __FUNCTION__ << __LINE__ << endl;
            if (comparePeptide_I_equls_L(gtPeptide, gtinfo.peptideseq)) {
                // it is not significantly labeled, but the sequence is the same as gt sequence,
                // we should make better annotation, maybe this is significant by another search engine, other than msfragger..
                cout << "WOW------------inSignificant annotation was correct!! improve the annotation!!" << endl;
                ostringstream oss;
                gtinfo.toOstringStreamNoId(oss);
                string tmp = to_string(
                        "Insignificantly annotated Candidates in spectral archive But CONSIST TO another Search engine",
                        "\t", "\nGROUNDTRUTH", gtPeptide, "OurResult", oss.str(),
                        gtinfo.getModifiedPeptide(false),
                        m_annNodes->m_anns[j].dist,
                        "pvalue:",
                        m_annNodes->m_anns[j].pvaluePartial,
                        "mass",
                        gtinfo.precursorNeutralMass()
                );
                promisingAlternativeSolutionPeptides.push_back(tmp);
            }
        }
    }
    cout << "END of for loop" << endl;
    return found;
}




SAnnGTSummary::SAnnGTSummary() {
    m_recallOfTrueNeighbor = false;
    totalnum = 0;
    correctnum = 0;

    m_nprobe = 0;
    m_indexNum = 0;
    m_recallTNNminDP = 0;
    m_minPeakNumInSpec = 1;
    logFileName = "TNNRecall.log";

    ppp = nullptr;
    
    
    
}

void SAnnGTSummary::print() {
    double ratio = correctnum * 1.0 / totalnum;
    string tmp = to_string("", " ", "correct", correctnum, "totalnum", totalnum, "recall", ratio, "\n");
    writeToLogFile(tmp);
    spdlog::get("A")->info("TNNRecall found {} #TNN_total {} overall_recall {:.4f} nprobe {} indexNum {} minDP {} minPeak {} topK {}", correctnum, totalnum, correctnum * 1.0 / totalnum, m_nprobe, m_indexNum, m_recallTNNminDP, m_minPeakNumInSpec, m_recallTNNtopK);
    // cout << "correct " << correctnum << " totalnum " << totalnum << " ratio " << correctnum * 1.0 / totalnum << endl;
}

// todo: to be replaced
// The validation on pepxml file is never used?
void SAnnGTSummary::increase(bool iscorrect, long queryindex) {
    if (queryindex != -1 and ppp != nullptr and m_queryindex2scan.size() > 0 and m_queryindex2filename.size() > 0) {
        int scan = m_queryindex2scan[queryindex];
        string filename = m_queryindex2filename[queryindex];

        PSMInfo psminfo;
        ppp->getPSMInfobyScanFileName(filename, scan, psminfo);

        bool isSignificant = ppp->isPSMSignificant(psminfo);
        if (not isSignificant) {
            cout << "not significant scan: skip it" << endl;
            return;
        }
    }
    totalnum++;
    if (iscorrect) correctnum++;
}

void SAnnGTSummary::setvalidationfile(string peptideprophetpepxmlfile) {
    if (ppp == nullptr and File::isExist(peptideprophetpepxmlfile,true)) {
        ppp = make_shared<PeptideProphetParser>(peptideprophetpepxmlfile);
    } else {
        //cout << "ppp pointer can be only set once!" << endl;
    }
}

// not used, to be deleted!
void SAnnGTSummary::setqueryindexMap(CMzFileReader &querySpectra) {
    string filename;
    for (long i = 0; i < querySpectra.getSpecNum(); i++) {
        MzSpecInfo &meta = querySpectra.getSpecMetaInfo(i, filename);
        m_queryindex2scan[i] = meta.scan;
        m_queryindex2filename[i] = filename;
    }
}

void CSocketServerSummary::print() {
    cout << "----------------------------------------------------" << endl
         << "#XHR\tTime\t#Total(s)" << endl << setprecision(2)
         << num_of_search_done << "\t" << time_used_for_current_search << "\t"
         << total_time_used_for_search << endl
         << "------------------------------------------------------" << endl;
}
std::mutex summary_lock;
void CSocketServerSummary::update(double time_used_in_sec) {
    std::lock_guard<mutex> guard(summary_lock);
    time_used_for_current_search = time_used_in_sec;
    total_time_used_for_search += time_used_for_current_search;
}

CSocketServerSummary::CSocketServerSummary() {
    num_of_search_done = 0;
    time_used_for_current_search = 0;
    total_time_used_for_search = 0;
}
