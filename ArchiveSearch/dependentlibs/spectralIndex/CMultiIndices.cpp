//
// Created by wulong on 1/25/19.
//
#include <vector>
#include "../../dependentlibs/msbasic/commonfun.h"
#include "../../../librarymsms/CMzFileReader.h"
#include "../../../librarymsms/ConcreteQueries.h"
#include "../../../librarymsms/Util.h"
#include "../../../librarymsms/ProteomicsDataTypes.h"
#include "../../../librarymsms/XMLFileParser.h"
#include "CMultiIndices.h"
#include "CKMeans.h"
#include "ICIndexWrapper.h"

MultipleIndicesImpl::MultipleIndicesImpl() : m_useCpu(true) {
}

// create multiple indices, set the filename
// assign seeds from string "indexshufleseeds" to N indices, to make random algorithm reproducible.
// the indexshuffleseeds can be "default", which means using 0, to N-1 as seeds.
// the indexshuffleseeds can be "customized:1;2;3;4;5;6" which means assign 1 2 3 4 5 6 as seeds for N indices.
// Attention: in different version of gcc, the shuffle function changes. 
//            This might cause problem. So use the proper shuffle function.
// This function will create no indexes, but make the parameters ready.
void MultipleIndicesImpl::initialize(int numIndex, string indexpath, string basename, bool useMyOwn,
                                     shared_ptr<CPQParam> option, string indexshuffleseeds) {
    m_numIndexUsed = numIndex;
    m_ShuffleIndex.assign(numIndex, nullptr);

    for (int i = 0; i < m_ShuffleIndex.size(); i++) {
        m_ShuffleIndex[i] = IndexFactory(useMyOwn, option);
    }

    string filename_suffix = ".index";
    if (useMyOwn) filename_suffix = "_MY.index";

    for (int i = 0; i < numIndex; i++) {
        string filename = indexpath + "/" + basename + to_string(i + 1) + filename_suffix;
        getShuffledIndex(i)->setfilename(filename);
    }

    // seeds should be initialized in a better way.
    m_seeds.assign(numIndex, 0);
    if(indexshuffleseeds == "default"){
        iota(m_seeds.begin(), m_seeds.end(), 0);
    }else{
        int pos = indexshuffleseeds.find_first_of(":")    ;
        if(string::npos == pos) cout << "Error: invalid seeds" << endl;
        indexshuffleseeds = indexshuffleseeds.substr(pos+1);
        vector<string> values;
        split_string(indexshuffleseeds, values, ';');
        if(values.size() != m_seeds.size()) cout << "Error: invalid seeds" << endl;
        for (int i = 0; i < values.size(); i ++){
            m_seeds[i] = stringTo<int>(values[i]);
        }


    }
    cout << "Seeds\t" ;
    for(auto x: m_seeds) cout << x << "\t";
    cout << endl;

}

void MultipleIndicesImpl::trainIndex(int i, int batchSize, float *vecBatch) {
    // get a random data set from vec todo:
    spdlog::get("A")->info("Start training {}-th index...", i);
    getShuffledIndex(i)->train(batchSize, vecBatch);

    spdlog::get("A")->info("Finish training {}-th index!", i);
}
std::mutex search_lock_on_gpu;
void MultipleIndicesImpl::search(int i, int ret_num, vector<float> &dist, vector<long> &ind, float *shuffledVquery,
                                 int numQuery) {
    if (m_useCpu) {
//        cout << "[Info] Using CPU: " << endl;
        getShuffledIndex(i)->search(numQuery, shuffledVquery, ret_num, dist, ind);
    } else {
//        cout << "[Info] Using GPU thread safe: " << endl;
        std::lock_guard<std::mutex> guard(search_lock_on_gpu);
        getShuffledIndex(i)->search(numQuery, shuffledVquery, ret_num, dist, ind);
    }

}

void MultipleIndicesImpl::createEmptyIndex(int i, string indexstr, int dim) {
    getShuffledIndex(i)->createEmptyIndex(dim, indexstr);
}

void MultipleIndicesImpl::train(long specnum, int dim, float *vec) {
    cout << "Total number of index to create: " << getNum() << endl;
    int batchSize = specnum / getNum();
    spdlog::get("A")->info("batchSize = {} numIndex = {}", batchSize, getNum());
    float *vecBatch = new float[batchSize * dim];
    CShuffledVector shuffleTool(specnum);
    for (int i = 0; i < getNum(); i++) {
        int offset = i * batchSize;
        shuffleTool.getVector(vec, offset, dim, batchSize, vecBatch);
        shuffleVectors(i, batchSize, vecBatch);
        trainIndex(i, batchSize, vecBatch);
    }
    delete[] vecBatch;
}

void MultipleIndicesImpl::add(int i, int batchsize, float *shuffleVec) {
    getShuffledIndex(i)->add(batchsize, shuffleVec);
}

void MultipleIndicesImpl::addshuffle(int i, float *vec, long specnum) {

    int dim = getShuffledIndex(i)->dim();
    int smallbatch = 10000;
    float *shufflevec = new float[smallbatch * dim];
    // so it is safe to have zero spectrum,spectrum = 0, skip
    for (int k = 0; k < specnum; k += smallbatch) {
        int start = k, end = (k + smallbatch) > specnum ? specnum : k + smallbatch;

        copy(vec + start * dim, vec + end * dim, shufflevec);

        int numvec = end - start;
        shuffleVectors(i, numvec, shufflevec);
        add(i, numvec, shufflevec);
    }
    delete[] shufflevec;
}

void MultipleIndicesImpl::add(float *vec, long specnum) {
    cout << "append spectrum to " << getNum() << " indices: " ;
    for (int i = 0; i < getNum(); i++) {
        cout << i +1 << " " << flush;
//        cout << "append to " << i << "-th index, #spectra:  " << specnum << endl;
        addshuffle(i, vec, specnum);
    }
    cout << endl;
}

void MultipleIndicesImpl::shufflevector(int seed, float *p, int dim) {
    gcc5shuffle(p, p + dim, mt19937(seed));
}

void MultipleIndicesImpl::shuffleVectors(int i, int numVec, float *p) {
    int dim = getDim(i);

    for (int offset = 0; offset < numVec * dim; offset += dim) {
        shufflevector(m_seeds[i], p + offset, dim);
    }
}

void MultipleIndicesImpl::createEmptyIndices(vector<string> &indexstrs, int dim) {
    //cout << "getNUM " << getNum() << endl;
    for (int i = 0; i < getNum(); i++) {
        createEmptyIndex(i, indexstrs[i], dim);
    }
}

void MultipleIndicesImpl::toGPU() {
    if (m_useCpu) {
        for (int i = 0; i < getNum(); i++) {
            toGPU(i);
        }
        m_useCpu = false;
    }
}

void MultipleIndicesImpl::toGPU(int i) {
    if (m_useCpu) {
        getShuffledIndex(i)->toGPU();
    }
}

void MultipleIndicesImpl::toCPU() {
    if (not m_useCpu) {
        for (int i = 0; i < getNum(); i++) {
            getShuffledIndex(i)->toCPU();
        }
        m_useCpu = true;
    }
}

void MultipleIndicesImpl::display(int i) {
    getShuffledIndex(i)->display();
}

void MultipleIndicesImpl::display() {
    // cout << "--------Multi-Indices Info-----begin----" << endl;
    for (int i = 0; i < getNum(); i++) {
        display(i);
        spdlog::get("A")->info("index: d = {}, is_trained = {}, ntotal ={}", getDim(i), istrained(i), ntotal(i));
    }
    // cout << "--------Multi-Indices Info----end-----" << endl;
}

bool MultipleIndicesImpl::istrained(int i) {
    return getShuffledIndex(i)->istrained();
}

void MultipleIndicesImpl::read() {
    for (int i = 0; i < getNum(); i++) {
        cout << "\t" << i ;
        readIndex(i);
    }
    cout << endl;
}

void MultipleIndicesImpl::readIndex(int i) {
    getShuffledIndex(i)->ICIndexWrapper::read();
//    spdlog::get("A")->info("Index loaded from file: {}", getIndexFileName(i));
}

void MultipleIndicesImpl::readIndex(int i, string filename) {
    getShuffledIndex(i)->read(filename);
}

void MultipleIndicesImpl::write() {
    cout << "Attention: writing to index in use may result in inconsistence between the indice in use and not in use" << endl;
    for (int i = 0; i < getNum(); i++) {
        writeIndex(i);
    }
}

void MultipleIndicesImpl::writeIndex(int i, string filename) {
    getShuffledIndex(i)->write(filename);
}

void MultipleIndicesImpl::writeIndex(int i) {
    getShuffledIndex(i)->ICIndexWrapper::write();
    spdlog::get("A")->info("Index is written to file: {}", getIndexFileName(i));
}

bool MultipleIndicesImpl::ifAllFileExist() {
    bool pass = true;
    for (int i = 0; i < getNum(); i++) {
        if (not ifFileExist(i)) pass = false;
    }
    return pass;
}

bool MultipleIndicesImpl::ifFileExist(int i) {
    return File::isExist(getIndexFileName(i),true);
}

long MultipleIndicesImpl::ntotal(int i) {
    return getShuffledIndex(i)->total();
}

std::string MultipleIndicesImpl::getIndexFileName(int i) {
    return getShuffledIndex(i)->getfilename();;
}

int MultipleIndicesImpl::getDim(int i) {
    return getShuffledIndex(i)->dim();
}

ICIndexWrapper *MultipleIndicesImpl::getShuffledIndex(int i) {
    return m_ShuffleIndex[i].get();
}

int MultipleIndicesImpl::getNum() const {
    return m_numIndexUsed;
}

void MultipleIndicesImpl::setNprobe(int nprobe) {
    m_nprobe = nprobe;
    for (int i = 0; i < getNum(); i++) {
        getShuffledIndex(i)->setnProb(nprobe);
    }
}

void MultipleIndicesImpl::getANNByShuffleQueries(int numQuery, int ret_num, const float *vquery,
                                                 vector<vector<long>> &multipleInd,vector<vector<double>> &results_dist, int i, vector<float> &dist,
                                                 vector<long> &ind) {
    int dim = getDim(i);
    long totalLen = numQuery * dim;
    float *shuffledVquery = new float[numQuery * dim];
    copy(vquery, vquery + totalLen, shuffledVquery);

    shuffleVectors(i, numQuery, shuffledVquery);
    search(i, ret_num, dist, ind, shuffledVquery, numQuery);

    multipleInd.push_back(vector<long>(ind.begin(), ind.end()));
    results_dist.push_back(vector<double>(dist.begin(), dist.end()));

    delete[] shuffledVquery;
}

void MultipleIndicesImpl::removeIds(vector<long> &idx) {
    for(int i = 0; i < getNumOfIndexBuilt(); i ++){
        // for each index, remove those ids.
        m_ShuffleIndex[i]->removeIds(idx);
    }
}

void MultipleIndicesImpl::setNum(int indexNum) {
    if (indexNum > 0 and indexNum <= getNumOfIndexBuilt()) m_numIndexUsed = indexNum;
//    cout << "[Info] only using " << m_numIndexUsed << " indexing files " << endl;
}

bool MultipleIndicesImpl::usingCPU() {
    return m_useCpu;
}

CMultiIndices::CMultiIndices(string indexstring, string indexshuffleseeds,string indexpath, string libname, bool doShuffle,
                             double tolerance, bool useMyOwn, shared_ptr<CPQParam> cpqParam, const int topPeakNum,
                             bool removePrecursor, bool useFlankingBins, int dim): TRAINING_SPEC_NUM(100000) {

    m_topPeakNum = topPeakNum;
    option = cpqParam;
    m_tolerance = tolerance;
    m_doShuffle = doShuffle;
    split_string(indexstring, m_multiIndicesStr, ';');
//    printIndicesStr();
    m_indexNum = m_multiIndicesStr.size();
    // custom seeds from here.
    m_impl.initialize(m_indexNum, indexpath, libname, useMyOwn, option, indexshuffleseeds);
    m_debug = false;
    m_removeprecursor = removePrecursor;
    m_useFlankingBins = useFlankingBins;
    m_dim = dim;
}

// todo: this function is not widely used!
void CMultiIndices::create( long &specnum, DataFile &splib,
                           bool removeprecursor, bool useFlankingBins) {
    if (m_impl.ifAllFileExist()) {
        cout << "Loading index from disk: " << flush;
        m_impl.read();
    } else {
        trainOnSingle(m_dim, specnum, splib, removeprecursor, useFlankingBins);
        float *vec = splib.toFloatVector(m_dim, specnum, removeprecursor, useFlankingBins, m_tolerance, 0, -1,
                                         m_topPeakNum); // ok
        m_impl.add(vec, specnum);
        delete[] vec;
        m_impl.write();
    }
}

// todo: train and add should be two step.
void CMultiIndices::trainOnFileList(string datafilelist) {
    vector<string> files = readlines(datafilelist);
    long total_specnum =0;
    float *vec =  collectTrainingSpectra(files, total_specnum);
    train(m_dim, total_specnum, vec);

    delete[] vec;


}

float * CMultiIndices::collectTrainingSpectra(const vector<string> &files, long & total) const {

    total= 0;
    // m_dim=4096; training_spec_num=100K; in total: 4096 X 100K = 4K X 100 K = 400 M
    float * vec= new float[TRAINING_SPEC_NUM * m_dim];
    vector<int> file_idx(files.size(), 0);

    iota(file_idx.begin(), file_idx.end(), 0);
    int seed = 42;
    gcc5shuffle(file_idx.begin(), file_idx.end(), mt19937(seed));

    cout << "[Info] Start to collect " << TRAINING_SPEC_NUM << " MS2 spectra for training.. " << endl;
    for (int i = 0; i < file_idx.size(); i++) {
        long count = 0;
        string filename = files[file_idx[i]];

        DataFile df(filename,0,-1);
        float *p = df.toFloatVector(m_dim, count, m_removeprecursor, m_useFlankingBins, m_tolerance, 0, -1,
                                    m_topPeakNum);
        long newSpecNum = total + count <= TRAINING_SPEC_NUM ? count : TRAINING_SPEC_NUM - total;
        copy(p, p + m_dim * newSpecNum, vec + m_dim * total);
        total += newSpecNum;
        spdlog::get("A")->info("Processing file {} / {} : {} : training set size {} / {}, new spectra {}", i, file_idx.size(), filename, total, TRAINING_SPEC_NUM, newSpecNum);

        delete[] p;
        if (total == TRAINING_SPEC_NUM) {
            break;
        }
    }
    return vec;
}

// todo: not very important!
void CMultiIndices::trainOnSingle(int dim, long &specnum, DataFile &splib, bool removeprecursor, bool useFlankingBins) {
    float *vec = splib.toFloatVector(dim, specnum, removeprecursor, useFlankingBins, m_tolerance); // not very important
    train(dim, specnum, vec);
    delete[] vec;
}

// index saved after train
void CMultiIndices::train(int dim, long &specnum, float *vec) {
    //cout << "Start Training" << endl;
    m_impl.createEmptyIndices(m_multiIndicesStr, dim);
    m_impl.train(specnum, dim, vec);
    m_impl.write(); // save index after train
    //cout << "End of training" << endl;
}

void CMultiIndices::toGpu() {
    m_impl.toGPU(); // todo
}

void CMultiIndices::toCpu() {
    m_impl.toCPU();
}


void CMultiIndices::recallOfAnn(DataFile &df, string ipropepxmlfilename, CMzFileReader &compactRawData,
                                bool useflankingbin, int dim, DataFile &splib, bool use_gpu, int indexChoice,
                                int batchsize, string outputpath,
                                int ret_num, long spec_start, long spec_end) {
    // This function is only used once. to be removed! in the future.
    if (not File::isExist(ipropepxmlfilename)) {
        spdlog::get("A")->error("File does not exist: {}", ipropepxmlfilename);
        throw runtime_error("Fail to open file!");
    }
    string outputfile = getOutFileName(ipropepxmlfilename, outputpath);
    ofstream fout(outputfile, ios::out);

    PeptideProphetParser ppp(ipropepxmlfilename);
    double probThreshold = ppp.getThresholdForFDR(0.01, false);

    int k = spec_start >= 0 and spec_start <= compactRawData.getSpecNum() ? spec_start : 0;
    spec_end = k <= spec_end and spec_end <= compactRawData.getSpecNum() ? spec_end : compactRawData.getSpecNum();
    cout << "[Info] Using default spec range [" << spec_start << ", " << spec_end << "]" << endl;

    if (use_gpu) {
        m_impl.toGPU();
    }

    CRecallStat recallANNs;
    Progress ps(spec_end - spec_start, "ANN-Recall Analysis");
    while (k < spec_end) {
        vector<vector<long>> results;
        vector<vector<double>> results_dist;
        int numQuery = batchsize > spec_end - k ? spec_end - k : batchsize;
        vector<long> queries(numQuery, 0);
        iota(queries.begin(), queries.end(), k);
        k += numQuery;
        CMzQuery mzquery(compactRawData, dim, queries, useflankingbin);
        getAnns(mzquery, ret_num, results,results_dist,-1);
        mzquery.print(8);
        verifyEachQuery(df, splib, indexChoice, ret_num, ppp, probThreshold, ps, results, results_dist,
                        queries, m_impl.getNum(), fout, recallANNs);
    }
    fout.close();
    spdlog::get("A")->info("Validation of ANNs has been written into: {}", outputfile);
    spdlog::get("A")->info("ANN recall: {} / {} = {}", recallANNs.m_correct, recallANNs.m_total, recallANNs.recall());
    m_impl.toCPU();
}

void CMultiIndices::verifyEachQuery(DataFile &df, DataFile &splib, int indexChoice, int ret_num,
                                    PeptideProphetParser &ppp, double probThreshold, Progress &ps,
                                    const vector<vector<long>> &results,const vector<vector<double>> &result_dists,
                                    const vector<long> &queries, int indexNum, ofstream &fout,
                                    CRecallStat &recallANNs) const {
    // this function is only used once, to be removed in the future .
    int numQuery = queries.size();
    for (int i = 0; i < numQuery; i++) {
        ps.increase();
        int idx = queries[i];

        CSpectrum *spec = df.getSpectrum(df.getIdx(idx));
        if (spec and spec->getMSLevel() == 2) {
            int scan = spec->getScanNum();
            PSMInfo psminfo;
            bool found_psm = ppp.getPSMInfobyScan(scan, psminfo);

            if (found_psm and ppp.isPSMSignificant(psminfo) and
            not psminfo.isDecoy(true, 0)) {
                recallANNs.verify();
                vector<long> int_ind;
                vector<double> dist_ann;
                collectANNs(indexChoice, ret_num, results, result_dists,i, int_ind, dist_ann, false);
                if (scan == 27504 and m_debug) {
                    cout << "id" << i << "\tQuery: " << idx << endl << "Neighbors: ";
                    for (auto e : int_ind) cout << e << "\t";
                    cout << endl;
                }

                int top_recall = findPeptide(splib, int_ind, psminfo.searchhits[0]->m_modified_peptide);

                if (top_recall >= 0) {
                    recallANNs.correct();
                }

                if (top_recall < 0) { fout << "NO\t"; }
                else { fout << "YES\t"; }
                fout << i << "\t" << scan << "\t" << psminfo.searchhits[0]->m_modified_peptide << "\t" << top_recall
                     << "\t" << psminfo.getProtein_UseAlterProteinIfItsNotDecoy(true, 0) << endl;
            }
        }
    }
}


int CMultiIndices::findPeptide(DataFile &splib, const vector<long> &int_ind, const string &peptide) const {
    /*
     * Looking for a peptide from a spectral library with index list: int_ind
     * This is used to verify whether a list of idx(=scan-1) values are refering to spectra annotated with "peptide"
     * This function should moved out of this CMultiIndices class because it is working alone.
     *
     * todo:
     *  * Split the function into two.
     *  1. retrieve a peptide list from a spectral library based on a list of index or (scan-1)
     *  2. Look up for peptide from the list
     * */
    int top_recall = -1;
    for (int j = 0; j < int_ind.size(); j++) {
        long neighbouridx = int_ind[j];
        CSpectrum *libspec = splib.getSpectrum(neighbouridx);
        if (libspec == nullptr) {
            cout << "Invalid neighbour: " << neighbouridx << endl;
            break;
        }
        string lib_peptide = libspec->getSpectrumName();
        lib_peptide = lib_peptide.substr(0, lib_peptide.find_last_of("/"));

        if (lib_peptide == peptide) {
            top_recall = j + 1;
            break;
        }
    }
    return top_recall;
}

void CMultiIndices::getAnns(ICQuery &q, int ret_num, vector<vector<long>> &results,vector<vector<double>> &results_dist,int indexNum) {
    m_impl.setNum(indexNum);
    float *vquery = q.L2Normalization().get();
    getAnnsForQueries(q.size(), ret_num, vquery, results, results_dist,false);
//    cout << "size: " << q.size() << "\t" << results.size() << endl;
}

void CMultiIndices::getAnnsForQueries(int numQuery, int ret_num, float *vquery, vector<vector<long>> &multipleInd,vector<vector<double>> &results_dist,
                                      bool verbose) {
    if(verbose)cout << "Number of ANNs in each index: " ;
//    cout <<
    SimpleTimer st;
    string tmp = to_string("ANNTime\t","\t","useCPU", m_impl.usingCPU(), "archiveSize",this->total(),"queryNum", numQuery, "indexNum", m_impl.getNum(), "nprobe", m_impl.getNprobe());
    for (int i = 0; i < m_impl.getNum(); i++) {
        double now = st.secondsElapsed();
        vector<float> dist(numQuery * ret_num, 0);
        vector<long> ind(numQuery * ret_num, 0);
        m_impl.getANNByShuffleQueries(numQuery, ret_num, vquery, multipleInd, results_dist, i, dist, ind);
        if(verbose) cout << i << ": " << multipleInd[i].size() <<" ";
        double time_used = st.secondsElapsed()-now;
        tmp = to_string("","\t", tmp, time_used);
    }
    spdlog::get("A")->info(tmp);
    if(verbose)cout << endl;
}

string CMultiIndices::getOutFileName(const string &ipropepxmlfilename, const string &outputpath) const {
    string ipropath, iproname;
    File::parent(ipropepxmlfilename, ipropath, iproname);
    string outputfile = outputpath + "/" + iproname + "_ANN.recall";
    return outputfile;
}

// plan to replace the CollectANN function
class CAssembleIdx {
    int m_indexNum; // 6
    int m_retNumPerIndex; // 1024
    int m_indexChoice; // -1, all
public:
    CAssembleIdx(int indexNum, int retnumperindex, int indexchoice) {
        m_indexChoice = indexchoice;
        m_retNumPerIndex = retnumperindex;
        m_indexNum = indexNum;
    }

    void collectANNs(const vector<vector<long>> &results, int idx,
                     vector<long> &int_ind) const {

        if (m_indexChoice == -1) {
            for (int j = 0; j < m_indexNum; j++) {
                int_ind.insert(int_ind.end(), results[j].begin() + idx * m_retNumPerIndex,
                               results[j].begin() + (idx + 1) * m_retNumPerIndex);
            }
            stableUniqueVector_deprecated(int_ind, false); // not used
        } else {
            int j = m_indexChoice;
            int_ind.insert(int_ind.end(), results[j].begin() + idx * m_retNumPerIndex,
                           results[j].begin() + (idx + 1) * m_retNumPerIndex);
        }
//    cout << "With choice of INDEX " << indexChoice << " #candidates "  << int_ind.size() << endl;
    }

};

// 40% fo time spend here. tried reserve space, but not working. related to the map...
void CMultiIndices::collectANNs(int indexChoice, int ret_num, const vector<vector<long>> &results, const vector<vector<double>> &results_dist, int idx,
                                vector<long> &int_ind, vector<double> &dist, bool verbose) const {
    if (indexChoice == -1) {
        for (int j = 0; j < m_impl.getNum(); j++) {
            int_ind.insert(int_ind.end(), results[j].begin() + idx * ret_num,
                           results[j].begin() + (idx + 1) * ret_num);
            dist.insert(dist.end(), results_dist[j].begin() + idx * ret_num,
                        results_dist[j].begin() + (idx + 1) * ret_num);
        }
//        stableUniqueVector_deprecated(int_ind, true);
        vector<long> tmp(int_ind.begin(), int_ind.end()); // keep a record of the ids.
        stableUniqueVectorSet(int_ind,verbose,ret_num);
        // now keep the corresponding distances.
        map<long, double> idx2dist;

        // this key value pair is used to collect the scores. There should be some better method.
        for(int i = 0; i < tmp.size(); i ++){
            if(idx2dist.count(tmp[i])==0){ // new key. add
                idx2dist[tmp[i]] = dist[i];
            }
        }
        vector<double> tmp_dist;
        tmp_dist.reserve(int_ind.size());

        for(int i = 0; i < int_ind.size(); i++){
            tmp_dist.push_back(idx2dist[int_ind[i]]);
        }
        tmp_dist.swap(dist);

    } else {
        int j = indexChoice;
//        int_ind.reserve(1024);
//        dist.reserve(1024);
        int_ind.insert(int_ind.end(), results[j].begin() + idx * ret_num,
                       results[j].begin() + (idx + 1) * ret_num);
        dist.insert(dist.end(), results[j].begin() + idx * ret_num, results[j].begin() + (idx + 1) * ret_num);
    }
//    cout << "With choice of INDEX " << indexChoice << " #candidates "  << int_ind.size() << endl;
}

long CMultiIndices::total() {
    if (0 == m_impl.getNum()) {
        throw "No index found!";
    }
    long ret = m_impl.ntotal(0);
    for (int i = 1; i < m_impl.getNum(); i++) {
        long x = m_impl.ntotal(i);
        if (x != ret) {
            throw logic_error("multiple indices inconsistent!");
        }

    }
    return ret;
}

void CMultiIndices::write() {
    m_impl.write();
}

void CMultiIndices::setNprobe(int nprobe, bool verbose) {
    if (nprobe >= 1 and nprobe <= 256) {
        if(verbose)cout << "set nprobe of multi-indices as " << nprobe << endl;
        m_impl.setNprobe(nprobe);
    } else {
        if(verbose)cout << "nprobe is out of range!" << endl;
        spdlog::get("A")->error("fail to set nprobe as {}", nprobe);
    }
}

void CMultiIndices::printIndicesStr() {
    for_each(m_multiIndicesStr.begin(), m_multiIndicesStr.end(),
             [](const string &keystr) { cout << keystr << endl; });
}

void CMultiIndices::display() {
    m_impl.display();
}

void CMultiIndices::append(DataFile &df) {
    int batchsize = 10000;
    for (long i = 0; i < df.getSpectrumNum(); i += batchsize) {
        long newspecnum = 0;
        long start_spec_id = i, end_spec_id = i + batchsize > df.getSpectrumNum() ? df.getSpectrumNum() : i + batchsize;
        cout << "Scan:  " << start_spec_id << " - " << end_spec_id << "\t";
        float *vec = df.toFloatVector(m_dim, newspecnum, m_removeprecursor, m_useFlankingBins, m_tolerance, start_spec_id,
                                      end_spec_id);
        m_impl.add(vec, newspecnum);
        delete[] vec;
    }
}

int CMultiIndices::getNum() { return m_indexNum; }


CMultiIndices::~CMultiIndices() {}

// append init files
void CMultiIndices::appendList(string datafilelist) {
    vector<string> files = readlines(datafilelist);
    for (int i = 0; i < files.size(); i++) {
        string eachfile = files[i];
        cout << "\nprocessing file " << i + 1 << " / " << files.size() << endl;
        DataFile df(eachfile,0,-1);
        append(df);
        m_impl.display();
    }
    m_impl.write(); // append init files

}

CShuffledVector::CShuffledVector(long specnum) {
    m_shuffledIdx.assign(specnum, 0);
    iota(m_shuffledIdx.begin(), m_shuffledIdx.end(), 0);
    gcc5shuffle(m_shuffledIdx.begin(), m_shuffledIdx.end(), mt19937(12345));
}

void CShuffledVector::getVector(float *src, int offset, int d, int n, float *dest) {
    for (int k = 0; k < n; k++) {
        long position = m_shuffledIdx[offset + k] * d;
        copy(src + position, src + position + d, dest + k * d);
    }
}

CRecallStat::CRecallStat() : m_total(0), m_correct(0) {}

void CRecallStat::verify() { m_total++; }

void CRecallStat::correct() { m_correct++; }

double CRecallStat::recall() {
    if (m_total == 0) { return 0; }
    return m_correct * 1.0 / m_total;
}

void CRecallStat::show() {
    cout << "recall = " << recall() << endl;
}
