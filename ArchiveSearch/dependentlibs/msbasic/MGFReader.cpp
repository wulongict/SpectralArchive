//
// Created by wulong on 4/11/17.
//

#include <mutex>
#include <string>
#include <fstream>
#include <utility>
#include <vector>
#include <random>
#include <spdlog/spdlog.h>
#include "MGFReader.h"
#include "CDebugMode.h"
#include "BasePSMFeature.h"
#include "ConcretePSMFeatures.h"
#include "PeakInfo.h"
#include "peak_annotation.h"
#include "CPeakPair.h"
#include "../../../librarymsms/CThreadsPool.h"
#include "../../../librarymsms/Util.h"
#include "../../../librarymsms/XMLFileParser.h"
#include "../../../External/SpectraST/SpectraSTPeakList.hpp"

class SpectraSTPeakListWrapper;

class SpectraSTPeakListWrapper{
    SpectraSTPeakList *pSPL;
    char *m_mgf_spec_buf;
    bool m_isLowMassAcc;
    string m_title;
    bool m_isParsed;
public:
    SpectraSTPeakListWrapper(){
        pSPL = nullptr;
        m_mgf_spec_buf = nullptr;
        m_isLowMassAcc = false;
        m_isParsed = false;
    }
    ~SpectraSTPeakListWrapper(){
        delete pSPL;
        pSPL = nullptr;
    }
    void set(char *mgf_spec_buf, bool isLowMassAcc){
        m_mgf_spec_buf = mgf_spec_buf;
        m_isLowMassAcc = isLowMassAcc;
    }
    string getM_title(){
        parse();
        return m_title;
    }
    SpectraSTPeakList * getSpectrumPtr(){
        parse();
        return pSPL;
    }

    void parse(){
        if(!m_isParsed){
            parseAsNeeded();
            m_isParsed = true;
        }
    }

private:
    void parseAsNeeded(){
        char *rest = m_mgf_spec_buf;
        char *p = strtok_r(m_mgf_spec_buf, "\r\n",&rest); // potential error memory leaking
        while (p != nullptr) {
            if(strstr(p, "BEGIN IONS")){
                pSPL = new SpectraSTPeakList(1.0, 1); // dummy spectra; parent mz=1.0; parent charge=1 todo
                string fragType = m_isLowMassAcc? "CID": "HCD";
                pSPL->setFragType(fragType);
            } else if(strstr(p, "PEPMASS")){
                char *q = strstr(p,"=");
                double pepmass = strtod(q+1, nullptr);
                pSPL->setParentMz(pepmass);
            } else if (strstr(p, "TITLE")){
                char *q = strstr(p, "=");
                m_title=string(q+1);
            } else if (strstr(p, "CHARGE")){
                char *q = strstr(p, "=");
                pSPL->setParentCharge(atoi(q+1));
            }else if (strstr(p, "END IONS")){
                if (pSPL == nullptr) {
                    cout << "Error:" << "Spectra without  (BEGIN IONS)" << endl;
                    exit(-1);
                }
                break;
            } else if(isdigit(*p))        {
                char *endptr;
                double mz=strtod(p,&endptr);
                double intensity = strtod(endptr, nullptr);

                pSPL->insert(mz, intensity, "", "");
            }
            p = strtok_r(rest, "\r\n",&rest);
        }
    }
};

class MGFLoaderThread: public ICThreadTask{
    SpectraSTPeakListWrapper *pSPLW;
public:
    explicit MGFLoaderThread(SpectraSTPeakListWrapper *p){
        pSPLW = p;
    }
    void run () override {
        pSPLW->parse();
    }
};

class MGFReaderX{
    vector<SpectraSTPeakListWrapper> m_spectra;
    char * m_buf;

public:
    string getTitle(int i ){return m_spectra[i].getM_title();}
    int getSpectrumNum(){return m_spectra.size();}
    SpectraSTPeakList * getSpectrumPtr(int i ){return m_spectra[i].getSpectrumPtr();}
    void releaseBuf(){
        if(m_buf!=nullptr) {
            delete [] m_buf;
            m_buf = nullptr;
        }
    }
    MGFReaderX(const string& filename, bool m_isLowMassAcc){
        SimpleTimer st("MGF Reader");
        string line;
        long len = 0;

        File::getfilesize(filename, len);
        FILE *pfile = fopen(filename.c_str(), "r");

        m_buf = new char[len+1];
        m_buf[len]='\0'; // fixed Dec 13 2020
        if(m_buf == nullptr){
            cout << "Fail to create buffer of " << len/1024/1024 << " MB" << endl;
            throw runtime_error("Fail to create buffer for mgf reader");
        }
        long nsize = fread(m_buf, sizeof(char), len, pfile);
        cout << "[Info] file loaded into buffer: "   << nsize*1.0/1024/1024 << "MB"<< endl;
        fclose(pfile);

        char *p="BEGIN IONS";
        for(char *q = strstr(m_buf, p); q!=nullptr; q=strstr(q+strlen(p),p)){
            m_spectra.emplace_back();
            m_spectra.back().set(q,m_isLowMassAcc);
        }
        loadSpectraMultThreads();
        cout << "[Info] number of spectra is " << m_spectra.size() << endl;
    }

    void loadSpectraMultThreads() {
        vector<ICThreadTask*> task_vec;
        for(auto & i : m_spectra){
            task_vec.push_back(new MGFLoaderThread(&i));
        }
        CTaskPool::issueTasks(task_vec,false, true, getProperThreads(),"parsing mgf ");
        for(auto & i : task_vec)  {
            delete i;
            i = nullptr;
        }
        releaseBuf();
    }

    ~MGFReaderX(){
       releaseBuf();
    }
};

//#include "../../../External/pepnovo/predicted_spectra.h"
//#include "../../../pepnovo/predicted_spectra.h"

class CAnnotationTask: public ICThreadTask{
    MGFReader *m_mgfr;
    int m_idx;
    string m_peptide;
    int m_charge;
public:
    CAnnotationTask(MGFReader *mgfr, int idx, int charge, string peptide): m_mgfr(mgfr){
        m_peptide = std::move(peptide);
        m_idx = idx;
        m_charge = charge;
    }
    void run() override{
        m_mgfr->annotate_one_spectrum(m_idx,m_peptide, m_charge);
    }
};

class PeptideStr{
    string m_pepstring;
public:
    explicit PeptideStr(string &pepstring):m_pepstring(pepstring){}
    PeptideStr & I_to_L(){
        replace(m_pepstring.begin(), m_pepstring.end(), 'I','L');
        return *this;
    }
    PeptideStr & strip_mod(){
        const regex e("\[[0-9\.]*\]");
        m_pepstring=regex_replace(m_pepstring, e, "");
        return *this;
    }

    PeptideStr & strip_charge(){
        m_pepstring = m_pepstring.substr(0,m_pepstring.find_first_of('/'));
        return *this;
    }
    string getPeptideString(){return m_pepstring;}
};


MGFReader::MGFReader(const string &filename, bool isLowMassAcc, bool is_truth_known) {
    m_truthKnown = is_truth_known;
    m_fixMz = false;
    m_isLowMassAcc = isLowMassAcc;
    m_filename = filename;

    m_lazy_mgfreader = make_shared<MGFReaderX>(filename, m_isLowMassAcc);
    cout << "[Info] " << getSpectraNum() << " spectra loaded from MGF file" << endl;
    peptides.assign(getSpectraNum(), nullptr);
    initPep2Idx();
}


vector<char *> getBeginPositions(char *buf, char *p){
    vector<char *> ret;
    for(char *q = strstr(buf, p); q!=nullptr; q=strstr(q+strlen(p), p)){
        ret.push_back(q);
    }
    return ret;
}

int countSubStringNum(char *buf, char *p){
    int count = 0;
    return getBeginPositions(buf, p).size();

}


MGFReader::~MGFReader() {
    releasePeptidePtrs();
}

//// the following function is not used!
//// It is not ready also. No return value.
//void MGFReader::generate_theoretical_spectra_from_title() {
//    cout << "To be used!!" << endl;
//
//    predicted_spectra ps;
//    for (int i = 0; i < m_titles.size(); ++i) {
//        cout << "[Info] Spectrum " << i << endl;
//        int seq_pos = getTitle(i).find_first_of("/");
//        if (seq_pos == string::npos) {
//            cout << "This is not a valid sequence format. Correct format is like: "
//                 << endl << " ABCDEFK/3" << endl;
//        }
//        string pepseq = getTitle(i).substr(0, seq_pos);
//
//        if (pepseq.find_first_of("n[]0123456789.") != string::npos) {
//            cout << "[Info] we found something " << endl;
//            continue;
//        }
//        int charge = Spectra[i]->getParentCharge();
//        cout << "peptide/charge: " << pepseq << "/" << charge << endl;
//        vector<Peaks> pks = ps.predict(pepseq, charge, 30);
//        cout << "---- theoretical fragmentation -----" << endl;
//        for (auto eachpeak: pks) eachpeak.print();
//    }
//}

void MGFReader::annotate_with_title(bool verbosity, bool debug, int scan, const string& debug_peptide, bool isTrainingSet) {
    cout << "[Info] Start annotating the spectra with title " << endl;
    if (CDebugMode::callDebug()->getMdebug()) verbosity = true;

    if (peptides.size() != getSpectraNum()) {
        cout << "Error: incorrect title and peptide size " << endl;
        exit(0);
    }
    releasePeptidePtrs();
    Peptide::defaultTables();
    vector<tuple<int, string, int>> workingindex; // index: pepseq, charge
    for(int i = 0; i < getSpectraNum(); ++i)  {
        string thisTitle = getTitle(i);
        string pepseq = PeptideStr(thisTitle).strip_charge().getPeptideString();
        if (CDebugMode::callDebug()->getMdebug()) {
            if (scan != i && string::npos == pepseq.find(debug_peptide)) // both not good scan and peptide
                continue;
            cout << "[DEBUG] Spectrum " << i << "-- peptide: " << getTitle(i) << endl;
        }
        int charge = getSpectrumPtr(i)->getParentCharge();
        setPeptidePtr(i,new Peptide(pepseq, charge));

        workingindex.emplace_back(i,pepseq, charge);
    }
    // refine
    if(isTrainingSet){
        refineWorkingIndex(workingindex);
    }

    // create the task
    vector<ICThreadTask *> task_vec;
    cout << "Total number of tasks "  << workingindex.size()<< endl;
    for(auto & i : workingindex){
        task_vec.push_back(new CAnnotationTask(this, get<0>(i), get<2>(i),get<1>(i)));
    }
    CTaskPool::issueTasks(task_vec,false,true,getProperThreads(),"Annotate with TITLE");
    for(int i = 0; i < workingindex.size(); i ++)    {
        delete task_vec[i];
        task_vec[i] = nullptr;
    }
    Peptide::deleteTables();
}

class CPsmCount{
    int incorrect;
    int correct;
public:
    CPsmCount(){
        incorrect = 0;
        correct=0;
    }
    void check(PSMInfo &psminfo, PSMInfo &psminfo2){
        if(psminfo.charge==psminfo2.charge and psminfo.end_scan == psminfo2.end_scan ){
            correct ++;
        }else{
            incorrect ++;
        }
    }
    void print() const{
        cout << "Correct | Incorrect " << endl << correct << " " << incorrect << endl;
    }
};

// For each psm in decoyfile,
// Attention: fdr_threshold not used!
void MGFReader::annotation_with_search_result(const string &decoyfile, int rank, bool debug,
                                              int scan, double fdr_threshold, const string& debug_peptide,bool isTrainingSet) {
    vector<tuple<int, string, int>> annotationTasks;
    releasePeptidePtrs();
    Peptide::defaultTables();

    vector<int> missing_scans;
    int num_spectra_with_no_hits_of_given_rank = 0;
    int num_decoy_PSMs_with_shared_peptides = 0;
    regex e("\[[0-9\.]*\]");
    CPsmCount psmStats;
    vector<int> ranks;
    if(rank>=0 and rank <= 4){
        ranks.assign(getSpectraNum(),rank);
    }
    if(rank==-1){
        ranks.assign(getSpectraNum(),0);
        int seed=42;
        std::mt19937_64 g(seed);
        int max_rank=3, min_rank=0;
        std::uniform_int_distribution<std::mt19937::result_type> dist(min_rank,max_rank);
        for(auto & x: ranks) x=dist(g);
    }
    {
        cout << "[Info] Annotation on search result " << decoyfile << endl;
        shared_ptr<CometPepXMLParser> crp=make_shared<CometPepXMLParser>(decoyfile);
        for(int i = 0; i < getSpectraNum(); i ++)   {
            PSMInfo psminfo;
            bool found = false;
            int scan_num = i + 1;
            if(CDebugMode::callDebug()->useScanCharge){
                int charge=getSpectrumPtr(i)->getParentCharge();
                found = crp->getPSMInfobyScanCharge(scan_num, charge, psminfo);
            }  else{
                found = crp->getPSMInfobyScan(scan_num,psminfo); // fixed: Dec 2020
            }

            if (found) {
                if (psminfo.searchhits.size() > ranks[i]) { // rank=0; get tophit!
                    string peptide_withmod = psminfo.searchhits[ranks[i]]->m_modified_peptide;
                    if (debug) {
                        if (scan != i && string::npos == peptide_withmod.find(debug_peptide)) {
                            continue;
                        } else {
                            cout << "[DEBUG] Spectrum " << i << "--peptide " << peptide_withmod << endl;
                        }
                    }
                    if(peptide_withmod=="NULL"){
                        cout << "[OK] Peptide == NULL " << i << endl;
                        psminfo.print();
                    }
                    // check this peptide against the gt
                    string plain_peptide = PeptideStr(peptide_withmod).strip_mod().I_to_L().getPeptideString();

                    if(not isTruthPeptide(plain_peptide, false)){
                        setPeptidePtr(i,new Peptide(peptide_withmod, psminfo.charge));

                        annotationTasks.emplace_back(i, peptide_withmod, psminfo.charge);
                    } else{
                        num_decoy_PSMs_with_shared_peptides ++;

                    }
                } else {
                    num_spectra_with_no_hits_of_given_rank ++;
                }
            } else { // not found
                missing_scans.push_back(scan_num);
            }
        }
        psmStats.print();
        // the crp object will be released here
        cout << "[Info] releasing CometPepXML Parser..." << endl;
        crp = nullptr;
    }

    cout << "[Info] number of spectra without peptide hit (with rank "<< rank << "): " << num_spectra_with_no_hits_of_given_rank << endl;
    cout << "[Info] number of decoy PSMs matched with shared peptides: " << num_decoy_PSMs_with_shared_peptides << endl;
    cout << "\n[Info] "<< missing_scans.size() <<" scans skipped by Search Engine (Comet): " << endl;
    std::copy(begin(missing_scans), end(missing_scans), std::ostream_iterator<int>(std::cout, " "));
    cout << endl;

    // refine working index for annotation
    if(isTrainingSet){
        // todo: if we do not annotate all of the spectra,it will fail to collect all of the fatures.
        refineWorkingIndex(annotationTasks);
    }

    vector<ICThreadTask *> task_vec;
    for(auto & annotationTask : annotationTasks)    {
        task_vec.push_back(new CAnnotationTask(this,get<0>(annotationTask),get<2>(annotationTask),get<1>(annotationTask)));
    }
    CTaskPool::issueTasks(task_vec,false, true, getProperThreads(), "Annotate with search results");
    for(int j = 0; j < annotationTasks.size(); j ++){
        delete task_vec[j];
        task_vec[j] = nullptr;
    }
    Peptide::deleteTables();
}


void MGFReader::annotate_one_spectrum(int i, const string &pepseq, int charge) {
    getSpectrumPtr(i)->setPeptidePtr(peptides[i]);
    getSpectrumPtr(i)->annotate(true, m_fixMz);
}

void MGFReader::print() {
    cout << "# Spectrum " << getSpectraNum() << endl;
    for (int i = 0; i < getSpectraNum(); ++i) {
        cout << "Spectrum " << i + 1 << endl;
        getSpectrumPtr(i)->printPeaks();
    }
}

void MGFReader::getWorkingIndex(const string& debug_peptide, int scan, vector<int> &workingIndex ){
    for (int j = 0; j < getSpectraNum(); ++j) {
        if (peptides[j] == nullptr) {
            continue;
        }

        if (CDebugMode::callDebug()->getMdebug()) {
            if (scan != j && string::npos == peptides[j]->fullWithCharge().find(debug_peptide)) {
                continue;
            }
            cout << "debug peptide: " << peptides[j]->fullWithCharge() << endl;
        }
        workingIndex.push_back(j);
    }
}

class FeatureExtractionTask: public ICThreadTask{
    int j;
    MGFReader * mgfr;
    vector<Feature*> *m_features;
    CMyMatrix *m_f_matrix;
public:
    FeatureExtractionTask(int _j, MGFReader *_mgfr,vector<Feature*> *features,CMyMatrix *f_matrix):
    j(_j),mgfr(_mgfr),m_features(features), m_f_matrix(f_matrix){}
    void run() override{
        int current_k ;
        try {
            PSMFeature pPSM(mgfr->getSpectrumPtr(j), mgfr->getPeptidePtr(j), j, mgfr->getMgfFileName());
            for (int k = 0; k < m_features->size(); ++k) {
                current_k = k;
                double feature_k = m_features->at(k)->calculate_feature(&pPSM);
                m_f_matrix->set(j, k, feature_k);
            }
        } catch (exception &e) {
            cout << "error caught: " << e.what() << endl << "j="<< j << "\t" << "mgfr size " << mgfr->getSpectraNum() << "\t feature "
            << m_features->at(current_k)->m_feature_name << endl;
        }
    }
};

struct range {
    double m_min;
    double m_max;

    range() {
        m_min = 1;
        m_max = 0;
    }
};

void MGFReader::export_features_with_PSM(vector<Feature *> features, const string& outputfilename, bool debug, int scan,
                                         const string& debug_peptide) {
    CMyMatrix mm(getSpectraNum(),features.size());
    mm.initialize(0);

    vector<int> workingIndex;
    getWorkingIndex(debug_peptide, scan, workingIndex);
    cout << "get working index: " << workingIndex.size() << endl;


    vector<ICThreadTask *> task_vec;
    for(int i = 0; i < workingIndex.size(); i ++){
        task_vec.push_back(new FeatureExtractionTask(workingIndex[i],this,&features,&mm));
    }
    cout << "task_vec init done" << endl;
    CTaskPool::issueTasks(task_vec,false, true, getProperThreads(),"Calcualte PSM Features");
    for(auto & task: task_vec){
        delete task;
    }
    cout << "Export feature to file: " << outputfilename << endl;
    ofstream fout;
    fout.open(outputfilename.c_str(), ios::out);
    // output header
    fout << "PeptideSeq\t";
    for (auto & feature : features) {
        fout << feature->m_feature_name << "\t";
    }
    fout << endl;
    vector<range> min_max(features.size());
    Progress ps2(getSpectraNum(), "Save feature to file");
    for (int spec_idx = 0; spec_idx < getSpectraNum(); spec_idx++) {
        ps2.increase();
        string pepstring="NULL";
        if (peptides[spec_idx] != nullptr) {
            if (CDebugMode::callDebug()->getMdebug()) {
                if (scan != spec_idx && string::npos == peptides[spec_idx]->fullWithCharge().find(debug_peptide)) continue;
            }
            pepstring=peptides[spec_idx]->full();
        }  else {
            if (CDebugMode::callDebug()->getMdebug()) {
                if (scan != spec_idx) continue;
            }
        }
        if(pepstring!="NULL" or CDebugMode::callDebug()->outputNullPep){
            fout << pepstring << "\t";
            for (int i = 0; i < features.size(); ++i) {
                double feature_i = mm.get(spec_idx, i);
                if(min_max[i].m_min>min_max[i].m_max){
                    min_max[i].m_min = feature_i;
                    min_max[i].m_max = feature_i;
                } else{
                    if(feature_i< min_max[i].m_min){
                        min_max[i].m_min = feature_i;
                    }
                    if(feature_i > min_max[i].m_max){
                        min_max[i].m_max = feature_i;
                    }
                }
                fout << feature_i << "\t";
            }
            fout << endl;
        }
    }
    fout.close();

    cout << "[Info] PSM feature saved to  " << outputfilename << endl;

    string featureRangeFile = outputfilename+".range";

    // In this part we try to rescale the feature values. from min max to 0 to 1. or -1 to 1.
    // we may do the same thing for the  features of another model. liblinear.
    fout.open(featureRangeFile.c_str(), ios::out);
    for(int i = 0; i < features.size(); i ++){
        fout << "Feature #" << i << "\t" << min_max[i].m_min << "\t" << min_max[i].m_max << "\t" << features[i]->m_feature_name << endl;
    }
    fout.close();

    // todo: dot product can be a better score.
}


void MGFReader::output_pepnovo_training_set(bool verbosity) {
    cout << "[Info] Start to write PepNovo training mgf file" << endl;
    string fileoutput = m_filename.substr(0, m_filename.size() - 4) + "_pepnovo_training.mgf";
    ofstream fout(fileoutput, ios::out);
    int export_counts = 0;
    Progress outputprogress( getSpectraNum());
    for (int i = 0; i < getSpectraNum(); ++i) {
        outputprogress.increase();
        string peptideseq = getTitle(i).substr(0, getTitle(i).size() - 2);
        if (regex_match(peptideseq, regex(".*[0-9]+.*"))) { // true if there is a match
            if (verbosity) { cout << "[Info] skip peptide " << peptideseq << endl; }
            continue;
        }
        if (verbosity) {
            cout << "[Info] loading spectrum of peptide " << peptideseq << endl;
        }

        fout << "BEGIN IONS" << endl;
        fout << "SEQ=" << peptideseq << endl;
        fout << "PEPMASS=" << setprecision(6) << fixed << getSpectrumPtr(i)->getParentMz() << endl;
        fout << "CHARGE=" << getSpectrumPtr(i)->getParentCharge() << endl;
        fout << "TITLE=" << getTitle(i) << endl;
        for (unsigned long j = 0; j < getSpectrumPtr(i)->getNumPeaks(); ++j) {
            Peak peak;
            getSpectrumPtr(i)->getPeak(j, peak);
            fout << setprecision(4) << fixed << peak.mz << " " << setprecision(1) << fixed << peak.intensity << endl;
        }
        fout << "END IONS" << endl;
        export_counts++;
    }
    fout.close();
    cout << "[Info] Number of spectra exported to training set: " << export_counts << endl;
}

void MGFReader::print_with_index(int i) { // from 1 to N  . NOT from zero
    if (i < 1 or i >getSpectraNum() ) {
        cout << "Error: MGFReader index out of bound: index range 1 ~ N" << endl;
        exit(0);
    }
    cout << "Spectrum " << i << endl;
    getSpectrumPtr(i-1)->printPeaks();
//    Spectra[i - 1]->printPeaks();
}

// todo: We could have a version without charge state.
void MGFReader::output_new_mgf(bool verbosity) {
    cout << "[Info] Start to write new mgf file" << endl;
    string fileoutput = m_filename.substr(0, m_filename.size() - 4) + "_new.mgf";
    ofstream fout(fileoutput, ios::out);
    int export_counts = 0;
    Progress outputprogress(getSpectraNum());
    for (int i = 0; i < getSpectraNum(); ++i) {
        outputprogress.increase();
        fout << "BEGIN IONS" << endl;
        fout << "PEPMASS=" << setprecision(6) << fixed << getSpectrumPtr(i)->getParentMz() << endl;
        fout << "CHARGE=" << getSpectrumPtr(i)->getParentCharge() << endl;
        fout << "TITLE=" << getTitle(i)  << endl;
        for (int j = 0; j < getSpectrumPtr(i)->getNumPeaks(); ++j) {
            Peak peak;
            getSpectrumPtr(i)->getPeak(j, peak);
            fout << setprecision(4) << fixed << peak.mz << " " << setprecision(1) << fixed << peak.intensity << endl;
        }
        fout << "END IONS" << endl;
        export_counts++;
    }
    fout.close();
    cout << "[Info] Number of spectra exported to training set: " << export_counts << endl;
}

void MGFReader::filter_min_mz(double min_mz) {
    Progress filterprogress(getSpectraNum());
    for (unsigned long i = 0; i < getSpectraNum(); ++i) {
        filterprogress.increase();
        //cout << "Spectrum " << i + 1 << endl;
        getSpectrumPtr(i)->simplify(9999, 999999, false, false, min_mz);
    }
}

void MGFReader::export_fragmentation_pattern(const string& outfile, bool zero_intensity_fragment_ions) {
    SimpleTimer st("Exporting fragmentation pattern");
    ofstream fout;
    fout.open(outfile.c_str(), ios::out);
    Progress ps(getSpectraNum());
    int counts = 0;

    cout << "[Info] #Peptide " << peptides.size() << endl;
    for (int i = 0; i < getSpectraNum() ; i++) {
        ps.increase();
        if (peptides[i] == nullptr) {
            continue;
        }
        counts += 1;
        export_feature_for_one_spectrum(zero_intensity_fragment_ions, fout, i);
    }
    fout.close();
    cout << "[Info] #Spectra processed: " << getSpectraNum() << endl;
    cout << "[Info] #Spectra annotated: " << counts << endl;

}

void MGFReader::getPeakPairs(bool useGhostPeaks, vector<CPeakPair> &peakPairs) {
    const int MAX_NUM_PSM_FOR_FRAG_PEAKPAIR = 20000;
    spdlog::get("A")->info("Max number of PSMs used for training fragmentation score: {}", MAX_NUM_PSM_FOR_FRAG_PEAKPAIR);
    vector<int> workingindex;
    workingindex.reserve(getSpectraNum());
    for(int i = 0; i < getSpectraNum(); i ++){
        if(peptides[i] != nullptr){
            workingindex.push_back(i);
        }
    }
    cout << "[Info] "<< workingindex.size() << " out of " << getSpectraNum() << " spectra are annotated"<< endl;
    if(workingindex.size() > MAX_NUM_PSM_FOR_FRAG_PEAKPAIR)
    {
        int seed = 42;
        gcc5shuffle(workingindex.begin(), workingindex.end(), mt19937(seed));
        workingindex.resize(MAX_NUM_PSM_FOR_FRAG_PEAKPAIR);
    }

    cout << "[Info] refined to " << workingindex.size() << " spectra" << endl;

    Progress ps(workingindex.size(), "Extract Peak Pair Features");
    for (int psm_id : workingindex) {
        ps.increase();
        get_feature_for_one_spectrum(useGhostPeaks, peakPairs, psm_id, peptides[psm_id], getSpectrumPtr(psm_id) );
    }
    cout << "[Info] #Spectra processed: " << workingindex.size() << endl;

}


// todo: with pointer of pka, copy n paste is faster.
// 2018-01-15
void collect_peak_info_for_single_spectrum(SpectraSTPeakList *pkl, vector<PeakInfo> &pifs) {
    int pknum = pkl->getNumPeaks();
    for (int j = 0; j < pknum; j++) {
        _peak pk;
        pkl->getPeak(j, pk);
        vector<peak_annotation> pka;
        parse_annotation_as_structX(pk.annotation, pka);
//        parse_annotation_as_struct(pk.annotation, pka);
        pifs.emplace_back(pk.mz, pk.intensity, pka);
    }
}

void MGFReader::export_feature_for_one_spectrum(bool zero_intensity_fragment_ions, ofstream &fout, int i) {
    string pepstr;
    int charge;
    vector<PeakInfo> base_type_no_NL_charge1_pifs;
    vector<bool> bionFound;
    vector<bool> yionFound;
    beforePeakPairFeatureExtraction(pepstr, charge, base_type_no_NL_charge1_pifs, bionFound, yionFound,
                                    peptides[i], getSpectrumPtr(i));


    // todo: to make this code be symmetric.
    // a--b   ---> b--a
    //
//    cout << "[Info] Exporting PeakPairs of spectrum " << i   << flush << endl;
    for (int k = 0; k < base_type_no_NL_charge1_pifs.size(); k++) {
        for (int l = k + 1; l < base_type_no_NL_charge1_pifs.size(); l++) {
            // compare a and b
            PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];
            PeakInfo pif2 = base_type_no_NL_charge1_pifs[l];
            export_for_peakinfo_pair(fout, i, pepstr, charge, pif1, pif2);
        }
    }

    // consider GHOST peaks
    bool consider_non_existed_peak = zero_intensity_fragment_ions;

    if (consider_non_existed_peak) {
        export_feature_of_ghost_peaks(fout, i, pepstr, charge, base_type_no_NL_charge1_pifs, bionFound, yionFound);
    }
}

void
beforePeakPairFeatureExtraction(string &pepstr, int &charge, vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                vector<bool> &bionFound, vector<bool> &yionFound, Peptide *peptidePtr,
                                SpectraSTPeakList *specpkl) {
    pepstr = peptidePtr->stripped;
    charge = peptidePtr->charge;
    if (bionFound.empty()) { bionFound.assign(pepstr.length(), false); }
    if (yionFound.empty()) { yionFound.assign(pepstr.length(), false); }
    vector<PeakInfo> pifs;
    collect_peak_info_for_single_spectrum(specpkl, pifs);
    collect_charge1_b_y_ion_info(pifs, base_type_no_NL_charge1_pifs);
    get_each_b_y_ion_existence(pifs, bionFound, yionFound);
}

void get_feature_for_one_spectrum(bool useGhostPeaks, vector<CPeakPair> &peakPairs, int psm_id, Peptide *ipeptide,
                                  SpectraSTPeakList *specpkl) {
    string pepstr;
    int charge;
    vector<PeakInfo> base_type_no_NL_charge1_pifs;
    vector<bool> bionFound;
    vector<bool> yionFound;
    beforePeakPairFeatureExtraction(pepstr, charge, base_type_no_NL_charge1_pifs, bionFound, yionFound,
                                    ipeptide, specpkl);

    for (int k = 0; k < base_type_no_NL_charge1_pifs.size(); k++) {
        for (int l = k + 1; l < base_type_no_NL_charge1_pifs.size(); l++) {
            // compare a and b
            PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];
            PeakInfo pif2 = base_type_no_NL_charge1_pifs[l];
            get_peakinfo_pair(peakPairs, psm_id, pepstr, charge, pif1, pif2);
        }
    }

    // consider GHOST peaks
    if (useGhostPeaks) {
        get_feature_of_ghost_peaks(peakPairs, psm_id, pepstr, charge, base_type_no_NL_charge1_pifs, bionFound, yionFound);
    }
}

void MGFReader::export_feature_of_ghost_peaks(ofstream &fout, int i, const string &pepstr, int charge,
                                              const vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                              const vector<bool> &bionFound,
                                              const vector<bool> &yionFound) const {// Here let's add those peak with intensity 0 ZERO
    for (int k = 0; k < base_type_no_NL_charge1_pifs.size(); k++) {
        for (int l = 1; l < bionFound.size(); l++) {
            if (bionFound.at(l) and yionFound.at(l)) continue;
            if (not bionFound.at(l)) {
                export_feature_for_by_ion_ghost(i, charge, k, l, fout, pepstr, base_type_no_NL_charge1_pifs, "b");


            }
            if (not yionFound.at(l)) {
                export_feature_for_by_ion_ghost(i, charge, k, l, fout, pepstr, base_type_no_NL_charge1_pifs, "y");

            }
        }
    }
}

void get_feature_of_ghost_peaks(vector<CPeakPair> &PPF, int psm_id, const string &pepstr, int charge,
                                const vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                const vector<bool> &bionFound,
                                const vector<bool> &yionFound) {
    for (int k = 0; k < base_type_no_NL_charge1_pifs.size(); k++) {
        for (int l = 1; l < bionFound.size(); l++) {
            if (bionFound.at(l) or yionFound.at(l)) continue;
            if (not bionFound.at(l)) {
                get_feature_for_by_ion_ghost(psm_id, charge, k, l, PPF, pepstr, base_type_no_NL_charge1_pifs, "b");
            }
            if (not yionFound.at(l)) {
                get_feature_for_by_ion_ghost(psm_id, charge, k, l, PPF, pepstr, base_type_no_NL_charge1_pifs, "y");
            }
        }
    }

}

void MGFReader::export_feature_for_y_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                               const vector<PeakInfo> &base_type_no_NL_charge1_pifs) const {
    PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];
//                    PeakInfo pif2 = base_type_no_NL_charge1_pifs[l];
    int pos1 = pif1.m_pka[0].pos;
    int pos2 = l;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = "y";
    if (base_ion_type1 == "y") { pos1 = pepstr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = pepstr.length() - pos2; }
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " "
         << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " " << (pif1.m_intensity > 0.0) << endl;

    // on the other way
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " "
         << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " " << (pif1.m_intensity < 0.0) << endl;
}

void MGFReader::export_feature_for_b_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                               const vector<PeakInfo> &base_type_no_NL_charge1_pifs) const {
    PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];
//                    PeakInfo pif2 = base_type_no_NL_charge1_pifs[l];
    int pos1 = pif1.m_pka[0].pos;
    int pos2 = l;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = "b";
    if (base_ion_type1 == "y") { pos1 = pepstr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = pepstr.length() - pos2; }
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " "
         << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " " << (pif1.m_intensity > 0.0) << endl;
    // on the other way
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " "
         << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " "

         << (pif1.m_intensity < 0.0) << endl;
}

// todo : missed another direction in normal mode! fix it
void MGFReader::export_feature_for_by_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                                const vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                                string iontype) const {
    PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];
//                    PeakInfo pif2 = base_type_no_NL_charge1_pifs[l];
    int pos1 = pif1.m_pka[0].pos;
    int pos2 = l;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = iontype; // iontype should be b or y
    if (iontype != "b" and iontype != "y") {
        cout << "Error: iontype must be either b or y" << endl;
        exit(0);
    }
    if (base_ion_type1 == "y") { pos1 = pepstr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = pepstr.length() - pos2; }
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " "
         << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " " << (pif1.m_intensity > 0.0) << endl;
    // on the other way
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << "0.0" << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2 << " "
         << l << " "
         << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " "
         << base_ion_type1 << " " << pif1.m_pka[0].pos << " "

         << (pif1.m_intensity < 0.0) << endl;
}

void get_feature_for_by_ion_ghost(int psm_id, int charge, int k, int l, vector<CPeakPair> &peakPairs,
                                  const string &peptideStr,
                                  const vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                  const string& iontype) {
    PeakInfo pif1 = base_type_no_NL_charge1_pifs[k];

    int pos1 = pif1.m_pka[0].pos;
    int pos2 = l;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = iontype; // iontype should be b or y
    if (iontype != "b" and iontype != "y") {
        cout << "Error: iontype must be either b or y" << endl;
        exit(0);
    }
    if (base_ion_type1 == "y") { pos1 = peptideStr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = peptideStr.length() - pos2; }

    CPeakPair x, y;
    x.m_PSM_ID = psm_id;
    x.m_peptide = peptideStr;
    x.m_precursor_charge = charge;
    x.m_intensity_of_peakA = pif1.m_intensity;
    x.m_left_aa_of_peakA = peptideStr[pos1 - 1];
    x.m_right_aa_of_peakA = peptideStr[pos1];
    x.m_Aby = base_ion_type1[0];
    x.m_Apos = pif1.m_pka[0].pos;
    x.m_intensity_of_peakB = 0.0;
    x.m_left_aa_of_peakB = peptideStr[pos2 - 1];
    x.m_right_aa_of_peakB = peptideStr[pos2];
    x.m_Bby = base_ion_type2[0];
    x.m_Bpos = l;//pif2.m_pka[0].pos;
    x.m_y = (pif1.m_intensity > 0);
    x.m_sample_id = peakPairs.size();
    peakPairs.push_back(x);
    // on the other way.
    y.m_PSM_ID = psm_id;
    y.m_peptide = peptideStr;
    y.m_precursor_charge = charge;
    y.m_intensity_of_peakA = 0.0;//pif2.m_intensity;
    y.m_left_aa_of_peakA = peptideStr[pos2 - 1];
    y.m_right_aa_of_peakA = peptideStr[pos2];
    y.m_Aby = base_ion_type2[0];
    y.m_Apos = l;//pif2.m_pka[0].pos;
    y.m_intensity_of_peakB = pif1.m_intensity;
    y.m_left_aa_of_peakB = peptideStr[pos1 - 1];
    y.m_right_aa_of_peakB = peptideStr[pos1];
    y.m_Bby = base_ion_type1[0];
    y.m_Bpos = pif1.m_pka[0].pos;
    y.m_y = (0.0 > pif1.m_intensity);
    y.m_sample_id = peakPairs.size();
    peakPairs.push_back(y);
}

// here is it!
void MGFReader::export_for_peakinfo_pair(ofstream &fout, int i, const string &pepstr, int charge, const PeakInfo &pif1,
                                         const PeakInfo &pif2) const {
    int pos1 = pif1.m_pka[0].pos;
    int pos2 = pif2.m_pka[0].pos;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = pif2.m_pka[0].ion_base_type;
    if (base_ion_type1 == "y") { pos1 = pepstr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = pepstr.length() - pos2; }
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " " << base_ion_type1
         << " " << pif1.m_pka[0].pos << " "
         << pif2.m_intensity << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2
         << " " << pif2.m_pka[0].pos << " " << (pif1.m_intensity > pif2.m_intensity) << endl;

    // on the other way.
    fout << i << " "; // spectrum ID unique
    fout << pepstr << " ";
    fout << charge << " ";
    fout << pif2.m_intensity << " " << pepstr[pos2 - 1] << " " << pepstr[pos2] << " " << base_ion_type2
         << " " << pif2.m_pka[0].pos << " "
         << pif1.m_intensity << " " << pepstr[pos1 - 1] << " " << pepstr[pos1] << " " << base_ion_type1
         << " " << pif1.m_pka[0].pos << " " << (pif2.m_intensity > pif1.m_intensity) << endl;

}


void get_peakinfo_pair(vector<CPeakPair> &peakPairs, int psm_id, string peptideStr, int charge, PeakInfo &pif1,
                       PeakInfo &pif2) {
    int pos1 = pif1.m_pka[0].pos;
    int pos2 = pif2.m_pka[0].pos;
    string base_ion_type1 = pif1.m_pka[0].ion_base_type;
    string base_ion_type2 = pif2.m_pka[0].ion_base_type;
    if (base_ion_type1 == "y") { pos1 = peptideStr.length() - pos1; }
    if (base_ion_type2 == "y") { pos2 = peptideStr.length() - pos2; }

    CPeakPair x, y;
    x.m_PSM_ID = psm_id;
    x.m_peptide = peptideStr;
    x.m_precursor_charge = charge;
    x.m_intensity_of_peakA = pif1.m_intensity;
    x.m_left_aa_of_peakA = peptideStr[pos1 - 1];
    x.m_right_aa_of_peakA = peptideStr[pos1];
    x.m_Aby = base_ion_type1[0];
    x.m_Apos = pif1.m_pka[0].pos;
    x.m_intensity_of_peakB = pif2.m_intensity;
    x.m_left_aa_of_peakB = peptideStr[pos2 - 1];
    x.m_right_aa_of_peakB = peptideStr[pos2];
    x.m_Bby = base_ion_type2[0];
    x.m_Bpos = pif2.m_pka[0].pos;
    x.m_y = (pif1.m_intensity > pif2.m_intensity)?1:-1;
    x.m_sample_id = peakPairs.size();
    peakPairs.push_back(x);
    // on the other way.
    y.m_PSM_ID = psm_id;
    y.m_peptide = peptideStr;
    y.m_precursor_charge = charge;
    y.m_intensity_of_peakA = pif2.m_intensity;
    y.m_left_aa_of_peakA = peptideStr[pos2 - 1];
    y.m_right_aa_of_peakA = peptideStr[pos2];
    y.m_Aby = base_ion_type2[0];
    y.m_Apos = pif2.m_pka[0].pos;
    y.m_intensity_of_peakB = pif1.m_intensity;
    y.m_left_aa_of_peakB = peptideStr[pos1 - 1];
    y.m_right_aa_of_peakB = peptideStr[pos1];
    y.m_Bby = base_ion_type1[0];
    y.m_Bpos = pif1.m_pka[0].pos;
    y.m_y = (pif2.m_intensity > pif1.m_intensity)?1:-1;
    y.m_sample_id = peakPairs.size();
    peakPairs.push_back(y);


}

void get_each_b_y_ion_existence(const vector<PeakInfo> &pifs, vector<bool> &bionFound,
                                vector<bool> &yionFound) {
    for (auto x: pifs) {
        if (x.m_pka.size() == 1
            and x.m_pka[0].ion_NL == "0"
            and x.m_pka[0].charge == 1
            and x.m_pka[0].isotopic == false
            and x.m_pka[0].ion_base_type != "others") {
            if (x.m_pka[0].ion_base_type == "b") {
                int b_pos = x.m_pka[0].pos;
                bionFound[b_pos] = true;
            } else if (x.m_pka[0].ion_base_type == "y") {
                int y_pos = x.m_pka[0].pos;
                yionFound[y_pos] = true;
            }
        }

    }
}

void collect_charge1_b_y_ion_info(const vector<PeakInfo> &pifs, vector<PeakInfo> &base_type_no_NL_charge1_pifs) {
    for (auto x: pifs) {
        if (x.m_pka.size() == 1
            and x.m_pka[0].ion_NL == "0"
            and x.m_pka[0].charge == 1
            and x.m_pka[0].isotopic == false
            and x.m_pka[0].ion_base_type != "others")// charge 1+; single match; and not neutral loss
        {
            base_type_no_NL_charge1_pifs.push_back(x);
        }
    }
}


string MGFReader::getMgfFileName() {
    return m_filename;
}

void MGFReader::getSpectra(int spec_id, vector<tuple<double, double, string>> &onespec) {
    int peaknum = getSpectrumPtr(spec_id)->getNumPeaks();
    for (int i = 0; i < peaknum; ++i) {
        _peak peak;
        getSpectrumPtr(spec_id)->getPeak(i, peak);
        onespec.emplace_back(peak.mz, peak.intensity, peak.annotation);
    }

}

int MGFReader::getSpectraNum() {
    return m_lazy_mgfreader->getSpectrumNum();

}

#include <random>
// function not used!!!
void MGFReader::refineWorkingIndex(vector<tuple<int, string, int>> &workingindex) {
    cout << "Working index is NOT to be refined. " << endl;
    return;
    // todo: make sure the working index is OK...
    int seed = 42;
    const int MAX_SAMPLE_NUM = 100000;
    if(workingindex.size() > MAX_SAMPLE_NUM)    {
        gcc5shuffle(workingindex.begin(), workingindex.end(), mt19937(seed));
        workingindex.resize(MAX_SAMPLE_NUM);
        cout << "Number of spectra to be annotated: " << workingindex.size() << endl;
    }
}

void MGFReader::releasePeptidePtrs() {
    for(auto & peptide : peptides)
    {
        if (peptide != nullptr) {
            delete peptide;
            peptide = nullptr;
        }
    }

}


void MGFReader::initPep2Idx() {
    if(m_truthKnown)    {
        cout << "[Info] Creating peptide2Idx mapping. 1) remove modifications; 2) replace 'I' with 'L'" << endl;
        regex e("\[[0-9\.]*\]");

        for (int i = 0; i < getSpectraNum(); i ++)   {
            string thisTitle = getTitle(i);
            string plain_peptide = PeptideStr(thisTitle).strip_charge().strip_mod().I_to_L().getPeptideString();

            if(m_pep2idx.count(plain_peptide)==0){
                m_pep2idx[plain_peptide]={i};
            }else{
                m_pep2idx[plain_peptide].push_back(i);
            }
        }
        int num_shared_psm = 0;
        for(auto & it : m_pep2idx){
            num_shared_psm += it.second.size();
        }
        cout << "Indexed #peptide: " << m_pep2idx.size() << " #PSM " << num_shared_psm << endl;
    }
}

bool MGFReader::isTruthPeptide(string &plain_peptide, bool verbose) {
    bool ret = false;
    if(m_truthKnown and m_pep2idx.count(plain_peptide)>0){
        ret = true;

        if(verbose){
            vector<int> &idx = m_pep2idx[plain_peptide];

            cout << "[Shared Peptides] Found " << plain_peptide << " in ground truth set : "<< endl;
            for(auto eachIdx : idx){
                cout << eachIdx << " " << getTitle(eachIdx) << endl;
            }
        }

    }
    return ret;
}

SpectraSTPeakList *MGFReader::getSpectrumPtr(int i) {
    return m_lazy_mgfreader->getSpectrumPtr(i);
}

string MGFReader::getTitle(int i) {
    return m_lazy_mgfreader->getTitle(i);

}
