//
// Created by wulong on 11/25/17.
//

#include "CPeakPairsImporter.h"
#include "CPeakPair.h"
#include <fstream>
#include <algorithm>
#include "MGFReader.h"
#ifdef LIBLINEAR_MT
#include "../../../External/liblinear-multicore-2.42/linear.h"
#else
#include "../../../External/liblinear-2.11/linear.h"
#endif

#include <boost/timer/timer.hpp>
#include <utility>
#include "../../../librarymsms/Util.h"

using namespace std;
using namespace boost::timer;


const vector<CPeakPair> CPeakPairsImporter::getsample() const {
    return m_peakpairs;
}

CPeakPairsImporter::CPeakPairsImporter(const CPeakPairsImporter &other):MIN_PEAK_INTENSITY(other.MIN_PEAK_INTENSITY) {
    m_peakpairs = other.getsample();
    startpos = -1;
    m_sourceFileName = other.m_sourceFileName;
    m_minIntensityFC = other.m_minIntensityFC;
    m_outFileStr = other.m_outFileStr;
}

CPeakPairsImporter::CPeakPairsImporter():MIN_PEAK_INTENSITY(200) {
    startpos = -1;
    m_peakpairs = vector<CPeakPair>();
    m_sourceFileName = "";
    m_minIntensityFC = 1.0; // equal size
    m_outFileStr = "";
}

CPeakPairsImporter::CPeakPairsImporter(string filename):CPeakPairsImporter() {
    m_sourceFileName = std::move(filename);
    cout << "[Info] Loading file " << m_sourceFileName << endl;
    startpos = -1;
    auto_cpu_timer t;
    ifstream fin;
    fin.open(m_sourceFileName.c_str(), ios_base::in);
    if (!fin) {
        cout << "[Info] Error! Can not open file: " << m_sourceFileName << endl;
        exit(0);
    }
    int sample_counts = 0;
    while (!fin.eof()) {
        if (sample_counts % 1000000 == 0) {
            cout << t.format() << "\r";
        }
        CPeakPair ts;
        ts.setSampleID(sample_counts);
        sample_counts++;
        fin >> ts;
        if (fin.eof()) break;
        m_peakpairs.push_back(ts);
    }
    cout << "read sample number: " << m_peakpairs.size() << endl;
    fin.close();
}

void CPeakPairsImporter::swap_as(vector<CPeakPair> inputsample) {
    m_peakpairs.swap(inputsample);
    startpos = -1;
}

CPeakPairsImporter::~CPeakPairsImporter() = default;

void CPeakPairsImporter::getSubsetOnMinIntenFoldChange(CPeakPairsImporter &ds, double minIntensityFC) {
    m_minIntensityFC = minIntensityFC;
    m_outFileStr = m_sourceFileName.substr(0, m_sourceFileName.length() - 4) + "_" + to_string(minIntensityFC) + "_";
    double min_intensity = 200;
    // ToDo: ds and sample_subset keep one of them
    // ToDo: There is a bug here: what if m_intensity_of_peakA is ZERO
    vector<CPeakPair> sample_subset(m_peakpairs.size());
    auto it = copy_if(m_peakpairs.begin(), m_peakpairs.end(), sample_subset.begin(),
                      [minIntensityFC, min_intensity](const CPeakPair& ts) {
                          if ((ts.m_intensity_of_peakA > ts.m_intensity_of_peakB * minIntensityFC and ts.m_intensity_of_peakA > min_intensity) ||
                              (ts.m_intensity_of_peakB > ts.m_intensity_of_peakA * minIntensityFC and ts.m_intensity_of_peakB > min_intensity))
                              return true;
                          else return false;
                      });
    sample_subset.resize(distance(sample_subset.begin(), it));
    ds.swap_as(sample_subset);
    ds.setMinIntensityFC(m_minIntensityFC);
    ds.setSourceFileName(m_sourceFileName);
    ds.setOutFileStr(m_outFileStr);
}

void CPeakPairsImporter::getSubsetOnPepLen(CPeakPairsImporter &ds, int peplen) {
    m_outFileStr = m_outFileStr + "sparse_len" + to_string(peplen);
    vector<CPeakPair> sample_subset(m_peakpairs.size(), CPeakPair());

    auto it = copy_if(m_peakpairs.begin(), m_peakpairs.end(), sample_subset.begin(),
                      [peplen](const CPeakPair& ts) { return ts.m_peptide.length() == peplen; });
    sample_subset.resize(distance(sample_subset.begin(), it));
    cout << "Subset: " << sample_subset.size() << endl;
    ds.swap_as(sample_subset);// I guess this is where the problem comes from
    ds.setOutFileStr(m_outFileStr);
    ds.setMinIntensityFC(m_minIntensityFC);
    ds.setSourceFileName(m_sourceFileName);

}
struct Range{
    double m_min;
    double m_max;
    Range(){
        m_min = 1;
        m_max = 0;
    }
    void refresh(double value){
        if(m_min>m_max){
            m_min = m_max = value;
        } else{
            if(value < m_min) m_min = value;
            else if (value > m_max) m_max = value;
        }

    }

};

void CPeakPairsImporter::toLiblinearProb(problem *prob, bool takeRange, COneHotEncodingAA &oheAA) {
    vector<string> featureNames;

    int max_index = 0;
    for (int i = 0; i < prob->l; i++) {
        CPeakPair &peakPair = m_peakpairs[i];

        if (peakPair.m_y > 0)
            prob->y[i] = 1; // use 0,1 instead of -1 +1
        else
            prob->y[i] = -1;// ion type of peak A

        int k = 1;
        vector<feature_node> elements;
        try{
            convertSampleToFeatureNodes(peakPair, k, elements, featureNames, oheAA);
        } catch (exception &ex) {
            cout << "caught error: " << ex.what() << endl;
            peakPair.print(); // Error because of the unknown amino acids: UZXOBJ; to continue;
            exit(0);
        }

        if(elements.empty())  {
            cout << "[Error] empty elements list" << endl;
        }

        prob->x[i] = nullptr;
        void * ptr = nullptr;

        ptr =  malloc(
                (elements.size() + 1) * sizeof(struct feature_node));// new feature_node[elements.size()+1];
        if (ptr == nullptr) {
            cout << "invalid pointer prob->x[i] : " << __FILE__ << ":" << __LINE__ << endl;
            exit(0);
        }
        else{
            prob->x[i] = (feature_node *)ptr;
        }

        for (int e = 0; e < elements.size(); e++) {
            prob->x[i][e].index = elements[e].index; // copy construct
            prob->x[i][e].value = elements[e].value;
        }

        prob->x[i][elements.size()].index = -1; // the last one should be -1

        if (k > max_index)
            max_index = k;

    }

    prob->n = max_index;

    if(takeRange){
        // find the range of the features
        cout << "[Info] Number of peak pair features collected: " << FormatWithCommas(prob->l) << endl;
        vector<Range> vec_ranger(prob->n);
        const int MAX_SAMPLE_INSPECTED = 20000;
        for(int i = 0; i < prob->l and i < MAX_SAMPLE_INSPECTED; i ++) {
            int previousid = 0;
            for(feature_node *a = prob->x[i]; a->index!=-1; a++)  {
                vec_ranger[a->index].refresh(a->value);
                if(a->index>previousid+1){
                    for(int k = previousid +1; k < a->index; k ++){
                        vec_ranger[k].refresh(0);
                    }
                }
            }

        }

        cout << "[Info] min max summary" << endl;
        for(int i = 1; i < vec_ranger.size(); i ++){
            cout << "Feature #" << i << "\t" << flush <<  featureNames[i-1] << "\t"<< vec_ranger[i].m_min << "\t" << vec_ranger[i].m_max << endl;
        }

    }

}

void
CPeakPairsImporter::convertSampleToFeatureNodes(CPeakPair &peakPair, int &k, vector<feature_node> &elements,
                                                vector<string> &featureNames, COneHotEncodingAA &oheAA) const {
    const bool doInitFeatureNames = featureNames.empty();

    string legalAAs = oheAA.getlegalAAs();
    if(doInitFeatureNames) featureNames.push_back("is_Peak_A_of_type_Y");
    // peak A is_Y_ion
    if (peakPair.m_Aby == 'y') {
        add_feature_node(k, elements, 1);
    }
    k++;

    if(doInitFeatureNames) featureNames.push_back("position_Of_Peak_A");
    // A Cleavage Site position
    // position of peak A
    double z = peakPair.m_Apos * 1.0 / peakPair.m_peptide.length();
    add_feature_node(k, elements, z);
    k++;

    // Amino acid to the left of A
    // one-hot encoding of amino acid to the left of A
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_left_aa_of_peakA); // aa_to_binary[peakPair.m_left_aa_of_peakA];
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_Peak_A_with_left_AA_"+string(1,legalAAs[i]));
        if (y != 0) {
            add_feature_node(k, elements, y);
        }
        k++;
    }

    // Amino acid to the right of A
    // one-hot encoding of amino acid to the left of A
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_right_aa_of_peakA); //aa_to_binary[peakPair.m_right_aa_of_peakA];
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_Peak_A_with_right_AA_"+string(1,legalAAs[i]));
        if (y != 0) {
            add_feature_node(k, elements, y);
        }
        k++;
    }

    if(doInitFeatureNames) featureNames.push_back("is_Peak_B_of_type_Y");
    // ion type of peak B
    // Peak B is_Y_ion
    if (peakPair.m_Bby == 'y') {
        add_feature_node(k, elements, 1);
    }
    k++;


    if(doInitFeatureNames) featureNames.push_back("position_Of_Peak_B");
    // B Amino Acid position
    // relative position of peak B
    add_feature_node(k, elements, peakPair.m_Bpos * 1.0 / peakPair.m_peptide.length());
    k++;

    // Amino acid to the left of B
    // one-hot encoding of amino acid to the left of B
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_left_aa_of_peakB);// aa_to_binary.at(peakPair.m_left_aa_of_peakB);
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_Peak_B_with_left_AA_"+string(1,legalAAs[i]));
        if (y != 0) {
            add_feature_node(k, elements, y);
        }
        k++;
    }

    // Amino acid to the right of B
    // one-hot encoding of amino acid to the right of B
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_right_aa_of_peakB);//aa_to_binary.at(peakPair.m_right_aa_of_peakB);
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_Peak_B_with_right_AA_"+string(1,legalAAs[i]));
        if (y != 0) {
            add_feature_node(k, elements, y);
        }
        k++;
    }


    // terminal AA
    // N-terminal AA
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_peptide.at(peakPair.m_peptide.size() - 1));// aa_to_binary.at(peakPair.m_peptide.at(peakPair.m_peptide.size() - 1));
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_nterm_AA_"+string(1,legalAAs[i]));
        if (y != 0) {
            add_feature_node(k, elements, y);
        }
        k++;
    }
    // C-terminal AA
    for(int i = 0; i < legalAAs.size(); i ++){
        auto encoding_vec = oheAA.getEncoding(peakPair.m_peptide.at(0));//aa_to_binary.at(peakPair.m_peptide.at(0));
        auto y = encoding_vec[i];

        if(doInitFeatureNames) featureNames.push_back("is_cterm_AA_"+string(1,legalAAs[i]));
        if (y != 0) { add_feature_node(k, elements, y); }
        k++;
    }


    // RHK the mobility related AAs
    string mobility_related_aas = "RHK";
    for(char &aa: mobility_related_aas){
        if(doInitFeatureNames) featureNames.push_back("proportion_Of_AA_"+string(1,aa));
        add_feature_node(k, elements, peakPair.getAAcounts(aa) * 1.0 / peakPair.m_peptide.length());
        k++;
    }

    if(doInitFeatureNames) featureNames.push_back("proportion_Of_RHK");
    add_feature_node(k, elements, peakPair.getAAcounts('R') + peakPair.getAAcounts('H') + peakPair.getAAcounts('K'));
    k++;


}

void CPeakPairsImporter::add_feature_node(int k, vector<feature_node> &elements, double value) const {
    feature_node u;
    u.index = k;
    u.value = value;
    elements.push_back(u);
}

string CPeakPairsImporter::exportTofile(string outfile) {
    ofstream fout;
    string sep = "\t";
    fout.open(outfile.c_str(), ios::out);

    ofstream psmid_fout;
    string psmid_filename = outfile + ".psm_id";
    psmid_fout.open(psmid_filename.c_str(), ios::out);

    auto_cpu_timer t;
    map<char, vector<int>> aa_to_binary;
    getOneHotEncodingAA(aa_to_binary);

    cout << "Peak pairs to be exported: " << m_peakpairs.size() << endl;
    Progress ps(m_peakpairs.size(), "Extracting peak pairs");
    int real_counts = 0;
    for (auto x: m_peakpairs) {
        int k = 1;
        ps.increase();

        if (x.m_y == 1) { fout << "+1" << sep; }
        else { fout << "-1" << sep; }

        if (x.m_Aby == 'y') { fout << k << ":" << "1" << sep; }
        k++;

        fout << k << ":" << x.m_Apos * 1.0 / x.m_peptide.length() << sep;
        k++;

        for (auto y : aa_to_binary[x.m_left_aa_of_peakA]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }
        for (auto y : aa_to_binary[x.m_right_aa_of_peakA]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }


        if (x.m_Bby == 'y') { fout << k << ":" << "1" << sep; }
        k++;

        fout << k << ":" << x.m_Bpos * 1.0 / x.m_peptide.length() << sep;
        k++;

        for (auto y : aa_to_binary[x.m_left_aa_of_peakB]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }
        for (auto y : aa_to_binary[x.m_right_aa_of_peakB]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }

        // terminal AA
        // Nterminal
        for (auto y : aa_to_binary[x.m_peptide[x.m_peptide.size() - 1]]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }
        // Cterminal
        for (auto y : aa_to_binary[x.m_peptide[0]]) {
            if (y != 0) { fout << k << ":" << y << sep; }
            k++;
        }

        // RHK the mobility related AAs
        fout << k << ":" << x.getAAcounts('R') * 1.0 / x.m_peptide.length() << sep;
        k++;

        fout << k << ":" << x.getAAcounts('H') * 1.0 / x.m_peptide.length() << sep;
        k++;

        fout << k << ":" << x.getAAcounts('K') * 1.0 / x.m_peptide.length() << sep;
        k++;

        fout << k << ":" << x.getAAcounts('R') + x.getAAcounts('H') + x.getAAcounts('K') << sep;
        k++;

        fout << endl;
        psmid_fout << x.m_PSM_ID << " " << x.m_sample_id << endl;
        real_counts++;
    }
    fout.close();
    psmid_fout.close();

    return outfile;
}

int CPeakPairsImporter::size() {
    return m_peakpairs.size();
}

CPeakPair CPeakPairsImporter::linearSearch(int sample_id) {
    if (startpos == -1) { startpos = 0; }
    for (int i = startpos; i < m_peakpairs.size(); i++) {
        if (m_peakpairs[i].m_sample_id == sample_id) {
            startpos = i + 1;
            return m_peakpairs.at(i);
        }
    }
    cout << "Error Not found" << endl;
    exit(0);
}

CPeakPair CPeakPairsImporter::get_by_sample_id(int index) {
    cout << "# samples: " << m_peakpairs.size() << endl;
    auto it = find_if(m_peakpairs.begin(), m_peakpairs.end(), [index](const CPeakPair ts) {
        return (ts.m_sample_id == index);
    });
    if (it == m_peakpairs.end()) {
        cout << "not found " << endl;
        exit(0);
    } else {
        return *it;
    }
    for (auto & m_peakpair : m_peakpairs) {
        if (m_peakpair.m_sample_id == index)
            return m_peakpair;
    }
}

void CPeakPairsImporter::print() {
    cout << "startpos: " << startpos << endl;
    cout << "m_sourceFileName: " << m_sourceFileName << endl;
    cout << "m_outFileStr: " << m_outFileStr << endl;
    cout << "m_minIntensityFC: " << m_minIntensityFC << endl;
    cout << "sample size: " << m_peakpairs.size() << endl;

}

void CPeakPairsImporter::getPeakPairSampleFrom(MGFReader &mgf_reader, bool useGhostPeaks) {

    mgf_reader.getPeakPairs(useGhostPeaks, m_peakpairs);


}

int CPeakPairsImporter::getNumOfUniqScanCounts() {
    // todo: 03152020 for some for the spectra, no feature found
    for (int i = 0; i < m_peakpairs.size(); i ++)  { }
    return 0;
}

CPeakPair CPeakPairsImporter::get_by_index(int index) {
    return m_peakpairs[index];
}

void CPeakPairsImporter::setMinIntensityFC(double minIntensityFC) {
    m_minIntensityFC = minIntensityFC;
}

void CPeakPairsImporter::setOutFileStr(string outFileStr) {
    m_outFileStr = std::move(outFileStr);
}

void CPeakPairsImporter::setSourceFileName(string sourceFileName) {
    m_sourceFileName = std::move(sourceFileName);
}

string CPeakPairsImporter::getOutFileStr() {
    return m_outFileStr;
}
