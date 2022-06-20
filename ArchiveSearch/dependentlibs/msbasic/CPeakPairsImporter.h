//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_DATASETTESTING_H
#define MYTOOL_DATASETTESTING_H

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <vector>
using namespace std;
class MGFReader;
class problem;
class feature_node;
class CPeakPair;
class COneHotEncodingAA;


// parse *.frag.txt file
// file format is :
// PSMID   PeptideSequence        charge intensity AA1 AA2 iontype postion intensity AA3 AA4 iontype position CLASS/LABEL
// 0       AAAAAAAAAAAAAAAGAGAGAK 2      903       A   G   y       3       1106      G   A   y       4        0
class CPeakPairsImporter {
    vector<CPeakPair> m_peakpairs;
    int startpos;
    string m_sourceFileName;
    string m_outFileStr;
    double m_minIntensityFC;

    const double MIN_PEAK_INTENSITY;
public:
    CPeakPairsImporter();
    CPeakPairsImporter(string filename);
    CPeakPairsImporter(const CPeakPairsImporter & other);
    ~CPeakPairsImporter();
    const vector<CPeakPair>  getsample() const;
    void swap_as(vector<CPeakPair> inputsample);
    void print();
    void getSubsetOnMinIntenFoldChange(CPeakPairsImporter &ds, double minIntensityFC);
    void getSubsetOnPepLen(CPeakPairsImporter &ds, int peplen);
    string getOutFileStr();
    void setSourceFileName(string sourceFileName);
    void setOutFileStr(string outFileStr);
    void setMinIntensityFC(double minIntensityFC);
    string exportTofile(string outfile);
    void toLiblinearProb(problem *prob, bool takeRange, COneHotEncodingAA &oheAA);
    int size();

    CPeakPair linearSearch(int sample_id);
    CPeakPair get_by_sample_id(int index);
    CPeakPair get_by_index(int index);
    // todo:  to be deleted
    int getNumOfUniqScanCounts();
    void getPeakPairSampleFrom(MGFReader &mgf_reader, bool useGhostPeaks);
    void add_feature_node(int k, vector<feature_node> &elements, double value) const;
    void convertSampleToFeatureNodes(CPeakPair &peakPair, int &k, vector<feature_node> &elements,
                                     vector<string> &featureNames, COneHotEncodingAA &oheAA) const;
};


#endif //MYTOOL_DATASETTESTING_H
