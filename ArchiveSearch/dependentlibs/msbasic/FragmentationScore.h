//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_LIBLINEARPSMRESULTANALYSIS_H
#define MYTOOL_LIBLINEARPSMRESULTANALYSIS_H

#include <string>
#include <vector>
using namespace std;

class classifier;
class CPeakPairsImporter;
class CFragScore;
class PSMScoreCalculator;

void calculatePSMscore(const string &scoretype, CPeakPairsImporter &sample_Data, classifier &lr);

class FragmentationScore {
    string m_samplefilename;
    string m_psmid_filename;
    string m_predicted_result;
    vector<PSMScoreCalculator> m_PSMscore;
    vector<double> m_groundtruthlabel;
    vector<int> m_spectrum_id;
    vector<int> m_SampleID;
    CFragScore *m_predictionFile;
    string m_PSM_scoretype;
    void load_psmid();
    void load_ground_truth();
    void groupSampleByPSM(CPeakPairsImporter &input);
    void calcScore(const string& scoretype);
    void print();
public:

    FragmentationScore(CFragScore * predictionFile);
    FragmentationScore(CPeakPairsImporter &data, string & scoreType, classifier &lr);
    ~FragmentationScore();
    FragmentationScore(const FragmentationScore &other);
    FragmentationScore &operator=(const FragmentationScore &other);

    string compute(const string &scoretype, CPeakPairsImporter &SampleData);
    string exportScore(const string& scoretype);
    double getscore(int psmid);
};




#endif //MYTOOL_LIBLINEARPSMRESULTANALYSIS_H
