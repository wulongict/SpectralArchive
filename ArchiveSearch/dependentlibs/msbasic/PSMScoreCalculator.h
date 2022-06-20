//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_LIBLINEARPSMSCORETESTING_H
#define MYTOOL_LIBLINEARPSMSCORETESTING_H


#include <vector>
#include <string>
using namespace std;
class PeakPairScore;
class CPeakPair;

class PSMScoreCalculator {
public:
    PSMScoreCalculator();
    PSMScoreCalculator(const PSMScoreCalculator &other);
    PSMScoreCalculator & operator = (const PSMScoreCalculator & other);
    ~PSMScoreCalculator();
    int getPSMId() const;
    void setPSMId(int psmid);
    void append(const PeakPairScore& score, double groundtruthlabel, CPeakPair ts);
    void collect(vector<vector<double>> &fc_score_pos);
    double getcombinedscore() const;
    void compute(const string& scoretype);
    void print();
private:
    void weightedprobScore();
    void logprobScore();
    void probScore2();
    void minusprob();
    void tprob();
    void probScore();
    void countsScore();

private:
    vector<double> m_groundtruthlabel; // what is groundtruth label?
    int m_psm_id{};
    vector<PeakPairScore> m_singlescore;
    double m_combined_score{}; // could be many different scores calculated by different method
    vector<CPeakPair> m_trainingsamples;
};


#endif //MYTOOL_LIBLINEARPSMSCORETESTING_H
