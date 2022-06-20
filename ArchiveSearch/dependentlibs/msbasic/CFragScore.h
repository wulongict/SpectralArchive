//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_LIBLINEARRESULTREADERTESTING_H
#define MYTOOL_LIBLINEARRESULTREADERTESTING_H

#include <string>
#include <vector>
#include "PeakPairScore.h"
using namespace std;
class PeakPairScore;

class CFragScore {
private:
    string m_filename;
    string m_header;  // first line of the file
    vector<PeakPairScore> m_peakPairScore;
public:
    CFragScore();
    // why we need the two following functions
    CFragScore(const CFragScore & other);
    CFragScore & operator = (const CFragScore & other);
    CFragScore(const string& filename);
    void print();
    PeakPairScore getScore(int index);
    void add(const PeakPairScore& x);
    const string &getM_filename() const;
    int size();
};


void calculate_score_for_each_PSM(const string& testingfile);

#endif //MYTOOL_LIBLINEARRESULTREADERTESTING_H
