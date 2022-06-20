//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_TRAININGSAMPLETEST_H
#define MYTOOL_TRAININGSAMPLETEST_H

#include <iostream>
#include <vector>
#include <map>

using namespace std;

static const double epsilon = 1e-7;
static const double DIVIDED_BY_ZERO_INF = 2.1;

class CPeakPair {
public:
    CPeakPair();

    CPeakPair(const CPeakPair &other);

    CPeakPair &operator=(const CPeakPair &other);

    ~CPeakPair();


    int getSampleID();

    void setSampleID(int sample_id);

    double getFoldChange();

    friend ostream &operator<<(ostream &fout, CPeakPair &ts);// bug fixed
    friend istream &operator>>(istream &fin, CPeakPair &ts);

    void print();

    int getAAcounts(char aa);

    int m_PSM_ID;
    string m_peptide;
    int m_precursor_charge;
    double m_intensity_of_peakA;
    char m_left_aa_of_peakA;
    char m_right_aa_of_peakA;
    char m_Aby;
    int m_Apos;
    double m_intensity_of_peakB;
    char m_left_aa_of_peakB;
    char m_right_aa_of_peakB;
    char m_Bby;
    int m_Bpos;
    int m_y;
    int m_sample_id;
    vector<int> m_AAcounts;
};

class COneHotEncodingAA{
    map<char, vector<int>> aaToOneHotEncoding;
    string m_legalAAs;
    vector<int> m_invalidAA;
    map<char, int> m_invalid_attempts;
public:
    COneHotEncodingAA();
    vector<int>& getEncoding(char c, bool verbose=false);
    string getlegalAAs();
    ~COneHotEncodingAA(){
//        cout << "Invalid symbols and counts" << endl;
//        for(auto &each: m_invalid_attempts){
//            cout << each.first << "\t" << each.second << endl;
//        }
//        cout << endl;
    }

};
void getOneHotEncodingAA(map<char, vector<int>> &aaToOneHotEncoding);

void export_trainingsample(const vector<CPeakPair> &sample, const string &outfile);

#endif //MYTOOL_TRAININGSAMPLETEST_H
