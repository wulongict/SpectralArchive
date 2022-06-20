//
// Created by wulong on 11/25/17.
//

#include "gnuplot_functions.h"
#include "CFragScore.h"
#include "PeakPairScore.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

CFragScore::CFragScore() {
    m_filename = "";
    m_header="";
    m_peakPairScore=vector<PeakPairScore>();

}

CFragScore::CFragScore(const CFragScore &other) {
    m_filename = other.m_filename;
    m_header = other.m_header;
    m_peakPairScore = other.m_peakPairScore;
}

CFragScore &CFragScore::operator=(const CFragScore &other) {
    m_filename = other.m_filename;
    m_header = other.m_header;
    m_peakPairScore = other.m_peakPairScore;
    return *this;
}

CFragScore::CFragScore(const string& filename) {
    m_filename = filename;
    ifstream fin;
    fin.open(filename, ios::in);
    if (!fin){cout << "Fail to read file " << filename << endl; exit(0);}
    std::getline(fin, m_header);
    double result, score_pos, score_neg;
    while (fin >> result >> score_pos >> score_neg) {
        m_peakPairScore.emplace_back(result, score_pos, score_neg);
    }
    fin.close();
}

void CFragScore::print() {
    cout << m_filename << endl << m_header << endl;
    if(m_peakPairScore.empty())    {
        cout << "No score available" << endl;
        exit(0);
    }  else{
        cout << "example: label, pos, neg " << endl;
        m_peakPairScore[0].print();
    }
}

PeakPairScore CFragScore::getScore(int index) {
    return m_peakPairScore.at(index);
}

const string &CFragScore::getM_filename() const {
    return m_filename;
}

int CFragScore::size() {
    return m_peakPairScore.size();
}

void CFragScore::add(const PeakPairScore& x) {
    m_peakPairScore.push_back(x);
}


void calculate_score_for_each_PSM(const string& testingfile) {

    string predictfile = testingfile + ".predict";
    CFragScore LLR(predictfile);

    string psmidfile = testingfile + ".psm_id";
    string outfile = testingfile + ".psmscore";
    ifstream fin_predict, fin_spectrum_id, fin_testingdata;
    ofstream fout_PSMscore;

    fin_predict.open(predictfile.c_str(), ios::in);
    fin_spectrum_id.open(psmidfile.c_str(), ios::in);
    fin_testingdata.open(testingfile.c_str(), ios::in);

    fout_PSMscore.open(outfile.c_str(), ios::out);

    // skip the first line
    string label;
    int postive, negative;
    fin_predict >> label >> postive >> negative;

    int last_psmid = -1;
    double score = 0;
    int num_score = 0;
    int correct_counts = 0;
    int incorrect_counts = 0;
    vector<double> score_plot;

    while (!(fin_predict.eof() or fin_spectrum_id.eof() or fin_testingdata.eof())) {
        int result;
        int spectrum_id;
        double score_pos, score_neg;
        fin_predict >> result >> score_pos >> score_neg;
        fin_spectrum_id >> spectrum_id;
        if (last_psmid != spectrum_id) {
            if (num_score != 0) {
                double normalized_score = score / sqrt(num_score);
                fout_PSMscore << last_psmid << " " << score << " "
                              << num_score << " " << normalized_score
                              << " " << correct_counts << " " << incorrect_counts << endl;
//                score_plot.push_back(normalized_score);
                score_plot.push_back(score / num_score);
            }
            num_score = 0;
            score = 0;
            last_psmid = spectrum_id;
            correct_counts = 0;
            incorrect_counts = 0;
        }


        string groundtruth;
        std::getline(fin_testingdata, groundtruth);
        int testingset_label = atoi(groundtruth.substr(0, groundtruth.find_first_of(' ')).c_str());

        if (testingset_label == 1) {
            score -= log(score_pos);
        } else {
            score -= log(score_neg);
//            score += (score_neg*score_neg);
        }
        if (testingset_label == result) {
            correct_counts++;
        } else {
            incorrect_counts++;
        }

        num_score++;

    }
    fout_PSMscore << last_psmid << " " << score << " "
                  << num_score << " " << score / num_score
                  << " " << correct_counts << " " << incorrect_counts << endl;
    fout_PSMscore.close();
    fin_predict.close();
    fin_spectrum_id.close();
    fin_testingdata.close();

    // plot a figure here
    gnuplot_histogram(testingfile, score_plot);


}
