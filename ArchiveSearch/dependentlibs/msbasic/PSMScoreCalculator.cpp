//
// Created by wulong on 11/25/17.
//

#include "PSMScoreCalculator.h"
#include <cmath>
#include "PeakPairScore.h"
#include "CPeakPair.h"
using namespace std;

PSMScoreCalculator::~PSMScoreCalculator() = default;

int PSMScoreCalculator::getPSMId() const {
    return m_psm_id;
}

void PSMScoreCalculator::setPSMId(int psmid) {
    m_psm_id = psmid;
}

void PSMScoreCalculator::append(const PeakPairScore& score, double groundtruthlabel, CPeakPair ts) {
    m_groundtruthlabel.push_back(groundtruthlabel);
    m_singlescore.push_back(score);
    m_trainingsamples.push_back(ts);
}

PSMScoreCalculator::PSMScoreCalculator() = default;

PSMScoreCalculator::PSMScoreCalculator(const PSMScoreCalculator &other) {
    m_psm_id = other.m_psm_id;
    m_singlescore = other.m_singlescore;
    m_combined_score = other.m_combined_score;
    m_groundtruthlabel = other.m_groundtruthlabel;
    m_trainingsamples = other.m_trainingsamples;

}

PSMScoreCalculator &PSMScoreCalculator::operator=(const PSMScoreCalculator &other) {
    m_psm_id = other.m_psm_id;
    m_singlescore = other.m_singlescore;
    m_combined_score = other.m_combined_score;
    m_groundtruthlabel = other.m_groundtruthlabel;
    m_trainingsamples = other.m_trainingsamples;
    return *this;
//        cout << "Here we are: assignment" << endl;
//        exit(0);
}

void PSMScoreCalculator::print() {
    cout << "psm id " << m_psm_id << endl;
    cout << "groundtruth num " << m_groundtruthlabel.size() << endl;
    cout << "sub score num " << m_singlescore.size() << endl;
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
//            cout << "inside loop" << endl;
        cout << m_groundtruthlabel.at(i) << " ";
        m_singlescore.at(i).print();
        m_trainingsamples.at(i).print();

//            cout << flush;
    }
}

void PSMScoreCalculator::collect(vector<vector<double>> &fc_score_pos) {
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        vector<double> tmp={log(m_trainingsamples[i].getFoldChange()), m_singlescore[i].getpositive_score() - 0.5};
        fc_score_pos.push_back(tmp);
    }
}

void PSMScoreCalculator::compute(const string& scoretype)// one should not call this function directly
{
    if (scoretype == "counts") {
        countsScore();
    } else if (scoretype == "prob") {
        probScore();
    } else if (scoretype == "prob2") {
        probScore2();
    } else if (scoretype == "logprob") {
        logprobScore();
    } else if (scoretype == "weightedprob") {
        weightedprobScore();
    } else if (scoretype == "minus") {
        minusprob();
    } else if (scoretype == "tprob") {
        tprob();
    } else {
        cout << "Error" << endl;
        exit(0);
    }
}

void PSMScoreCalculator::weightedprobScore() {
    m_combined_score = 0;
    // simple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        double foldchange = m_trainingsamples[i].m_intensity_of_peakA / m_trainingsamples[i].m_intensity_of_peakB;
        if (label == 1) {
            m_combined_score += (m_singlescore[i].getpositive_score() * foldchange);
        } else {
            m_combined_score += (m_singlescore[i].getnegative_score() / foldchange);

        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::logprobScore() {
    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == 1) {
            double pos = m_singlescore[i].getpositive_score();
            m_combined_score += (-log(pos));
        } else {
            double neg = m_singlescore[i].getnegative_score();
            m_combined_score += (-log(neg));

        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::probScore2() {
    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == 1) {
            double pos = m_singlescore[i].getpositive_score();
            m_combined_score += (pos * pos);
        } else {
            double neg = m_singlescore[i].getnegative_score();
            m_combined_score += (neg * neg);

        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::minusprob() {
    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == 1) {
            if (m_groundtruthlabel[i] == label) {
                m_combined_score += m_singlescore[i].getpositive_score();
            } else {
                m_combined_score -= m_singlescore[i].getnegative_score();
            }

        } else {
            if (m_groundtruthlabel[i] == label) {
                m_combined_score += m_singlescore[i].getnegative_score();
            } else {
                m_combined_score -= m_singlescore[i].getpositive_score();
            }


        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::tprob() {
    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == 1) {
            if (m_singlescore[i].getpositive_score() > 0.7)
                m_combined_score += m_singlescore[i].getpositive_score();
        } else {
            if (m_singlescore[i].getnegative_score() > 0.7)
                m_combined_score += m_singlescore[i].getnegative_score();

        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::probScore() {
    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == 1) {

            m_combined_score += m_singlescore[i].getpositive_score();
        } else {
            m_combined_score += m_singlescore[i].getnegative_score();

        }

    }
    m_combined_score /= m_groundtruthlabel.size();
}

void PSMScoreCalculator::countsScore() {

    m_combined_score = 0;
    // simaple score
    for (int i = 0; i < m_groundtruthlabel.size(); i++) {
        int label = m_groundtruthlabel[i];
        if (label == m_singlescore[i].getlabel()) {
            m_combined_score += 1;
        }
    }
    if (m_groundtruthlabel.empty()) {cout << "Error: empty ground truth label vector" << endl; exit(0);}
    m_combined_score /= m_groundtruthlabel.size();
}

double PSMScoreCalculator::getcombinedscore() const {
    return m_combined_score;
}
