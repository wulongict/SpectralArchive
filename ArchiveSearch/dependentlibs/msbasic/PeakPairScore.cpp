//
// Created by wulong on 11/25/17.
//

#include "PeakPairScore.h"
#include <iostream>
using namespace std;

PeakPairScore::PeakPairScore(double label, double pos_score, double neg_score) {
    m_label = label;
    m_positive_score = pos_score;
    m_negative_score = neg_score;
}

PeakPairScore::PeakPairScore(const PeakPairScore &other) {
    m_label = other.m_label;
    m_positive_score = other.m_positive_score;
    m_negative_score = other.m_negative_score;

}

PeakPairScore &PeakPairScore::operator=(const PeakPairScore &other) {
    m_label = other.m_label;
    m_positive_score = other.m_positive_score;
    m_negative_score = other.m_negative_score;
    return *this;
}

void PeakPairScore::print() {

    cout << m_label << " " << m_positive_score << " " << m_negative_score << endl;
}

int PeakPairScore::getlabel() {
    return int(m_label);
}

double PeakPairScore::getpositive_score() {
    return m_positive_score;
}

double PeakPairScore::getnegative_score() const {
    return m_negative_score;
}
