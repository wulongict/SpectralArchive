//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_PEAKPAIRSCORE_H
#define MYTOOL_PEAKPAIRSCORE_H

class PeakPairScore {
    double m_label;
    double m_positive_score;
    double m_negative_score;
public:
    PeakPairScore(double label, double pos_score, double neg_score);
    PeakPairScore(const PeakPairScore &other);
    PeakPairScore & operator = (const PeakPairScore &other);
    void print();
    int getlabel();
    double getpositive_score();
    double getnegative_score() const;
};

#endif //MYTOOL_PEAKPAIRSCORE_H
