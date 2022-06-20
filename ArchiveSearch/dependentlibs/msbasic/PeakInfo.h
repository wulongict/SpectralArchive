//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_PEAKINFO_H
#define MYTOOL_PEAKINFO_H

#include <vector>
#include <iostream>

using namespace std;
class peak_annotation;

using namespace std;

class PeakInfo {
public:
    PeakInfo(double mz, double inten, vector<peak_annotation> pka);
    PeakInfo(const PeakInfo &other);
    void print();
    ~PeakInfo();
    double m_mz;
    double m_intensity;
    vector<peak_annotation> m_pka;
};

#endif //MYTOOL_PEAKINFO_H
