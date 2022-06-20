//
// Created by wulong on 11/25/17.
//

#include "PeakInfo.h"
#include "peak_annotation.h"
PeakInfo::PeakInfo(double mz, double inten, vector<peak_annotation> pka) {
    m_mz = mz;
    m_intensity = inten;
    m_pka = pka;
}

PeakInfo::PeakInfo(const PeakInfo &other) {
    m_mz = other.m_mz;
    m_intensity = other.m_intensity;
    m_pka = other.m_pka;
}

void PeakInfo::print() {
    cout << m_mz << " " << m_intensity << endl;
    peak_annotation::print(m_pka);
}

PeakInfo::~PeakInfo() = default;
