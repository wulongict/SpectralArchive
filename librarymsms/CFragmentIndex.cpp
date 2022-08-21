//
// Created by wulong on 8/2/19.
//

#include <numeric>
#include <random>
#include "CFragmentIndex.h"
#include "Util.h"
#include "ProteomicsDataTypes.h"
#include "CMzFileReader.h"

void CFragIndex::output(const string& filename, long maxlen) {
    ofstream fout(filename, ios::out);
    maxlen = maxlen<=0 or maxlen >= getSizeAllMzBin()  ? getSizeAllMzBin() : maxlen;
    for (int i = 0; i < getSizeAllMzBin()  and i < maxlen; i++)   {
        fout << i ;
        m_allMzBins_ptr->at(i).output(fout);
    }
}

CFragIndex::CFragIndex(int PeakPerSpec) {
    m_PeakNum = PeakPerSpec;
    m_id2rank = make_shared<CId2RankInten>(m_PeakNum, false);
    int len = UINT16_MAX;
    m_allMzBins_ptr=make_shared<vector<FragBin>>(len+1, FragBin());
    m_bgSpectraNum = 0;
}

void CFragIndex::buildIndex(shared_ptr<ICMzFile> &mzfile, vector<long> &bgSpectraIds) {
    m_bgSpectraNum = bgSpectraIds.size();
    m_norm.assign(m_bgSpectraNum, 0);
//    Progress ps(m_bgSpectraNum,"building fragment index");
    for(long i = 0; i < bgSpectraIds.size(); i ++)    {
//        ps.increase();
        float norm_i = m_id2rank->getNorm(mzfile->getSpecBy(bgSpectraIds.at(i)));
        setNorm(norm_i,i);
        CMzSpec spec = mzfile->getMzSpec(bgSpectraIds.at(i));
        int pkn = spec.getPeakNum();

        for(int j = 0; j < pkn; j ++)   {
            int intensity = m_PeakNum - j;
            int bin = spec.getMz(intensity);// to be changed
            m_allMzBins_ptr->at(bin).addPeak(i, 0.0, 0.0, intensity);
        }
    }
}

void CFragIndex::dpscores(uint16_t *spec, vector<int> &score, bool normalize, int tol) {
    score.assign(m_bgSpectraNum, 0);
    for(int i = 0; i < m_id2rank->getPeakPerSpec(); i ++)  {
        int mz_bin = spec[i];
        int rankIntensity = m_id2rank->toIntensity(i);
        int left_bin = mz_bin - tol >= 0 ? mz_bin - tol : 0;
        int right_bin = mz_bin + tol < getSizeAllMzBin() ? mz_bin + tol : getSizeAllMzBin() - 1;
        for(int k = left_bin; k <= right_bin; k ++)  {
            FragBin &fragEntry = (*m_allMzBins_ptr)[k]; // for all of the spec with frag k
            for (auto & j : fragEntry.entry)  {
                CPeakRecord *spIdx = &j; // spectra record
                double fragIntensity= spIdx->fragmentIntensity;
                score[spIdx->spectrumId] += fragIntensity * rankIntensity;
            }
        }
    }


    if(normalize) {
        float querynorm = m_id2rank->getNorm(spec);
        const int MAX_SCORE = m_id2rank->getMaxScore();
        for (long i = 0; i < score.size(); i++) {
            double s = score[i]/getNorm(i)/querynorm*MAX_SCORE;
            score[i] = floor(s) > MAX_SCORE ? MAX_SCORE : floor(s);
        }
    }
}

int CFragIndex::getSizeAllMzBin() {
    return m_allMzBins_ptr->size();
}

long CFragIndex::getSpecNum() const {return m_bgSpectraNum;}

shared_ptr<ICId2RankInten> CFragIndex::getId2Inten() {return m_id2rank;}

void CFragIndex::setNorm(float norm, long bgspecId) {
    m_norm[bgspecId]=norm;
}

void CFragIndex::setNonZeroPeakNum(shared_ptr<ICMzFile> &mzfile, vector<long> &bgSpectraIds) {
    m_nonZeroPeakNum.assign(bgSpectraIds.size(),0);
    for(long i = 0; i < bgSpectraIds.size(); i ++){
        m_nonZeroPeakNum[i]=mzfile->getPeakNum(bgSpectraIds[i]);
    }
}

float CFragIndex::getNorm(int i) {
    return m_norm[i];
}

void CPeakRecord::output(ofstream &fout) const {
    fout << spectrumId << "\t" << fragmentIntensity << endl;
}

CPeakRecord::CPeakRecord(const CPeakRecord &other) {
    spectrumId = other.spectrumId;
    precursorMz = other.precursorMz;
    precursorCharge = other.precursorCharge;
    fragmentIntensity = other.fragmentIntensity;
}

CPeakRecord::CPeakRecord(long specid, double parentMz, int chg, double fragInten) {
    spectrumId = specid;
    precursorMz = parentMz;
    precursorCharge = chg;
    fragmentIntensity = fragInten;
}
