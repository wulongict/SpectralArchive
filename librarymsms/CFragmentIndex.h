//
// Created by wulong on 8/2/19.
//

#ifndef MYTOOL_CFRAGMENTINDEX_H
#define MYTOOL_CFRAGMENTINDEX_H



#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
using namespace std;
class CMzFileReader;
class ICMzFile;
class ICId2RankInten;


struct CPeakRecord
{
    long spectrumId;
    double precursorMz;
    int precursorCharge;
    double fragmentIntensity;

    CPeakRecord(long specid, double parentMz, int chg, double fragInten);
    CPeakRecord(const CPeakRecord &other);

    void output(ofstream &fout) const;
};


struct FragBin {
    vector<CPeakRecord> entry;
    FragBin()= default;

    void addPeak(long specid, double parentMz, int chg, double fragInten) {
        entry.emplace_back(specid, parentMz, chg, fragInten);
    }
    void output(ofstream &fout)   {
        fout << "\t" << entry.size() << endl;
        for (auto & i : entry)        {
            i.output(fout);
        }
    }
};

class CFragBinConvertor {
    uint16_t m_bin;
    double m_fragbin_to_mz;
public:
    explicit CFragBinConvertor(uint16_t bin): m_bin(bin){
        m_fragbin_to_mz = 2000.0/UINT16_MAX;
    }
    double toMz() const{
        return m_bin*m_fragbin_to_mz;
    }
};


// FragIndex
// The mz range [0, 2000] of MS2 spectra was splited into 65535 bins [FragBin]
// For each FragBin, we create a list of peak records [CPeakRecords]
// each contains mz, intensity, precursorMz, spectrumId of each single peaks.
//  |---------------------------------------------------|
// 0|---------------------------------------------------|65535
// 0|---------------------------------------------------|2000.0 mz
//          mz1         mz2         mz3
//          |           |           |
//          |           V           |
//          V                       V
//
class CFragIndex
{
    shared_ptr<vector<FragBin>> m_allMzBins_ptr;
    vector<float> m_norm;
//    shared_ptr<vector<float>> m_norm_ptr;
    int m_PeakNum;
    long m_bgSpectraNum;
    vector<int> m_nonZeroPeakNum;
    shared_ptr<ICId2RankInten> m_id2rank;
public:
    explicit CFragIndex(int PeakPerSpec);

    void buildIndex(shared_ptr<ICMzFile> &mzfile, vector<long> &bgSpectraIds);
    void dpscores(uint16_t *spec, vector<int> &score, bool normalize, int tol);

    void output(const string& filename, long maxlen=-1);
    long getSpecNum() const;
    int getSizeAllMzBin();
    float getNorm(int i);
    void setNorm(float norm, long bgspecId);
    void setNonZeroPeakNum(shared_ptr<ICMzFile> &mzfile, vector<long> &bgSpectraIds);
    int getNonZeroPeakNum(int bgspecId){return m_nonZeroPeakNum[bgspecId];}
    shared_ptr<ICId2RankInten> getId2Inten();
};

#endif //MYTOOL_CFRAGMENTINDEX_H
