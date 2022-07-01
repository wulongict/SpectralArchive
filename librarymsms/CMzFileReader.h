//
// Created by wulong on 8/1/18.
//

#ifndef COMPACTSPECMZ_H
#define COMPACTSPECMZ_H

#include "ICMzFile.h"


using namespace std;
class DataFile;
class CSpectrum;
class PSMInfo;
class PeptideProphetParser;

class SpectraFilter{
	string m_pepxml;
	shared_ptr<PeptideProphetParser> m_ppp;
	double m_fdr_threshold;
	vector<int> m_scans;
	vector<PSMInfo> m_psms;
public:
	SpectraFilter(string pepxml);
	~SpectraFilter();
	bool skip(CSpectrum *spec);
};

struct MzSpecInfo{
    long residx;
    int peaknum;  // not included!
    float rt;
    float precursormz;
    int charge;
    int scan;
    int ms2idx;
    MzSpecInfo(long _residx, int _peaknum, float _rt, float _parentMz, int _chg, int _scan, int _ms2idx);
    MzSpecInfo(const MzSpecInfo & other){
        residx = other.residx;
        peaknum = other.peaknum;
        rt = other.rt;
        precursormz = other.precursormz;
        charge = other.charge;
        scan = other.scan;
        ms2idx = other.ms2idx;
    }
    MzSpecInfo();
    // todo: this input and output step do not has peaknum included
    friend ostream & operator<<(ostream &osm, const MzSpecInfo & sf);
    friend istream & operator>>(istream &ism, MzSpecInfo &sf);
    long getResIdx() const;
    double getPrecursorNeutralMass() const;
};


// *.mzXML.scan TSV file
// format
// IDX      ms2Idx    scan      preMz      charge     RT
// 0        0           2       698.417     1       1.5715
class CScanFile {
    string m_datafilename;
    string m_scanFilename;
    long m_resIdxOffset; // previous
public:
    vector<MzSpecInfo> infos;
    CScanFile(string _filename);
    CScanFile(const CScanFile & other);
    ~CScanFile(){}
    long getLastQueryIDX();
    MzSpecInfo & binarySearch(long queryidx);
    string getFileName();
    int getScanNum() const;
    bool isCreated();
    void create(long &cnts, bool verbose);
    void readLines(istream &fin, long len);
    void write_specinfo();
    void write_specinfo(ofstream &fout);
    // create scan header and export to filename
    void init();

    bool load();

    void fixResIdx();
};

// a file that are corresponding to  multiple scan files.
class CScanInfo {
    long m_lastResIdxOffset ;
    vector<CScanFile> m_readers;  // corresponding to more thatn one scan file
    string m_mzXMLListfile;
    string m_scanfilename;
public:
    CScanInfo();
    CScanFile & getScanFileReader(int fileid);
    MzSpecInfo &getMzSpecInfo(long queryidx, string &filename);
    ~CScanInfo();
    void load(string m_mzxml_filelist, bool isList, bool verbose);
    void appendFile(string file);
    int getFileNum();
    void appendSpecInfo(int _peaknum, float rt, float mz, int charge, int scan, int ms2count);
    long getTotalScanNum() const;
    void exportToCombinedFile();
    void clear();
private:
    void loadFiles(const vector<string>& files);
    // void setCombinedFilename(string mzxmlfiles);
    string getCombinedScanFilename(){return m_mzXMLListfile + ".scan";}
    void loadCombinedFile();
    bool isEmpty();
    int getFileId(long queryidx);
    string getFileName(int fileid);


    bool initWithCombinedFile();

};

class CMzFileReader: public ICMzFile {
    uint16_t *m_mzs;
    long m_specnum;
    int m_PeakNumPerSpectrum;
    string m_mzxml_filelist;
    vector<int> m_peaknum;
    CScanInfo m_scanInfo; // the only place where CScanInfo is used.
    double m_localMaxHalfWidth;
    int m_minPeakNum;
    shared_ptr<vector<float>> m_norm2_ptr;
    string m_mzFileName;

public:
    CMzFileReader(string mzxml_filelist, bool overwrite, bool islist, bool rmParentIon,
                  double localMaxHalfWidth, int minPeakNum, bool verbose, string mzFileName="");
	CMzFileReader(DataFile &df, bool overwrite, bool rmParentIon, SpectraFilter *sf,
                  double localMaxHalfWidth, int minPeakNum, bool verbose);            //CScanInfo

	// make shared
	static shared_ptr<CMzFileReader> makeShared(string mzxml_filelist, bool overwrite, bool islist, bool rmParentIon,
                                                double localMaxHalfWidth, int minPeakNum, bool verbose){
	    return make_shared<CMzFileReader>(mzxml_filelist, overwrite, islist, rmParentIon, localMaxHalfWidth, minPeakNum, verbose);
	}
	static shared_ptr<CMzFileReader> makeShared(DataFile &df, bool overwrite, bool rmParentIon, SpectraFilter *sf,
                                                double localMaxHalfWidth, int minPeakNum, bool verbose){
	    return make_shared<CMzFileReader>(df, overwrite, rmParentIon, sf, localMaxHalfWidth, minPeakNum, verbose);
	}
	virtual ~CMzFileReader();

	// functions to be removed!!
	string getListFile();
    string getMzFilename();

    
    void append(DataFile &df, bool rmParentIon, SpectraFilter *sf);  // CScanInfo
    void refresh_mz_object_from_disk();
    

    MzSpecInfo & getSpecMetaInfo(long queryidx, string &filename);   //CScanInfo

    long getSpecNum() override;
    string getClassName() override;
    int getPeakNum(long queryindex) override;
    float getSquaredNorm(uint16_t *p) override;
    float getSquaredNorm(long queryindex) override;
    uint16_t *getSpecBy(long queryindex) override;
    void calcDotProduct(int TopNPeak, int tol, uint16_t *queryspec, int blockSize, vector<long> &indexlist, vector<int> &scores) override;
    void dpscore(double tolerance, vector<vector<long>> &allRetIdx, int threadnum, vector<vector<float>>& accDist, ICQuery &query ,vector<vector<int>> &dpscores) override;



private:
    void toMzScanFilePair(ofstream &fout_mzFile, DataFile &df, bool rmParentIon, 
    SpectraFilter *sf, bool verbose);                  //CScanInfo

    void initPeakNum(bool verbose);
    
    //CScanInfo
    void loadScanMzToMemory(bool isList, bool verbose);                     //CScanInfo

    
    void createWithList(bool rmParentIon, const vector<string> &filelist);     //CScanInfo
    void init(bool redo, bool is_file_list, bool rmParentIon);
    int getPeakNumPerSpec() const override;
    void loadMzFile(bool verbose);


};

class CMzFileReaderFactory: public ICMzFactory
{
public:
    shared_ptr<ICMzFile> create() override;
};



#endif //COMPACTSPECMZ_H
