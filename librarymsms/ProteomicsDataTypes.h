//
// Created by wulong on 4/19/18.
//

#ifndef MYTOOL_PROTEOMICSDATATYPES_H
#define MYTOOL_PROTEOMICSDATATYPES_H

static const double y_ion_prefix_mass = 18.0106;

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <set>

using namespace std;
class PeakList;
class BinningPeakList;
class SPsmAnnotation;
class PeptideProphetParser;



// todo: should not be member function
bool comparePeptide_I_equls_L(string &pepA, string &pepB);

void L2Normalization(vector<float> &v, int dim);
void L2Normalization(float *v, int dim);
float *get_float_vector(BinningPeakList *bpl, int nLen, bool verbose=false);
void getFloatVecPaddingZeros(BinningPeakList *bpl, int nLen, bool verbose, vector<float> &v);
float * vpl_to_normalized_vec(vector<PeakList *> &vpl, int dim, bool useFlankingBins , const int topPeakNum);
void Preprocess(PeakList *pkl, const int PeakNum);
vector<string> readlines(const string &mzXMLList);

//
//def setModiCompound():
//mydic={'_Methyl_K_14.015': [1,2,0,0,0,0], 'Nitro_Y': [0,-1,1,2,0,0], '_Carbamidomethyl_C_57.021': [2,3,1,1,0,0],
//'tri-Methylation_K': [3,6,0,0,0,0], 'di-Methylation_K': [2,4,0,0,0,0], 'Methyl_R': [1,2,0,0,0,0],
//'Sulfation_Y': [0,0,0,3,1,0], 'di-Methylation_R': [2,4,0,0,0,0], 'Acetyl_K': [2,2,0,1,0,0],
//'Deamidated[Q]':[0,-1,-1,1,0,0],'Carbamidomethyl[C]':[2,3,1,1,0,0],'Gly->Pro[G]':[3,4,0,0,0,0],'Acetyl[AnyN-term]':[2,2,0,1,0,0],
//'Acetyl[ProteinN-term]':[2,2,0,1,0,0],'Carbamyl[AnyN-term]':[1,1,1,1,0,0],'Gln->pyro-Glu[AnyN-termQ]':[0,-3,-1,0,0,0],
//'Deamidated[N]':[0,-1,-1,1,0,0],'Oxidation[M]':[0,0,0,1,0,0],'Methyl[K]': [1,2,0,0,0,0], 'Nitro[Y]': [0,-1,1,2,0,0], 'Carbamidomethyl[C]': [2,3,1,1,0,0],
//'Trimethyl[K]': [3,6,0,0,0,0], 'Dimethyl[R]': [2,4,0,0,0,0], 'Methyl[R]': [1,2,0,0,0,0],
//'Sulfo[Y]': [0,0,0,3,1,0], 'Dimethyl[K]': [2,4,0,0,0,0], 'Acetyl[K]': [2,2,0,1,0,0]}
//return mydic
//
//
//# Sulfo[Y]   1.0   0.0010775862069
//# Nitro[Y]   78.0   0.0840517241379
//# Acetyl[K]   73.0   0.0786637931034
//# Methyl[K]   111.0   0.119612068966
//# Dimethyl[K]   87.0   0.09375
//# Trimethyl[K]   78.0   0.0840517241379
//# Methyl[R]   47.0   0.0506465517241
//# Dimethyl[R]   56.0   0.0603448275862
//
//        def setAminoAcid():
//        dic={}


// a module for any amino acid.
class CAminoAcid {
    double m_residue_mass_Da;
    string m_residue_name;
public:
    CAminoAcid(string name, double mass);
    ~CAminoAcid();
    double getMass() const;
    string getName();
};

class CAAConsts {
    static CAAConsts *pConst;
    CAAConsts();
public:
	~CAAConsts();
    vector<CAminoAcid *> AAs;
    map<string, CAminoAcid *> name2AA;

    static CAAConsts *get();

private:
	void Create20AAs(vector<CAminoAcid *> &AAs);
	void CreateAAMap(vector<CAminoAcid *> &AAs, map<string, CAminoAcid *> &nameToAA);
};

// this is a struct
class CPeakAnnotation {
public:
    string m_ion_type;
    int m_charge_status;
    int m_postion;
    string m_neutral_loss;
    string m_other;
public:
    CPeakAnnotation();
    ~CPeakAnnotation();
    void print() const;
};

class CPeak {
public:
    double m_mass;
    double m_inten;
    // one peak could have multiple annotations
    vector<CPeakAnnotation *> m_anno;
public:
    CPeak(double mass, double intensity);
    CPeak(const CPeak & other){
//        cout << "coloning peaks" << endl;
        m_mass = other.m_mass;
        m_inten = other.m_inten;
        for(auto each :  other.m_anno){
//            cout << "coloning annotations " << each << endl;
            m_anno.push_back(new CPeakAnnotation(*each));
        }
    }
    ~CPeak();  // safe release
    void AddAnnotation(CPeakAnnotation *anno);
    void print();
};

class TheoSpec {
    vector<CPeak *> m_ions;
    void add_ions(const string& iontype, double ion_mass, vector<CAminoAcid *> AAs, int from_charge, int to_charge);
    void sortPeaks();

public:
    TheoSpec(vector<CAminoAcid *> AAs, int max_charge_allowed);
    ~TheoSpec(); // safe release
    void printSpec();
    CPeak *getCPeak(int i);

    int getPeakNum();

};

class CPeptide {
    string m_seq;
    vector<CAminoAcid *> m_AAs;
    TheoSpec *m_theo_spec;
    int m_precursor_charge;
public:
    CPeptide(string peptide_seq, map<string, CAminoAcid *> allAAs, int charge);
	CPeptide(string peptide_seq, int charge);
    ~CPeptide();
    string getPepSeq();
    TheoSpec *generateTheoreticalSpectra();
};

// todo: for annotaiton, I could do something more: parent ion and immonium ions there. (also internal ions)
class DataFile;

class CSpectrum {
public:
    vector<CPeak *> m_PeakList;
    int m_precursor_charge;
    string other;
    string m_spectrum_name; // mgf title
    DataFile *m_data_source;
    double m_parentMz;
    CPeptide *m_pep;
    double m_rt_in_sec; // RT
    int m_MS_Level; // 1 or 2
    int m_scanNum;
    // protein info
    string m_protein;
    string m_collision_energy;

public:
    CSpectrum();

    ~CSpectrum();
    CSpectrum(const CSpectrum &other_spec);
    double getPrecursorNeutralMass() const;
    string getProtein() const;
    string getCollisionEnergy() const;
    void setProtein(string protein);
    void setCollisionEnery(string collisionEnergy);
    int getScanNum() const;
    void setScanNum(int scanNum);
    void setMSLevel(int ms_level);
    void setRTinSeconds(double RT_In_Sec);
    void setSpectrumName(string spectrum_name);
    void setParentMz(double parentMz);
    void setParentCharge(int parentCharge);
    void addOnePeak(double mz, double intensity);
    int getPeakNum() const;
    int getMSLevel() const;
    string getSpectrumName() const;
    void getOnePeak(double &mz, double &intensity, int peak_idx);

    void print();
    void annotate(string pepstr, int charge);
    void annotate(CPeptide *pep);

    // lower bound; upper bound
    int find_max_peak_within(double center,double halfWidthTol);

    vector<double> keepTopNinRange(int N, double width);
    vector<bool> keepTopNinRangeBool(int N, double width, bool verbose);

    vector<double> deIsotope(bool rmIsotopicPeaks=false);
    void getAllPeaks(vector<double> &mz, vector<double> &intensity, bool removeLowIntensePeaks, bool rmParentIon,
                     bool rmIsotopicPeaks, double localMaxTolWidth);
    bool isLocalMaximum(int i, double halfWidthTol);
    double maxIntensity(bool rmParentIon, double left, double right);
    vector<double> top2Intensity(bool rmParentIon, double left, double right);
};

class ICFileParser {
public:
    virtual ~ICFileParser();

    // this function is so difficult to write.
    // Shall we load everything to memory?
    // Of course not!
    // We do not have that much momery in any time.
    // So, we should get a handle for that file or spectrum
    // We could access the spectrum anytime very fast.
    // We should maintain a cache to keep some of the spectra
    // how about this! Each time we keep one file in memory?
    // Great!
    virtual void load(string filename, vector<CSpectrum *> &spec) =0;
    virtual void load(string filename, vector<CSpectrum *> &spec, int start, int end)=0;
};

shared_ptr<ICFileParser> FileParserFactory(string extension);

class SptxtParser: public ICFileParser{
public:
    SptxtParser();
    ~SptxtParser();
    void load(string filename, vector<CSpectrum*> &spec) override;
    void load(string filename, vector<CSpectrum*> &spec, int start, int end) override;
};

class MGFParser : public ICFileParser {
public:
    MGFParser();
    ~MGFParser();
    void load(string filename, vector<CSpectrum *> &spec) override;
    void load(string filename, vector<CSpectrum*> &spec, int start, int end) override;

};

class mzXMLParser : public ICFileParser {
public:
    mzXMLParser();

    ~mzXMLParser();

    void load(string filename, vector<CSpectrum *> &spec) override;
    void load(string filename, vector<CSpectrum*> &spec, int start, int end) override;
};

class mzMLParser : public ICFileParser {
public:
    mzMLParser();
    ~mzMLParser();
    void load(string filename, vector<CSpectrum *> &spec);
    void load(string filename, vector<CSpectrum*> &spec, int start, int end) override;
};

class scanListParser: public ICFileParser{
    shared_ptr<DataFile> m_df;
public:
    scanListParser(){
        m_df=nullptr;
    }
    ~scanListParser(){}

    void load(string filename, vector<CSpectrum *> &spec);
    void load(string filename, vector<CSpectrum *> &spec, int start, int end) override;
};

class Chromatogram{
    struct chromatogramPeak{
        double mz;
        double rt;
        double intensity;
        chromatogramPeak(double _mz, double _rt, double _intensity):mz(_mz), rt(_rt), intensity(_intensity){}
    };
    vector<chromatogramPeak> m_peaks;
    double m_mzTol;
public:
    double m_precursorMz;
    double m_ms2Scan;
    double m_ms2_rt_in_second;
    int m_precursor_charge;
    Chromatogram(vector<CSpectrum*> & spectra, double mzTol, int ms2Scan, double mzOffset);
    void print(double minRT, double maxRT, int starWeight);
    void toVector(vector<vector<double>> &x, double rt_tol);
    string toJsonStr(double rt_tol);
};

class DataFile {
    vector<CSpectrum *> m_Spectra;
    string m_sourcefile;
    map<int, int> m_ms2idx_to_idx;
    map<int, int> m_Scan_to_idx;
public:
    DataFile(const string& filename, int start, int end);
    ~DataFile();
    string getSourceFileName();
    CSpectrum *getSpectrum(long i);
    CSpectrum * getSpectrumByScan(int scan);
    int getIdxByScan(int scan);
    void print(int specnum);
    vector<PeakList *> toPeakList(double localMaxHalfWidth, int msLevel, bool remove_precursor, long start_spec_id=0, long end_spec_id=-1);

    float *toFloatVector(int dim, long &specnum, bool removeprecursor, bool useFlankingBins,
            double tolerance, long start_spec_id=0, long end_spec_id=-1, const int topPeakNum=50);
    void printSummary();
    int getSpectrumNum();
    int getMS2SpectrumNum();

    int getIdx(int ms2idx);

    shared_ptr<Chromatogram> extractChromatogram(double MzTol, int MS2scan, double mzOffset, int charge)
    {
        return make_shared<Chromatogram>(m_Spectra, MzTol, MS2scan, mzOffset);
    }
};


// Never used!
class DataArchive {
    vector<DataFile *> m_DataSource;
public:
    DataArchive();
    ~DataArchive();
    void addFile(const string& filename);
    void print();
    DataFile *getDataFile(int index);
    int getDataFileNum();
};


//
//return dic
//        def makeuplost():
//        dic={}
//dic['CO']=27.9949
//dic['H2O']=18.0106
//dic['NH3']=17.0265
//dic['H2ONH3']=35.0371
//return dic
//
//        def setAAcompound():
//        dic={}
//dic['A']=[3,5,1,1,0,0]
//dic['R']=[6,12,4,1,0,0]
//dic['N']=[4,6,2,2,0,0]
//dic['D']=[4,5,1,3,0,0]
//dic['C']=[3,5,1,1,1,0]
//dic['E']=[5,7,1,3,0,0]
//dic['Q']=[5,8,2,2,0,0]
//dic['G']=[2,3,1,1,0,0]
//dic['H']=[6,7,3,1,0,0]
//dic['I']=[6,11,1,1,0,0]
//dic['L']=[6,11,1,1,0,0]
//dic['K']=[6,12,2,1,0,0]
//dic['M']=[5,9,1,1,1,0]
//dic['F']=[9,9,1,1,0,0]
//dic['P']=[5,7,1,1,0,0]
//dic['S']=[3,5,1,2,0,0]
//dic['T']=[4,7,1,2,0,0]
//dic['W']=[11,10,2,1,0,0]
//dic['Y']=[9,9,1,2,0,0]
//dic['V']=[5,9,1,1,0,0]
//dic['X']=[1,1,1,1,1,0]
//dic['B']=[3,5,1,1,1,0]
//dic['J']=[5,9,1,1,1,0]
//dic['Z']=[1,1,1,1,1,0]
//dic['X']=[1,1,1,1,1,0]
//dic['U']=[1,1,1,1,1,0]
//dic['O']=[1,1,1,1,1,0]
//return dic
//
//        def setIsotopeDistribution():
//        dic={}
//dic['C']=[0.988930,0.011070,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
//dic['H']=[0.99985,0.00015,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
//dic['N']=[0.996337,0.003663,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
//dic['O']=[0.997590,0.000374,0.002036,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
//dic['S']=[0.9502,0.0075,0.0421,0.0002,0.0,0.0,0.0,0.0,0.0,0.0]
//return dic



class CMzSpec {
    uint16_t  *m_p;
    const int m_PeakNum;
public:
    CMzSpec(uint16_t *p, int PeakNum);
    CMzSpec(const CMzSpec & other);
    uint16_t *getPeakPtr();
    ~CMzSpec();
    string getpeakmsg();
    void display();
    void printMZStr();
    int getPeakNum();
    uint16_t getMz(int intensity);
    uint16_t getMzOfRank(int rank);
};

// id: 0 ~ 50
// PeakPerSpec = 50
class ICId2RankInten {
public:
    virtual int getPeakPerSpec() =0;
    virtual int toIntensity(int id) =0;
    virtual int getMaxScore()=0;
    float getSquaredNorm(const uint16_t *p);
    float getNorm(const uint16_t *p);
};
class CId2RankInten: public ICId2RankInten{
    const int PeakPerSpec;
    int maxScore;
public:
    CId2RankInten(int peak_per_spec, bool verbose);
    int getMaxScore() override;
    int getPeakPerSpec() override;
    int toIntensity(int id) override;
};


#endif //MYTOOL_PROTEOMICSDATATYPES_H
