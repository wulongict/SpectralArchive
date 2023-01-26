//
// Created by wulong on 4/19/18.
//

#include "ProteomicsDataTypes.h"
#include "Util.h"
#include "SpectraST_cramp.hpp"
#include "XMLFileParser.h"
#include "PeakList.h"
#include <spdlog/spdlog.h>

#include <utility>
using SpectraST_msms::cRamp;
using SpectraST_msms::rampRunInfo;
using SpectraST_msms::rampScanInfo;
using SpectraST_msms::rampPeakList;

bool comparePeptide_I_equls_L(string &pepA, string &pepB) {
    if(pepA.size()!=pepB.size()) return false;
    bool ret = true;
    for(int i = 0; i < pepA.size(); i ++)    {
        if(pepA[i] == pepB[i]){continue;}
        else if((pepA[i]=='I' and pepB[i] == 'L')   or (pepA[i]=='L' and pepB[i] == 'I') )     {
            continue;
        }  else{
            ret = false;
            break;
        }
    }
    return ret;
}

void L2Normalization(vector<float> &v, int dim) {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += (v[i] * v[i]);
    }
    if (sum > EPSILON) {
        sum = sqrt(sum);
        for (int i = 0; i < dim; i++) {
            v[i] /= sum;
        }
    }
}


void L2Normalization(float *v, int dim) {
    double sum = accumulate(v, v + dim, 0.0L, [](const double &s, const double &x) { return s + x * x; });
    if (sum > EPSILON) {
        sum = sqrt(sum);
        transform(v, v + dim, v, [&](double x) { return x / sum; });
    }
}

void getFloatVecPaddingZeros(BinningPeakList *bpl, int nLen, bool verbose, vector<float> &v) {
    int N = bpl->GetBinNum();
    if (N > nLen) {
        cout << "vector real length: " << N << " is larger than space : " << nLen << endl;
        throw out_of_range("vector size is shorter than real length");
    }
    if (verbose) cout << "getting memory" << endl;
    v.assign(nLen, 0);
    int nonzeroLen = bpl->nonzeros.size();
    for (int i = 0; i < nonzeroLen; i++) {
        int j = bpl->nonzeros[i];
        v[j] = bpl->GetIntensity(j);
    }
}

/// Get vector as a float array! This function is not very efficient! will be removed!
float *get_float_vector(BinningPeakList *bpl, int nLen, bool verbose) {
    int N = bpl->GetBinNum();
    if (N > nLen) {
        cout << "vector real length: " << N << " is larger than space : " << nLen << endl;
        throw out_of_range("vector size is shorter than real length");
    }
    if (verbose) cout << "getting memory" << endl;

    float *v = new float[nLen];

    if (v == nullptr) {
        cerr << v << " empty v" << endl;
        throw runtime_error("out of memory!");
    }

    if (verbose) {
        cout << "now memory v is " << v << endl;

    }
    if (verbose) cout << "Start to get intensity" << endl;
    for (int i = 0; i < N; i++) {
        v[i] = bpl->GetIntensity(i);
    }
    if (verbose) cout << "padding zeros at the end of the vector" << endl;
    for (int i = N; i < nLen; i++) {
        v[i] = 0;
    }
    return v;
}

float *vpl_to_normalized_vec(vector<PeakList *> &vpl, int dim, bool useFlankingBins,const int topPeakNum) {
    long long num_spectra_current_batch = vpl.size();

    float *results = new float[num_spectra_current_batch * dim];

    if (results == nullptr) {
        cout << "fail to get memory" << endl;
        throw runtime_error("out of memory!");
    }
    try {
        //Progress ps(vpl.size(), "Vectorize");
        // Attention: bug on overflow of index. i.
        for (long long i = 0; i < vpl.size(); i++) {
            //ps.increase();
            if (vpl[i] == nullptr) {
                for (int k = 0; k < dim; k++) {
                    results[k + i * dim] = 0;
                }
            } else {
                Preprocess(vpl[i],topPeakNum);
                BinningPeakList *bpl = vpl[i]->CreateBinningPeaks(useFlankingBins);
                if (bpl == nullptr) throw runtime_error("out of memory: bpl");
                float *vec = get_float_vector(bpl, dim, false);
                if (vec == nullptr) throw runtime_error("out of memory: vec");
                L2Normalization(vec, dim);
                for (int k = 0; k < dim; k++) {
                    results[k + i * dim] = vec[k];
                }
                delete[] vec;
            }
        }
    }
    catch (exception &ex) {
        cout << "unexpected error: " << ex.what() << endl;
    }
    catch (...) {
        cout << "Unknow error: " << __FILE__ << "\t" << __LINE__ << endl;
    }
    return results;
}

void Preprocess(PeakList *pkl,const int PeakNum) {
    pkl->KeepTopN(PeakNum);
    pkl->rankingAsIntensity(PeakNum);
}

vector<string> readlines(const string &mzXMLList) {
    vector<string> filelist;
    cout << "reading from file " << mzXMLList << endl;
    ifstream fin(mzXMLList.c_str(), ios::in);
    int fileid = 0;
    if (fin.is_open()) {
        string line;
        while (getline(fin, line)) {
            line.erase(remove(line.begin(), line.end(),'\r'), line.end());
            filelist.push_back(line);
        }
        fin.close();
    } else {
        cout << "[Error] List file does not exist! filename: \"" << mzXMLList << "\""<< endl;
        //throw logic_error("File does not exist!");
    }
    return filelist;
}

CAminoAcid::CAminoAcid(string name, double mass) {
    this->m_residue_mass_Da = mass;
    m_residue_name = std::move(name);
}

CAminoAcid::~CAminoAcid() = default;

double CAminoAcid::getMass() const {
    return m_residue_mass_Da;
}

string CAminoAcid::getName() {
    return m_residue_name;
}

CPeakAnnotation::CPeakAnnotation() {}

CPeakAnnotation::~CPeakAnnotation() = default;

void CPeakAnnotation::print() const {
    cout << m_ion_type << m_postion << "^" << m_charge_status << "NL:" << m_neutral_loss << "\t" << m_other << " ";
}

CPeak::CPeak(double mass, double intensity) {
    m_mass = mass;
    m_inten = intensity;
}

CPeak::~CPeak() {
    // Release the annotaitons
    for (auto x : m_anno) {
        delete x;
    }
}

void CPeak::AddAnnotation(CPeakAnnotation *anno) {
    m_anno.push_back(anno);
}

void CPeak::print() {
    cout << m_mass << "\t" << m_inten << "\t ";
    for (auto x: m_anno) x->print();
    cout << endl;
}

// TheoSpec can only be created from a CPeptide Object
TheoSpec::TheoSpec(vector<CAminoAcid *> AAs, int max_charge_allowed) {
    add_ions("b", 0, AAs, 1, max_charge_allowed);

    reverse(AAs.begin(), AAs.end());
    add_ions("y", y_ion_prefix_mass, AAs, 1, max_charge_allowed);

    sortPeaks();
}

TheoSpec::~TheoSpec() {
    for (auto x: m_ions) {
        delete x;
    }
}

void TheoSpec::add_ions(const string& iontype, double ion_mass, vector<CAminoAcid *> AAs, int from_charge, int to_charge) {
    const double proton_mass = 1.007276; // to be verified
    for (int i = 0; i < AAs.size() - 1; i++) {
        CAminoAcid *x = AAs[i];
        ion_mass += x->getMass();

        for (int chg = from_charge; chg <= to_charge; chg++) {
            CPeak *p = new CPeak(0, 0);
            p->m_mass = (ion_mass + chg * proton_mass) / chg;
            p->m_inten = 100; // default value

            CPeakAnnotation *peakAnno = new CPeakAnnotation();
            peakAnno->m_charge_status = chg;
            peakAnno->m_ion_type = iontype;
            peakAnno->m_neutral_loss = "0";
            peakAnno->m_postion = i + 1;
            peakAnno->m_other = "";

            p->AddAnnotation(peakAnno);

            m_ions.push_back(p);
        }
    }
}

void TheoSpec::printSpec() {
    for (auto x: m_ions) {
        x->print();
    }
}

void TheoSpec::sortPeaks() {
    sort(m_ions.begin(), m_ions.end(), [](const CPeak *x, const CPeak *y) -> bool { return x->m_mass < y->m_mass; });
}

CPeak *TheoSpec::getCPeak(int i) {
    return m_ions[i];
}

int TheoSpec::getPeakNum() {
    return m_ions.size();
}

// todo: do not copy allAAs, just use a pointer or reference
CPeptide::CPeptide(string peptide_seq, map<string, CAminoAcid *> allAAs, int charge) {
    m_precursor_charge = charge;

    m_seq = std::move(peptide_seq);
    m_theo_spec = nullptr;
    for (int i = 0; i < m_seq.size(); i++) {
        string name = m_seq.substr(i, 1);
        m_AAs.push_back(allAAs[name]);
    }
    cout << "Peptide: " << m_seq << " Length: " << m_seq.length() << " Charge: " << charge << endl;
}


CPeptide::CPeptide(string peptide_seq, int charge){
	m_precursor_charge = charge;
	m_seq = std::move(peptide_seq);
	m_theo_spec = nullptr;
	map<string, CAminoAcid*> &allAAs = CAAConsts::get()->name2AA;
	for (int i = 0; i < m_seq.size(); i++) {
		string name = m_seq.substr(i, 1);
		m_AAs.push_back(allAAs[name]);
	}
	cout << "Peptide: " << m_seq << " Length: " << m_seq.length() << " Charge: " << charge << endl;
}

CPeptide::~CPeptide() = default;

string CPeptide::getPepSeq() { return m_seq; }

TheoSpec *CPeptide::generateTheoreticalSpectra() {
    if (nullptr == m_theo_spec) {
        m_theo_spec = new TheoSpec(m_AAs, m_precursor_charge);
    }
    if (nullptr == m_theo_spec) {
        cout << "[Info] Fail to generate theoretical spectrum for peptide " << m_seq << endl;
    }
    return m_theo_spec;
}


shared_ptr<ICFileParser> FileParserFactory(string extension) {
    //string extension_lowercase;
    transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    //extension_lowercase = extension;
    if (extension == "mgf") {
        return make_shared<MGFParser>();//MGFParser();
    } else if (extension == "mzml") {
        return make_shared<mzMLParser>();//new mzMLParser();  //mz
    } else if (extension == "mzxml") {
        return make_shared<mzXMLParser>();//new mzXMLParser(); // m
    } else if (extension == "sptxt") {
        return make_shared<SptxtParser>();//new SptxtParser();  //mz
    }else if (extension == "scanlist") {
        return make_shared<scanListParser>();//new SptxtParser();  //mz
    } else {
        cout << "Unknown file format: " << extension << endl;
        return nullptr;
    }
}

// Warning: using start, end could break the indexing system. ms2idx, idx
DataFile::DataFile(const string& filename, int start, int end) {
    if(end < start)    {
        //cout << "Warning: invalid range [" << start << ", " << end << "] changed to [0, " <<INT_MAX << "]" << endl;
        end=INT_MAX;
        start=0;
    }

    string ext = File::CFile(filename).ext;
    shared_ptr<ICFileParser> fp = FileParserFactory(ext);

    // todo: we forget to release the pointer!
    if (fp) {
        //cout << "[Info] loading file: " << filename << endl;
//        cout << "start loading " << start << " to " << end<< endl;
         fp->load(filename, m_Spectra,start,end);
        m_sourcefile = filename;
    } else {
        cout << "[Error] File format error! ext = " << ext << endl;
    }
    //printSummary();
    for(int i = 0; i < getSpectrumNum(); i ++){

        int scan = getSpectrum(i)->getScanNum();
        m_Scan_to_idx[scan]=i;
    }
}

DataFile::~DataFile() {
    //cout << "Releasing spectrum " << m_Spectra.size() << endl;
    for (auto x: m_Spectra) {
            delete x;
    }
}

string DataFile::getSourceFileName() { return m_sourcefile; }

CSpectrum *DataFile::getSpectrum(long i) {
    if (i >= 0 and i < m_Spectra.size()) {
        return m_Spectra.at(i);
    } else {
        cout << "[Error] DataFile::getSpectrum() Invalid index " << i << endl;
        return nullptr;
    }
}

void DataFile::print(int specnum) {
    cout << "//Spectrum " << specnum << "  begins " << endl;
    if (specnum >= 0 and specnum < m_Spectra.size()) {
        if (m_Spectra[specnum] != nullptr) {
            m_Spectra[specnum]->print();
        }
    }
    cout << "//Spectrum " << specnum << " ends" << endl;
}


float *DataFile::toFloatVector(int dim, long &specnum, bool removeprecursor, bool useFlankingBins,
        double tolerance, long start_spec_id, long end_spec_id, const int topPeakNum) {
    if(end_spec_id<start_spec_id)    {
        end_spec_id = getSpectrumNum();
    }
    //cout << "Mz file contains peaks with lossy compression..." << endl;
    int msLevel = 2;
    vector<PeakList *> vpl = toPeakList(tolerance, msLevel, removeprecursor, start_spec_id, end_spec_id);
    specnum = vpl.size();

    float *ret = vpl_to_normalized_vec(vpl, dim, useFlankingBins,topPeakNum);
    releaseVectorPeakListPtr(vpl);
    return ret;
}

void DataArchive::addFile(const string& filename) {
    DataFile *p = new DataFile(filename,0,-1);
    if (nullptr != p) {
        m_DataSource.push_back(p);
    } else {
        cout << "Error: file format not supported! " << filename << endl;
    }
}

void DataArchive::print() {
    cout << "There are " << m_DataSource.size() << " files in this DataArchive!" << endl;
    for (int i = 0; i < m_DataSource.size(); i++) {
        if (m_DataSource[i] != nullptr) {
            cout << "File " << i << "\t" << m_DataSource[i]->getSourceFileName() << endl;
        }
    }
}

DataArchive::DataArchive() = default;

DataArchive::~DataArchive() {
    // add release of DataArchive
    cout << "Releasing space for datafile in Data Archive" << endl;
    for (auto & i : m_DataSource) {
        if (i != nullptr) {
            delete i;
            i = nullptr;
        }
    }
}

DataFile *DataArchive::getDataFile(int index) {
    if (index >= 0 and index < m_DataSource.size()) {
        return m_DataSource[index];
    } else {
        return nullptr;
    }
}

int DataArchive::getDataFileNum() {return m_DataSource.size();}

CSpectrum::~CSpectrum() {
    delete m_pep;
    for (auto x: m_PeakList) {

            delete x;
    }
}

void CSpectrum::setSpectrumName(string spectrum_name) {
    m_spectrum_name = spectrum_name;
}

void CSpectrum::setParentMz(double parentMz) {
    m_parentMz = parentMz;
}

void CSpectrum::setParentCharge(int parentCharge) {
    m_precursor_charge = parentCharge;
}

void CSpectrum::addOnePeak(double mz, double intensity) {
    CPeak *p = new CPeak(mz, intensity);
    if (p != nullptr) {
        m_PeakList.push_back(p);
    } else {
        cout << "[Info] fail to create peak " << endl;
    }
}

void CSpectrum::print() {
    cout << "SpecName\t" << m_spectrum_name << endl;
    cout << "ScanNum\t" << m_scanNum << endl;
    cout << "msLevel\t" << m_MS_Level << endl;
    cout << "retention time\t" << m_rt_in_sec << endl;
    cout << "HCD\t" << m_collision_energy << endl;
    cout << "Protein\t" << m_protein << endl;
    if (m_MS_Level == 2) {
        cout << "ParentMZ\t" << m_parentMz << endl;
        cout << "ParentChg\t" << m_precursor_charge << endl;
    }

    cout << "-- Peak Num --" << getPeakNum() << endl;
    for (auto x: m_PeakList) {
        if (x != nullptr)
            x->print();
    }
}

void CSpectrum::annotate(CPeptide *pep) {
    m_pep = pep;
    double half_tol = 0.5;
    TheoSpec *tSpec = pep->generateTheoreticalSpectra();
    tSpec->printSpec();
    int ePeakNum = m_PeakList.size();
    int tPeakNum = tSpec->getPeakNum();
    int i = 0, j = 0;
    while (i < ePeakNum and j < tPeakNum) {
        CPeak *ePeak = m_PeakList[i], *tPeak = tSpec->getCPeak(j);
        if (ePeak->m_mass - tPeak->m_mass > half_tol) {
            // experimental peak larger,
            j++;
        } else if (tPeak->m_mass - ePeak->m_mass > half_tol) {
            // theoretical peak is larger
            i++;
        } else {
            // the same anno; shared pointer
            ePeak->m_anno = tPeak->m_anno;
            i++;
            // todo: to be continued
        }
    }
}

void CSpectrum::annotate(string pepstr, int charge) {
    CPeptide *pep = new CPeptide(pepstr, charge);
    annotate(pep);
}

int CSpectrum::getPeakNum() const {
    return m_PeakList.size();
}

bool CSpectrum::isLocalMaximum(int i, double halfWidthTol) {
    int j =find_max_peak_within(m_PeakList[i]->m_mass,halfWidthTol);
    return j == i;
}

vector<double> CSpectrum::top2Intensity(bool rmParentIon, double left, double right) {
    vector<double> top2Intensity{0,0};
//    if(m_PeakList.size() == 0){
//        top2Intensity[0] = m_PeakList[0]->m_inten;
//    }
    for(int i = 0; i < m_PeakList.size(); i ++){
        auto &pk = m_PeakList[i];
        if(rmParentIon and
           (m_parentMz+left<pk->m_mass and pk->m_mass < m_parentMz + right
           or m_parentMz*2+left<pk->m_mass and pk->m_mass < m_parentMz*2 + right))        {
            continue;
        }
        if(pk->m_inten > top2Intensity[0]){
            // larger than first one.
            top2Intensity[1] = top2Intensity[0];
            top2Intensity[0] = pk->m_inten;
        } else if ( pk->m_inten > top2Intensity[1]){
            top2Intensity[1] = pk->m_inten;
        }
    }

    return top2Intensity;

}

// retrun peak with maximum intensity, skipping peaks related to precursor mz. 
double CSpectrum::maxIntensity(bool rmParentIon, double left, double right) {
    double maxinten = 0;
    for(auto & i : m_PeakList)    {
        if(rmParentIon and
        (m_parentMz+left<i->m_mass and i->m_mass < m_parentMz + right
        or m_parentMz*2+left<i->m_mass and i->m_mass < m_parentMz*2 + right))        {
            continue;
        } else if(maxinten< i->m_inten) {
            maxinten = i->m_inten;
        }
    }
    return maxinten;
}

// There are many preprocessing steps.
void CSpectrum::getAllPeaks(vector<double> &mz, vector<double> &intensity, bool removeLowIntensePeaks, bool rmParentIon,
                            bool rmIsotopicPeaks, double localMaxTolWidth) {

    int peaknum = m_PeakList.size();
    mz.assign(peaknum, 0);
    intensity.assign(peaknum, 0);
    int k = 0;
    double left = -17, right = 3;
    double inten_threshold = 0.02 * maxIntensity(rmParentIon, left, right);
    vector<double> top2Inten = top2Intensity(rmParentIon, left, right);
//    cout << "threshold from " << 0.02 * top2Inten[0] << " to " << 0.02*top2Inten[1] << endl;
    inten_threshold = 0.02 * top2Inten[1];
    // change to second most Intense peak.

    vector<double> isIsotopicPeak = deIsotope(rmIsotopicPeaks);
//    vector<double> isTopInRange = keepTopNinRange(6,100);
    vector<bool> isTopInRangeBool = keepTopNinRangeBool(6, 100, false);
    const double MAX_PEAK_MZ = 2000.0;
    vector<string> whyRemoved(peaknum,"PASS");
    bool verbose = false;
    for(int i = 0; i < peaknum; i ++)    {
        if(removeLowIntensePeaks and m_PeakList[i]->m_inten < inten_threshold)    {
            if(verbose) whyRemoved[i] = to_string("less than threshold: "," ", inten_threshold);
            continue;
        }
        else if (rmParentIon and m_parentMz+left<m_PeakList[i]->m_mass and m_PeakList[i]->m_mass < m_parentMz + right) {
            if(verbose) whyRemoved[i] = "is near precursor mz";
            continue;
        }
        else if (rmParentIon and m_parentMz*2-proton_mass+left<m_PeakList[i]->m_mass and m_PeakList[i]->m_mass < m_parentMz*2-proton_mass + right) {
            // assuming it is a charge state 2 precursor, remove the corresponding ion.
            if(verbose) whyRemoved[i] = "is near precursor mz*2";
            continue;
        }
        else if(rmIsotopicPeaks and isIsotopicPeak[i] > 0.9){
            //cout << "deleted: " << m_PeakList[i]->m_mass << endl;
            if(verbose) whyRemoved[i] = "is Isotopic peak";
            continue;
        }
        else if(not isLocalMaximum(i, localMaxTolWidth) )   {
            if(verbose) whyRemoved[i] = "not local maximal";
            continue;
        }
        else if(m_PeakList[i]->m_mass>MAX_PEAK_MZ){ // skip peak with mz greater than 2000
            if(verbose) whyRemoved[i] = "out of range [0,2000]";
            continue;
        }

        else if(not isTopInRangeBool[i])     {
            if(verbose) whyRemoved[i] = "not top K in range";
            continue;
        }
        getOnePeak(mz[k], intensity[k], i);
        k++;
    }
    mz.resize(k);
    intensity.resize(k);
    if(verbose){
        for(int i = 0; i < m_PeakList.size(); i ++){
            cout << m_PeakList[i]->m_mass << "\t"
            << m_PeakList[i]->m_inten << "\t" << whyRemoved[i] << endl;
        }
    }
}

vector<double> CSpectrum::deIsotope(bool rmIsotopicPeaks) {
    int peaknum = m_PeakList.size();
    vector<double> isIsotopicPeak(peaknum, 0);
    if(rmIsotopicPeaks)    {
        int max_charge = 7;
        double average_gap = 1.0032;

        for(int i = 0; i < peaknum; i ++)        {
            double monoisotope = m_PeakList[i]->m_mass;
            double halfWidthTol = monoisotope * 20/1e6; // 20ppm
            if(isIsotopicPeak[i]>0.9) {
                continue;
            }
            for(int chg = 1; chg <= max_charge; chg ++)    {
                vector<int> isoid={i};
                double mass_diff = 1.0032/chg;
                int max_num_isotopic_peaks = 4;
                for(int k = 1; k <max_num_isotopic_peaks;  k ++)    {
                    double newpeak = monoisotope + mass_diff *k;
                    int id = find_max_peak_within(newpeak,halfWidthTol);
                    isoid.push_back(id);
                }
                if(isoid[1]>0)                {
                    int monoid = isoid[0];
                    for(int k = 1; k < isoid.size(); k ++) {
                        int x = isoid[k];
                        if(x == -1) break;
                        bool isIso = m_PeakList[monoid]->m_inten > m_PeakList[x]->m_inten;
                        if(isIso) isIsotopicPeak[x] = 1;
                    }
                }
            }
        }
    }
    return isIsotopicPeak;
}
#include <queue>
vector<bool> CSpectrum::keepTopNinRangeBool(int N, double width, bool verbose) {
    // only keep top N peaks in a window of [a, a+ width]
    int peaknum = getPeakNum();
    vector<bool> isTopN(peaknum, true);
    if(N<=0 or width<=0 or peaknum <=1) {
        if(verbose)cout << "topN " << N << " window width: " << width << " number of peaks to be processed: " << peaknum << endl;
        // invalid input
//        cout << "Invalid parameter for retaining top N peaks: " << N << " " << width << endl;
        return isTopN;
    }

    // lambda expression
    auto cmp = [&](int a, int b){
            return m_PeakList[a]->m_inten>m_PeakList[b]->m_inten;
    };
    std::priority_queue<int, std::vector<int>, decltype(cmp)> peakInWin(cmp), tmpQ(cmp);
    auto ptrQ = &peakInWin, ptrTmp = &tmpQ;

    for(int i = 0; i < peaknum; i ++){
        double mz = m_PeakList[i]->m_mass;
        ptrQ->push(i);
        if( ptrQ->size()>N){
//            printf("\npeakInWinSize: %d; topPeak: %d; peakAdded: %d\n", ptrQ->size(), ptrQ->top(), i);
            while(not ptrQ->empty()){
                int k = ptrQ->top();
                if(m_PeakList[k]->m_mass > mz-width) {ptrTmp->push(k); }
                ptrQ->pop();
            }
            std::swap(ptrQ, ptrTmp);
            if(ptrQ->size()>N){
                isTopN[ptrQ->top()]=false;
                ptrQ->pop();
            }
        }
    }
    return isTopN;
}


vector<double> CSpectrum::keepTopNinRange(int N, double width) {
    // only keep top N peaks in a windows of [a , a + width]
    int peaknum = m_PeakList.size();
    vector<double> isTopN(peaknum, 1);

    for(int i = 0; i < peaknum -N; i ++)    {
        if(isTopN[i]<0.5) continue;
        vector<int> v;
        for(int j = i; j < peaknum and m_PeakList[j]->m_mass - m_PeakList[i]->m_mass < width; j ++)        {
            if(isTopN[j]<0.5) continue;
            v.push_back(j);
        }
        if(v.size() < N) continue;
        sort(v.begin(), v.end(),[&](const int &a, const int &b)->bool{return m_PeakList[a]->m_inten < m_PeakList[b]->m_inten;});
        v.resize(v.size() -N);
        for(int j : v)        {
            isTopN[j] = 0;
        }
    }
    return isTopN;
}


int CSpectrum::find_max_peak_within(double center, double halfWidthTol) {

    CPeak left(center-halfWidthTol,1), right(center+halfWidthTol,1);
    // left open right open interval: (left, right)
    auto it_left = upper_bound(m_PeakList.begin(), m_PeakList.end(),&left, [](const CPeak *a , const CPeak *b){return a->m_mass < b->m_mass;});
    auto it_right = upper_bound(m_PeakList.begin(), m_PeakList.end(),&right, [](const CPeak *a , const CPeak *b){return a->m_mass < b->m_mass;});

    if(std::distance(it_left, it_right)>=1)
    {
        auto it = max_element(it_left,it_right,[](const CPeak *a, const CPeak *b){return a->m_inten<b->m_inten;});
        return std::distance(m_PeakList.begin(), it);
    }    else    {
        return -1;
    }

}

void CSpectrum::getOnePeak(double &mz, double &intensity, int peak_idx) {
    mz = m_PeakList[peak_idx]->m_mass;
    intensity = m_PeakList[peak_idx]->m_inten;
}

CSpectrum::CSpectrum() {
    m_collision_energy="";
    m_data_source = nullptr;
    m_MS_Level = -1;
    m_parentMz = 0.0;
    m_pep = nullptr;
    m_precursor_charge = 0;
    m_protein="";
    m_rt_in_sec = -1;
    m_scanNum = -1;
    m_spectrum_name="";
}

double CSpectrum::getPrecursorNeutralMass() const {
    if(m_MS_Level == 2) {
        return m_precursor_charge * m_parentMz - m_precursor_charge * proton_mass;
    }
    else {
        return -1;
    }
}

string CSpectrum::getProtein() const {return m_protein;}

string CSpectrum::getCollisionEnergy() const {return m_collision_energy;}

void CSpectrum::setProtein(string protein) {m_protein = std::move(protein);}

void CSpectrum::setCollisionEnery(string collisionEnergy) {m_collision_energy = std::move(collisionEnergy);}

int CSpectrum::getScanNum() const {return m_scanNum;}

void CSpectrum::setScanNum(int scanNum) {m_scanNum = scanNum;}

void CSpectrum::setMSLevel(int ms_level) {m_MS_Level = ms_level;}

void CSpectrum::setRTinSeconds(double RT_In_Sec) {m_rt_in_sec = RT_In_Sec;}

int CSpectrum::getMSLevel() const {return m_MS_Level;}

string CSpectrum::getSpectrumName() const {
    return m_spectrum_name;
}

CSpectrum::CSpectrum(const CSpectrum &other_spec) {
//    cout << "cloning new peak list from another one" << endl;
    for (auto  each: other_spec.m_PeakList){
        if(each==nullptr){
//            cout << "empty peak!" << endl;
            m_PeakList.push_back(nullptr);
        }else{
            m_PeakList.push_back(new CPeak(*each));
        }
    }
//    cout << "peak list copied: " << m_PeakList.size() << " <-" << other_spec.m_PeakList.size() << endl;
    m_precursor_charge = other_spec.m_precursor_charge;
    other = other_spec.other;
    m_spectrum_name = other_spec.m_spectrum_name;
    m_data_source = other_spec.m_data_source;
    m_parentMz = other_spec.m_parentMz;
    if(other_spec.m_pep != nullptr){
        m_pep = new CPeptide(*other_spec.m_pep);

    }else{
        m_pep = nullptr;
    }
    m_rt_in_sec = other_spec.m_rt_in_sec;
    m_MS_Level = other_spec.m_MS_Level;
    m_scanNum = other_spec.m_scanNum;
    m_protein = other_spec.m_protein;
    m_collision_energy = other_spec.m_collision_energy;

}


mzMLParser::mzMLParser() : ICFileParser() {}

void mzMLParser::load(string filename, vector<CSpectrum *> &spec) {
    mzMLReader mzML(filename);
//        mzML.get_all_spectra(); // loading
    vector<Spectrum> *all = mzML.getSpectra();
    for (auto & i : *all) {
        CSpectrum *onespec = new CSpectrum();
        onespec->setRTinSeconds(i.m_RT_in_seconds);
        onespec->setMSLevel(i.m_ms_level);
        onespec->setParentMz(i.m_precursormz);
        onespec->setParentCharge(i.m_precursorcharge);
        onespec->setScanNum(i.spectrum_scan);
        onespec->setSpectrumName(i.spectrum_name);
        vector<double> mzs = i.pkl->getM_mzList();
        vector<double> intens = i.pkl->getM_intensityList();
        for (int j = 0; j < mzs.size(); j++) {
            onespec->addOnePeak(mzs[j], intens[j]);
        }
        spec.push_back(onespec);
    }
}

mzMLParser::~mzMLParser() = default;

void mzMLParser::load(string filename, vector<CSpectrum *> &spec, int start, int end) {
    cout << "Attention: parameter start and end (" << start << " " << end << ") are not used! All spectra will be loaded" << endl;
    load(filename,spec);
}

void MGFParser::load(string filename, vector<CSpectrum *> &spec) {
    bool warningIndexAsScanOnce = false;
    SimpleTimer st("MGF Reader");
    string line;
    if(not File::isExist(filename,true)){
        throw runtime_error(filename + ": file not exist");

    }
    ifstream fin(filename, ios_base::in);
    if(not fin.good()){
        cout << "Fail to open file!!!" << endl;
//        cout << "File opened successfully " << endl;
    }
    CSpectrum *pSpec = nullptr;
    struct mgfSummary{
        int num_comment_lines;
        int unexpected_lines;
        mgfSummary(){num_comment_lines=0; unexpected_lines=0;}
        void print() const{
            cout << "--mgf summary--" << endl;
            cout << "comment lines: " << num_comment_lines << endl;
            cout << "unexpected lines: " << unexpected_lines << endl;
        }
    } mgf_summary;

    CountProgress cp(10000, "Parsing MGF");
    while (fin.good() and getline(fin, line)) {
        if(line.empty()) {
            continue; // fix the bug of an empty line.
        }
        int startpos = line.find_first_not_of(' ');
        int endpos = line.find_last_not_of(' ');
//        cout << "start pos and end pos " << startpos << "\t" << endpos << endl;
        line = line.substr(startpos, endpos-startpos + 1);
        if(not line.empty() and line[0] == '#'){
            mgf_summary.num_comment_lines++;
        }
        else if (string::npos != line.find("BEGIN IONS")) {
            pSpec = new CSpectrum();
            cp.increase();
            continue;
        } else if (string::npos != line.find("PEPMASS")) {
            double pepmass = stod(line.substr(line.find('=') + 1));
            pSpec->setParentMz(pepmass);
        } else if (string::npos != line.find("TITLE")) {
            string title = line.substr(line.find('=') + 1);
            pSpec->setSpectrumName(title);
            vector<string> tokens;
            split_string(title, tokens, '.');
            if (tokens.size() >= 4 && string::npos == title.find_first_of('/'))    {
                pSpec->setScanNum(atoi(tokens[1].c_str())); // start scan
            } else {
                if (not warningIndexAsScanOnce) {
                    cout << "Warning: Making scan number with index +1! " << endl;
                    warningIndexAsScanOnce = true;
                }
                pSpec->setScanNum(spec.size() + 1);
            }
        } else if (string::npos != line.find("CHARGE")) {
            string charge = line.substr(line.find('=') + 1);
            pSpec->setParentCharge(atoi(charge.c_str()));
        }  else if (string::npos != line.find("RTINSECONDS")) {
            string rtStr = line.substr(line.find('=') + 1);
            pSpec->setRTinSeconds(stod(rtStr.c_str()));
        } else if (string::npos != line.find("END IONS")) {
            if (pSpec == nullptr) {
                cout << "Error:" << "Spectra without  (BEGIN IONS)" << endl;
                exit(-1);
            }
            pSpec->setMSLevel(2);
            spec.push_back(pSpec);
            pSpec = nullptr;
            continue;
        } else if(not line.empty() and isdigit(line[0])) {
            string::size_type space_pos;
            double mz = stod(line, &space_pos);
            line.erase(0, space_pos + 1);
            double intensity = stod(line, &space_pos);
            pSpec->addOnePeak(mz, intensity);
        } else{
            mgf_summary.unexpected_lines++;
        }
    }
    fin.close();
}

MGFParser::MGFParser() : ICFileParser() {
}

MGFParser::~MGFParser() = default;

// TODO: start and end not used properly.
void MGFParser::load(string filename, vector<CSpectrum *> &spec, int start, int end) {
    cout << "[Warning] MGF reader do not support selecting spectra range, will load all" << endl;
//    cout << "Do not use this function: " << __FUNCTION__ << " as it is not efficient " << start << " " << end << " not defined!" << endl;
    load(filename,spec);
}

CAAConsts *CAAConsts::pConst = nullptr;

CAAConsts::CAAConsts() {
    Create20AAs(AAs);
    CreateAAMap(AAs, name2AA);
}

CAAConsts::~CAAConsts(){
	for (auto x : AAs)	{
		if (x)		{
			delete x;
			x = nullptr;
		}
	}
}

CAAConsts *CAAConsts::get() {
    if (nullptr == CAAConsts::pConst) {
        CAAConsts::pConst = new CAAConsts();
    }
    return pConst;
}


void CAAConsts::Create20AAs(vector<CAminoAcid *> &AAs){
	AAs.push_back(new CAminoAcid("A", 71.03711));
	AAs.push_back(new CAminoAcid("R", 156.10111));
	AAs.push_back(new CAminoAcid("N", 114.04293));
	AAs.push_back(new CAminoAcid("D", 115.02964));
	AAs.push_back(new CAminoAcid("C", 103.00919));
	AAs.push_back(new CAminoAcid("E", 129.04259));
	AAs.push_back(new CAminoAcid("Q", 128.05858));
	AAs.push_back(new CAminoAcid("G", 57.02146));
	AAs.push_back(new CAminoAcid("H", 137.05891));
	AAs.push_back(new CAminoAcid("I", 113.08406));
	AAs.push_back(new CAminoAcid("L", 113.08406));
	AAs.push_back(new CAminoAcid("K", 128.09496));
	AAs.push_back(new CAminoAcid("M", 131.04049));
	AAs.push_back(new CAminoAcid("F", 147.06841));
	AAs.push_back(new CAminoAcid("P", 97.05276));
	AAs.push_back(new CAminoAcid("S", 87.03203));
	AAs.push_back(new CAminoAcid("T", 101.04768));
	AAs.push_back(new CAminoAcid("W", 186.07931));
	AAs.push_back(new CAminoAcid("Y", 163.06333));
	AAs.push_back(new CAminoAcid("V", 99.06841));
	AAs.push_back(new CAminoAcid("X", 101.0));
	AAs.push_back(new CAminoAcid("B", 103.00919));
	AAs.push_back(new CAminoAcid("J", 131.04049));
	AAs.push_back(new CAminoAcid("Z", 0.0));
}

void CAAConsts::CreateAAMap(vector<CAminoAcid *> &AAs, map<string, CAminoAcid *> &nameToAA){
	for (auto & AA : AAs) {
		nameToAA[AA->getName()] = AA;
	}
}

void mzXMLParser::load(string filename, vector<CSpectrum *> &spec) {
    cRamp *cramp = new cRamp(filename.c_str());
    if (!cramp->OK()) {
        cout << "[Error] Cannot open file \"" << filename << "\". File skipped." << endl;
        delete (cramp);
        throw runtime_error("File does not exist!");
    }

    rampRunInfo *runInfo = cramp->getRunInfo();

    if (!runInfo) {
        cout << "Cannot open file \"" << filename << "\". File skipped." << endl;
        exit(0);
    }
    int emptySpecNum = 0;
    //Progress ps(cramp->getLastScan(), filename);
    for (int k = 1; k <= cramp->getLastScan(); k++) {
        //ps.increase();
        rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);
        CSpectrum *p = new CSpectrum();

        if (!scanInfo) {
            cout << "----- No scan Info----: Memory Leaking !!! " << endl;
            continue;
        }
        if (!p) {
            cout << "Could not create new spectrum!" << endl;
            throw runtime_error("Not enough memory in heap!");
        }
        p->setRTinSeconds(scanInfo->getRetentionTimeSeconds());

        int scan_num = scanInfo->m_data.acquisitionNum;
        p->setSpectrumName(File::CFile(filename).basename + "." + to_string(scan_num) + "." + to_string(scan_num));
        p->setMSLevel(scanInfo->m_data.msLevel);
        p->setParentMz(scanInfo->m_data.precursorMZ);
        p->setParentCharge(scanInfo->m_data.precursorCharge);
        p->setScanNum(scanInfo->m_data.acquisitionNum);

        rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);
        delete scanInfo;
        if (peaks) {
            for (int j = 0; j < peaks->getPeakCount(); j++) {
                p->addOnePeak(peaks->getPeak(j)->mz, peaks->getPeak(j)->intensity);
            }
            delete peaks;
        } else {
            emptySpecNum++;
        }
        spec.push_back(p);
    }
    if(emptySpecNum>0) cout << "[Info] Empty spectra: " << emptySpecNum << endl;

    delete cramp;
}

mzXMLParser::mzXMLParser() : ICFileParser() {}

mzXMLParser::~mzXMLParser() = default;

void mzXMLParser::load(string filename, vector<CSpectrum *> &spec, int start, int end) {
    //cout << "Warning: mzXML Reader do not support selecting spectra range, will implement later." << endl;
//    cout << "Warning: do not used this function!!! to be implemented" << endl;
    load(filename,spec);
}

void SptxtParser::load(string filename, vector<CSpectrum *> &spec) {
    SimpleTimer st("sptxt Reader");
    CountProgress cp(1000, "Parsing sptxt");
    string line;
    ifstream fin(filename, ios_base::in);
    if (not fin.good()) {
        cout << "Fail to open file " << filename << endl;
        exit(0);
    }
    CSpectrum *pSpec = nullptr;
    while (getline(fin, line)) {
        // bugfix: maybe there is a carriage return in the empty line!! CR!!!
        int spacepos = line.find_last_not_of("\r\n");
        if (spacepos != string::npos or !line.empty()) {
            line.erase(spacepos + 1);
        }
        // remove zeros!
        int startpos = line.find_first_not_of(' ');
        int endpos = line.find_last_not_of(' ');

        if (startpos != string::npos) {
            line = line.substr(startpos, endpos-startpos + 1);

        } else {
            if (pSpec == nullptr) {
                cout << "Empty line and No spectrum read!" << endl;
                continue;
            } else {
                line = "";
            }
        }
//            trim_space_only(line);
        if (0 == line.find("Name: ")) { // find Name: KSLDFJSLAJFLDSJKFLSJDFLSDJKFLKSDKJFLAJWEIROUEIFJ/3
            pSpec = new CSpectrum();
            string title = line.substr(line.find(' ') + 1);
            int charge = atoi(title.substr(title.find("/") + 1).c_str());
            pSpec->setParentCharge(charge);
            pSpec->setSpectrumName(title);

            cp.increase();
            continue;

        } else if (0 == line.find("LibID: ")) { // fomd "LibID: 343"
            int libid = atoi(line.substr(line.find(':') + 1).c_str());

            pSpec->setScanNum(libid + 1);


        } else if (0 == line.find("Comment: ")){
            string commentstr = line.substr(strlen("Commet: "));
            CKeyValuesParser kvp(commentstr, "=", " ", false);
            pSpec->setCollisionEnery(kvp.getvalue("HCD"));
            pSpec->setProtein(kvp.getvalue("NISTProtein"));
        }

        else if (0 == line.find("PrecursorMZ: ")) {
            double pepmass = stod(line.substr(line.find(":") + 1).c_str());
            pSpec->setParentMz(pepmass);
        } else if (line.length() == 0) { // end of a spectra
            if (pSpec == nullptr) {
                cout << "Error:" << "Empty spectrum ? " << endl;
                exit(-1);
            }
            pSpec->setMSLevel(2);
            spec.push_back(pSpec);
            pSpec = nullptr;
            continue;
        } else if (not line.empty() and isdigit(line[0])) { // peaks number
            string::size_type space_pos;
            double mz = stod(line, &space_pos);
            line.erase(0, space_pos + 1);
            double intensity = stod(line, &space_pos);
            pSpec->addOnePeak(mz, intensity);
        } else {
            continue;
        }
    }
    fin.close();
}

// read all then release the ones from 0 to start...
// todo: this function can be improved by calling previous load function, without start end parameter first, then remove extra ones.
//  Anyway, this is not good solution. should be improved later.
void SptxtParser::load(string filename, vector<CSpectrum *> &spec, int start, int end) {
    if(end <= start)    {
        cout << "range is zero" << endl;
        return ;
    }

    SimpleTimer st("sptxt Reader");
    CountProgress cp(1000, "Parsing sptxt");
    string line;
    ifstream fin(filename, ios_base::in);
    if (not fin.good()) {
        cout << "Fail to open file " << filename << endl;
        exit(0);
    }
    CSpectrum *pSpec = nullptr;
    while (getline(fin, line)) {
        // bugfix: maybe there is a carriage return in the empty line!! CR!!!
        int spacepos = line.find_last_not_of("\r\n");
        if (spacepos != string::npos or !line.empty()) {
            line.erase(spacepos + 1);
        }

        int startpos = line.find_first_not_of(' ');
        int endpos = line.find_last_not_of(' ');

        if (startpos != string::npos) {
            line = line.substr(startpos, endpos -startpos+ 1);
        } else {
            if (pSpec == nullptr) {
                cout << "Empty line and No spectrum read!" << endl;
                continue;
            } else {
                line = "";
            }
        }

        if (0 == line.find("Name: ")) { // find Name: KSLDFJSLAJFLDSJKFLSJDFLSDJKFLKSDKJFLAJWEIROUEIFJ/3
            if(spec.size()==end){
                break;
            }
            pSpec = new CSpectrum();
            string title = line.substr(line.find(' ') + 1);
            int charge = atoi(title.substr(title.find('/') + 1).c_str());
            pSpec->setParentCharge(charge);
            pSpec->setSpectrumName(title);

            cp.increase();
            continue;

        } else if (0 == line.find("LibID: ")) { // fomd "LibID: 343"
            int libid = atoi(line.substr(line.find(':') + 1).c_str());
            pSpec->setScanNum(libid + 1);
        } else if (0 == line.find("Comment: ")){
            string commentstr = line.substr(strlen("Commet: "));

            CKeyValuesParser kvp(commentstr, "=", " ", false);

//            kvp.display();
            pSpec->setCollisionEnery(kvp.getvalue("HCD"));
            pSpec->setProtein(kvp.getvalue("NISTProtein"));
//            cout << "parsing protein" << endl;
            if (pSpec->getProtein() == ""){
                pSpec->setProtein(kvp.getvalue("Protein"));
            }
            string retentionTime_triplets =kvp.getvalue("RetentionTime");
            if(not retentionTime_triplets.empty()){
//                cout << retentionTime_triplets << "|<- what is this?" << endl;
                vector<string> threeRTs;
                split_string(retentionTime_triplets,threeRTs,',');
//                cout << "parsing RT" << endl;
                pSpec->setRTinSeconds(stringTo<double>(threeRTs.at(0)));
            }

        }
        else if (0 == line.find("PrecursorMZ: ")) {
            double pepmass = stod(line.substr(line.find(":") + 1).c_str());
            pSpec->setParentMz(pepmass);
        } else if (line.length() == 0) { // end of a spectra
            if (pSpec == nullptr) {
                cout << "Error:" << "Empty spectrum ? " << endl;
                exit(-1);
            }
            pSpec->setMSLevel(2);
            spec.push_back(pSpec);
            pSpec = nullptr;
            continue;
        } else if (line.length() > 0 && line[0] >= '0' && line[0] <= '9') { // peaks number
            string::size_type space_pos;
            double mz = stod(line, &space_pos);
            line.erase(0, space_pos + 1);
            double intensity = stod(line, &space_pos);
            pSpec->addOnePeak(mz, intensity);
        } else {
            continue;
        }
    }
    fin.close();
    // release 0 to start
    for(int i = 0; i < spec.size() and i < start; i ++)    {
        // release
        delete spec[i];
        spec[i]=nullptr;
    }
    spec.erase(spec.begin(),spec.begin()+start);

}

SptxtParser::SptxtParser() : ICFileParser(){

}

SptxtParser::~SptxtParser() = default;

/// todo: to be replaced by a new function
vector<PeakList *> DataFile::toPeakList(double localMaxHalfWidth, int msLevel, bool remove_precursor, long start_spec_id, long end_spec_id) {
    if (remove_precursor) {
        // comment out for clean output.
        // cout << "[Info] peaks related to precursor will be removed: [mz - 17, mz + 3]" << endl;
    }
    if(end_spec_id<start_spec_id )    {
        end_spec_id = getSpectrumNum();
    }

    //Progress ps(end_spec_id-start_spec_id, "PeakList");
    vector<PeakList *> vpl;
    for (long i = start_spec_id; i < end_spec_id; i++) {
        //ps.increase();
        CSpectrum *spec = getSpectrum(i);
        if (spec == nullptr or spec->getMSLevel() != msLevel) {
            continue;
        }

        vector<double> mz, intensity;
        bool removeLowIntensePeaks = true;
        bool rmIsotopicPeaks = true;
        spec->getAllPeaks(mz, intensity, removeLowIntensePeaks, remove_precursor, rmIsotopicPeaks, localMaxHalfWidth);

        PeakList *pl = new PeakList();
        if (pl == nullptr) {
            cout << "fail to get space from server" << endl;
            throw runtime_error("out of memory!");
        }
        pl->setM_mzList(mz);
        pl->setM_intensityList(intensity);

        vpl.push_back(pl);
    }

    return vpl;
}

void DataFile::printSummary() {
    int ms1counts = 0;
    int ms2counts = 0;
    int empty_spec_num = 0;
    for (auto x: m_Spectra) {
        if (x == nullptr) {
            empty_spec_num++;
            exit(0);
        }
        if (x->getMSLevel() == 1) ms1counts++;
        else if (x->getMSLevel() == 2) ms2counts++;
    }
    cout << endl << "Summary: #Spec:\t" << m_Spectra.size() << "\tMS1:\t" << ms1counts << "\tMS2:\t" << ms2counts
         << endl << endl;
}

int DataFile::getSpectrumNum() {return m_Spectra.size();}

int DataFile::getMS2SpectrumNum() {
    int ms2count = 0;
    for(int i = 0; i < m_Spectra.size();i ++)
    {
        if(m_Spectra[i] == nullptr){cout << "Empty spectrum: " << i << endl; continue;}
        if(m_Spectra[i]->getMSLevel() == 2)
        {
            ms2count ++;
        }
    }
    return ms2count;
}

int DataFile::getIdx(int ms2idx) {
    // build map : ms2idx --> idx
    if(m_ms2idx_to_idx.count(ms2idx) == 0)
    {
        int ms2counts = 0;
        for (int k = 0; k < this->getSpectrumNum(); k++) {
            if (this->getSpectrum(k)->getMSLevel() != 2) continue;
            m_ms2idx_to_idx[ms2counts] = k;
            ms2counts++;
        }
    }
    return m_ms2idx_to_idx[ms2idx];
}

CSpectrum *DataFile::getSpectrumByScan(int scan) {

    return getSpectrum(getIdxByScan(scan));
}

int DataFile::getIdxByScan(int scan) {
    if(m_Scan_to_idx.count(scan)==0){
        cout << "Error: Scan " << scan << " not found in current data file" << endl;
    }
    return m_Scan_to_idx[scan];
}


float ICId2RankInten::getSquaredNorm(const uint16_t *p) {
    float norm = 0;
    for (int i = 0; i < getPeakPerSpec(); i++) {
        if (p[i] == 0) {
            break;
        } else {
            norm += toIntensity(i) * toIntensity(i);
        }
    }
    return norm;
}

float ICId2RankInten::getNorm(const uint16_t *p) {
    return sqrt(getSquaredNorm(p));
}

CId2RankInten::CId2RankInten(int peak_per_spec, bool verbose) : PeakPerSpec(peak_per_spec){
    maxScore = 0;
    for(int i = 0; i < PeakPerSpec; i ++)    {
        maxScore += (i+1)*(i+1);
    }
    if(verbose)cout <<  "[Info] MAX score is (42925 if peaknum=50): " << maxScore <<" and peaknum=" << PeakPerSpec<< endl;
}

int CId2RankInten::getMaxScore() {return maxScore;}

int CId2RankInten::getPeakPerSpec() {return PeakPerSpec;}

int CId2RankInten::toIntensity(int id) {
    return getPeakPerSpec() - id;
}

string CMzSpec::getpeakmsg() {
    const double MAX_SPEC_MZ = 2000;
    vector<double> mz(m_PeakNum, 0);
    string pkmsg;
    for (int i = 0; i < m_PeakNum; i++) {
        mz[i] = m_p[i] * 1.0 / UINT16_MAX * MAX_SPEC_MZ;
        pkmsg += to_string(mz[i]);
        pkmsg += " ";
    }
    return pkmsg;
}

void CMzSpec::display() {
    const double MAX_SPEC_MZ = 2000;
    vector<double> mz(m_PeakNum, 0), intensity(m_PeakNum, 0);
    for (int i = 0; i < m_PeakNum; i++) {
        mz[i] = m_p[i] * 1.0 / UINT16_MAX * MAX_SPEC_MZ;
        intensity[i] = m_PeakNum - i;
    }
    vector<int> idx(m_PeakNum, 0);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&mz](const int &x, const int &y) { return mz[x] < mz[y]; });
    cout << "id\tmz\tintensity\n-------------------------" << endl;
    for (int i = 0; i < m_PeakNum; i++) {
        cout << i << "\t" << mz[idx[i]] << "\t" << intensity[idx[i]] << endl;
    }
    cout << "--------------------------" << endl;
}

CMzSpec::CMzSpec(uint16_t *p, int PeakNum) :m_PeakNum(PeakNum){
    m_p = p;
}

CMzSpec::CMzSpec(const CMzSpec &other) :m_PeakNum(other.m_PeakNum){
    m_p = other.m_p;
}

uint16_t *CMzSpec::getPeakPtr() {
    return m_p;
}

CMzSpec::~CMzSpec() = default;

void CMzSpec::printMZStr() {
    cout << "----------------------START------------" << endl;
    if(m_PeakNum > 0) cout << m_p[0];
    for(int i = 1; i < m_PeakNum; i ++)
    {
        cout <<"," <<  m_p[i];
    }
    cout << endl << "------------------END------------" << endl;
}

int CMzSpec::getPeakNum() {
    int num = 0;
    for(int i = 0; i < m_PeakNum; i ++)   {
        if(m_p[i] == 0)   {
            break;
        }
        num ++;
    }
    return num;
}

uint16_t CMzSpec::getMz(int intensity) {
    if(intensity >=1 and intensity <=m_PeakNum) {
        return m_p[m_PeakNum - intensity];
    } else{
        cout << "Error: invalid intensity" << endl;
        exit(0);
    }
}

uint16_t CMzSpec::getMzOfRank(int rank) {
    if(rank>=0 and rank < m_PeakNum)    {
        return m_p[rank];
    } else{
        cout << "[Error] Incorrect rank " << rank << endl;
        throw logic_error("invalid rank");
    }
}

Chromatogram::Chromatogram(vector<CSpectrum *> &spectra,  double mzTol, int ms2Scan, double mzOffset) {
//    m_precursorMz = precursorMz;
    m_mzTol = mzTol;
    m_ms2Scan = ms2Scan;
    // find all ms1 index
    vector<int> ms1Idx;
    int precedingMS1 = -1;
    for(int i = 0; i < spectra.size(); i ++)
    {
        int mslevel = spectra[i]->getMSLevel();
        int scanNum = spectra[i]->getScanNum();
        if(mslevel == 1){
            ms1Idx.push_back(i);
            if(scanNum<m_ms2Scan){
                precedingMS1 = scanNum;

            }
        }else if(scanNum == ms2Scan){
            m_ms2_rt_in_second = spectra[i]->m_rt_in_sec;
            m_precursorMz = spectra[i]->m_parentMz;
            m_precursor_charge = spectra[i]->m_precursor_charge;

        }
    }
    // find ms1 before ms2,
    if(precedingMS1 == -1){
        cout << "no MS1 spectra found!" << endl;
    } else{
        double searchMz = m_precursorMz+mzOffset;
        for (int i : ms1Idx) {
            int idx = spectra[i]->find_max_peak_within(searchMz,m_mzTol);
            chromatogramPeak chrPeak(0,0,0);
            if(idx == -1)
            {
                chrPeak.rt = spectra[i]->m_rt_in_sec;
                chrPeak.mz = searchMz;
                if(chrPeak.rt<1379.5 and chrPeak.rt > 1379)
                cout << "no peak found: RT=" << chrPeak.rt << " mz=" << chrPeak.mz << "mzOffset=" << mzOffset << " searchMz=" << searchMz << endl;
                // no peak in range
            } else{
                spectra[i]->getOnePeak(chrPeak.mz, chrPeak.intensity, idx);
                chrPeak.rt = spectra[i]->m_rt_in_sec;
            }
            m_peaks.push_back(chrPeak);
        }

    }
}

void Chromatogram::print(double minRT, double maxRT, int starWeight) {
    cout << "chromatogram" << endl;
    for(auto & m_peak : m_peaks)
    {
        double rt = m_peak.rt ;
        if(rt < maxRT and rt > minRT){
        int stars = m_peak.intensity/ starWeight;
        cout << m_peak.rt << "\t" << m_peak.intensity << "\t" << string(stars,'*') << endl;

        }
    }
}

void Chromatogram::toVector(vector<vector<double>> &x, double rt_tol) {
    double minRT = m_ms2_rt_in_second - rt_tol, maxRT = m_ms2_rt_in_second + rt_tol;
    for(auto & m_peak : m_peaks)
    {
        double rt = m_peak.rt;
        if(rt < maxRT and rt > minRT)
        {
            x.emplace_back(vector<double>{m_peak.rt,m_peak.intensity, m_peak.mz});
        }
    }
}

string Chromatogram::toJsonStr(double rt_tol) {
    ostringstream oss;
    oss << "{" << endl;
    oss << "\"precursorMz\": " << m_precursorMz << endl;
    oss << ",\"ms2scan\": " << m_ms2Scan << endl;
    oss << ",\"ms2RT\": " << m_ms2_rt_in_second << endl;
    oss << ",\"mzTol\": " << m_mzTol << endl;
    oss << ",\"RtTol\": " << rt_tol << endl;
    string mz_str, rt_str, intensity_str;
    double minRT = m_ms2_rt_in_second - rt_tol;
    double maxRT = m_ms2_rt_in_second + rt_tol;
    bool first_value = true;
    for(auto & m_peak : m_peaks){
        if(m_peak.rt > maxRT or m_peak.rt < minRT) continue;
        if(not first_value){
            mz_str += " ";
            intensity_str += " ";
            rt_str += " ";
        }
        mz_str += to_string(m_peak.mz);
        intensity_str += to_string(m_peak.intensity);
        rt_str += to_string(m_peak.rt);
        first_value = false;
    }
    oss << R"(,"mzs": ")" << mz_str << "\"" << endl;
    oss << R"(,"rt": ")" << rt_str << "\"" << endl;
    oss << R"(,"intensity": ")" << intensity_str << "\"" << endl;
    oss << "}" << endl;
    return oss.str();

}

void scanListParser::load(string filename, vector<CSpectrum *> &spec) {
    // get the filename and the scan number list
    ifstream fin(filename, ios::in);
    string datafilename;
    int scan;
    set<int> scanList;
    fin >> datafilename;
    while(not fin.eof()){
        fin >> scan;
        scanList.insert(scan);
    }
    // now we have a list of id stored.
    m_df = make_shared<DataFile> (datafilename, 0, -1);

    for(auto it=scanList.begin(); it!=scanList.end(); it ++){
        CSpectrum * cspec = m_df->getSpectrumByScan(*it);
        CSpectrum * p = new CSpectrum(*cspec);
        spec.push_back(p);
    }
}

void scanListParser::load(string filename, vector<CSpectrum *> &spec, int start, int end) {
    load(filename, spec);
}
