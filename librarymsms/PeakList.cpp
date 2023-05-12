#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <numeric>
#include <set>
#include "Util.h"
#include "PeakList.h"

using namespace std;

BinningPeakList::BinningPeakList() {
    peakcountOutOfRange = 0;
    m_intensityList.clear();
}

BinningPeakList::BinningPeakList(const vector<double> &x, const vector<double> &y, bool useFlankingBins, bool verbose) {
    peakcountOutOfRange = 0;
    const double BinSize = 0.5;
    const double MaxMZ = 2001;

    int BinNum = floor(MaxMZ / BinSize) + 1;

    m_intensityList.assign(BinNum, 0);

//    int peakCountsOutOfReange = 0;

    for (int i = 0; i < x.size(); i++) {
        int idx = floor(x[i] / BinSize);
        if (idx >= BinNum) {
            peakcountOutOfRange ++;
            idx = BinNum - 1;
        }
        m_intensityList[idx] += y[i];

        if (idx - 1 >= 0 && useFlankingBins) {
            m_intensityList[idx - 1] += y[i] / 2;
        }

        if (idx + 1 < BinNum && useFlankingBins) {
            m_intensityList[idx + 1] += y[i] / 2;
        }
    }
//    bool verbose = true;
    if(verbose and peakcountOutOfRange>0){
        cout << "Warning: peaks with m/z > 2000: " << peakcountOutOfRange << endl;
    }

    norm = -1;

    for (int j = 0; j < BinNum; ++j) {
        if (fabs(m_intensityList[j]) > EPSILON) nonzeros.push_back(j);
    }

}

BinningPeakList::~BinningPeakList() = default;

void BinningPeakList::Print() {
    for (int i = 0; i < GetBinNum(); i++) {
        cout << m_intensityList[i] << " ";
    }
    cout << endl;
}

int BinningPeakList::GetBinNum() {
    return m_intensityList.size();
}

double BinningPeakList::GetIntensity(int k) const {
    return m_intensityList.at(k);
}

double BinningPeakList::CalcDotProductByNonZeros(BinningPeakList &other) {
    double dot = 0;
    vector<int> v;
    int m = 0, n = 0;

    while (m < nonzeros.size() && n < other.nonzeros.size()) {
        if (nonzeros[m] == other.nonzeros[n]) {
            v.push_back(nonzeros[m]);
            m++;
            n++;
        } else if (nonzeros[m] < other.nonzeros[n]) m++;
        else n++;
    }

    for (int i : v) {
        double x = GetIntensity(i);
        double y = other.GetIntensity(i);
        // Richard: Tried Harmonic Mean
        // fixed by lwu
        dot += (x * y);// /(x+y);
    }
    double normx = GetNorm();
    double normy = other.GetNorm();
    //cout << normy << " " << normx << " " << dot << endl;
    double ret = 0;
    if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
        ret = dot / (normx * normy);
        cout << "normx, normy and dot" << normx << " " << normy << " " << dot << endl;
    }
    //cout << "dot product by nonzeros " << dot << endl;
    return ret;
}

// TODO:now just copy dotproduct one.
// this function seems never called?
double BinningPeakList::CalcSimilarByPnorm(BinningPeakList &other, double p) {
    cout << "--------------------We should not be here----------------" << endl;
    double dot = 0;
    double normx_p = this->CalcNorm();
    double normy_p = other.CalcNorm();
    if (normx_p < EPSILON || normy_p < EPSILON) {
        dot = 1;
    } else {
        for (int i = 0; i < m_intensityList.size(); i++) {
            double x = GetIntensity(i);
            double y = other.GetIntensity(i);
            if (x == 0) {
                dot += y / normy_p;
            } else if (y == 0)
                dot += fabs(other.GetIntensity(i) / normy_p - GetIntensity(i) / normx_p);
        }
    }

    double normx = GetNorm();
    double normy = other.GetNorm();

    double ret = 0;
    if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
        ret = dot / (normx * normy);
    }
    return ret;

}

// the norm here is normalized by power of 1/p
double BinningPeakList::CalcNorm() {
    double distance = 0;
    for (int i = 0; i < GetBinNum(); i++) {
        distance += (GetIntensity(i) * GetIntensity(i));
    }
    distance = sqrt(distance);
    return distance;
}

double BinningPeakList::GetNorm() {
    if (norm == -1) return CalcNorm();
    else return norm;
}


PeakList::PeakList() {
    m_binningPeakList = nullptr;
    m_RTinSeconds = 0;
}

PeakList::~PeakList() {
    delete m_binningPeakList;
}

void PeakList::InsertPeak(double mz, double intensity) {
    m_mzList.push_back(mz);
    m_intensityList.push_back(intensity);

}

BinningPeakList *PeakList::CreateBinningPeaks(bool useFlankingBins, bool verbose) {
    if (m_binningPeakList == nullptr) {
        m_binningPeakList = new BinningPeakList(m_mzList, m_intensityList, useFlankingBins, verbose);
        if(m_binningPeakList == nullptr){
            cerr << "Fail to create binning PeakList " << endl;
        }
    }
    return m_binningPeakList;
}


double PeakList::CalcDotProduct(PeakList &other) {
    BinningPeakList *x = CreateBinningPeaks();
    BinningPeakList *y = other.CreateBinningPeaks();
    double dot = x->CalcDotProductByNonZeros(*y);
    return dot;
}

void TestPeakList() {
    cout << "Test" << endl;
    PeakList x;
    x.InsertPeak(1, 2);
    x.InsertPeak(2, 4);
    x.InsertPeak(3, 6);
    PeakList y;
    y.InsertPeak(2, 4);
    PeakList z = y;

    cout << "should be (2.4) " << z.getM_mzList()[0] << " " << z.getM_intensityList()[0] << endl;

    cout << "DotProduct" << x.CalcDotProduct(y) << endl;
    cout << "Test Done" << endl;
}

PeakList PeakList::getPeaksWithin(double leftmz, double rightmz) {
    PeakList pl;

    if (!m_mzList.empty()) {
        std::vector<double>::iterator lower, upper;
        lower = std::lower_bound(m_mzList.begin(), m_mzList.end(), leftmz);
        upper = std::lower_bound(m_mzList.begin(), m_mzList.end(), rightmz);
        int start = lower - m_mzList.begin();
        int end = upper - m_mzList.begin();

        for (int i = start; i < end; ++i) {
            pl.InsertPeak(m_mzList[i], m_intensityList[i]);
        }
    }

    return pl;
}

PeakList::PeakList(const PeakList &pl) {
    bool useflankingbins = true; // added May 29 2020
    m_mzList = pl.getM_mzList();
    m_intensityList = pl.getM_intensityList();
    if (pl.getM_binningPeakList() == nullptr) {
        m_binningPeakList = nullptr;
    } else {
        m_binningPeakList = new BinningPeakList(m_mzList, m_intensityList, useflankingbins);
    }
    m_RTinSeconds = pl.getRTinSeconds();
}


void PeakList::addPeakList(PeakList &pl, int ms2_number) {
    vector<double> temp_mz = pl.getM_mzList();

    double mz_shift = 2000 * ms2_number;
    if (fabs(mz_shift) > EPSILON) {
        cout << "[Info] mzhshift=" << mz_shift << endl;
        //#pragma omp parallel for  schedule(dynamic)
        for (double & i : temp_mz) {
            i += mz_shift;
        }
    }

    vector<double> temp_intensity = pl.getM_intensityList();

    m_mzList.insert(m_mzList.end(), temp_mz.begin(), temp_mz.end());
    m_intensityList.insert(m_intensityList.end(), temp_intensity.begin(), temp_intensity.end());
}

double PeakList::getRTinSeconds() const {
    return m_RTinSeconds;
}

void PeakList::setRTinSeconds(double RTinSeconds) {
    m_RTinSeconds = RTinSeconds;
}

void PeakList::setM_mzList(const vector<double> &mzList) {
    m_mzList = mzList;
}

void PeakList::setM_intensityList(const vector<double> &intensityList) {
    m_intensityList = intensityList;
}

BinningPeakList *PeakList::getM_binningPeakList() const {
    return m_binningPeakList;
}

vector<double> PeakList::getM_mzList() const {
    return m_mzList;
}

vector<double> PeakList::getM_intensityList() const {
    return m_intensityList;
}

void PeakList::print() {
    cout << "// Peak Start " << endl;
    cout << "// i\tmz\tinten" << endl;
    for(int i = 0; i < m_mzList.size(); i ++)    {
        cout << i << "\t" << m_mzList[i] << "\t" << m_intensityList[i] << endl;
    }
    cout << "// Peak End" << endl;
}

 void PeakList::KeepTopN(int N) {
    if (N >= m_intensityList.size()) {
        return; // nothing to do
    }

    // Sort mzList and intensityList together based on intensity
    std::vector<std::pair<double, double>> mz_intensity_pairs(m_intensityList.size());
    for (std::size_t i = 0; i < m_intensityList.size(); ++i) {
        mz_intensity_pairs[i] = std::make_pair(m_mzList[i], m_intensityList[i]);
    }
    std::partial_sort(mz_intensity_pairs.begin(), mz_intensity_pairs.begin() + N, mz_intensity_pairs.end(),
                      [](const std::pair<double, double>& lhs, const std::pair<double, double>& rhs) { return lhs.second > rhs.second; });

    // Copy the top N values back into mzList and intensityList
    for (std::size_t i = 0; i < N; ++i) {
        m_mzList[i] = mz_intensity_pairs[i].first;
        m_intensityList[i] = mz_intensity_pairs[i].second;
    }

    // Resize the lists to keep only the top N elements
    m_mzList.resize(N);
    m_intensityList.resize(N);
}

// void PeakList::KeepTopN(int N) {
//     double threshold = 0;
//     vector<double> v = m_intensityList;
//     sort(v.begin(), v.end(), std::greater<double>());
//     if (N <= v.size()) {
//         threshold = v[N - 1];
//     }
//     vector<double> tmp_mz, tmp_intensity;
//     for (int i = 0; i < m_mzList.size(); ++i) {
//         if (m_intensityList[i] >= threshold) {
//             tmp_mz.push_back(m_mzList[i]);
//             tmp_intensity.push_back(m_intensityList[i]);
//         }
//     }

//     m_mzList = tmp_mz;
//     m_intensityList = tmp_intensity;
// }

void PeakList::rankingAsIntensity(int maxRanking) {
    vector<int> index(m_mzList.size(),0);
    iota(index.begin(), index.end(),0);

    sort(index.begin(), index.end(), [&](const int &x, const int &y) {return m_intensityList[x] > m_intensityList[y];} );

    for(int j : index)    {
        m_intensityList[j] = maxRanking > EPSILON ? maxRanking: 0;
        maxRanking --;
    }
}

void PeakList::removePeaksWithin(double mz_left, double mz_right) {
    vector<double> tmp_mz, tmp_intensity;
    tmp_mz.reserve(m_mzList.size());
    tmp_intensity.reserve(m_intensityList.size());
    for(int i = 0; i < m_mzList.size(); i ++)    {
        if(m_mzList[i] < mz_left or m_mzList[i] > mz_right)         {
            tmp_mz.push_back(m_mzList[i]);
            tmp_intensity.push_back(m_intensityList[i]);
        }
    }

    m_mzList.swap(tmp_mz);
    m_intensityList.swap(tmp_intensity);
}

void PeakList::setBinningPeakList(BinningPeakList *binlist) {
    m_binningPeakList = binlist;
}

void PeakList::NormalizedToSum() {
    double sum = 0;
    for (double i : m_intensityList) {
        sum += i;
    }
    for (double & j : m_intensityList) {
        j /= sum;
    }
}

void releaseVectorPeakListPtr(vector<PeakList *> &vpl) {
    for (auto & i : vpl) {
        if (i != nullptr) {
            delete i;
            i = nullptr;
        }
    }
//    cout << "End of releasing" << endl;
}

VectorFilter *VectorFilterFactoryMethod(const string& FilterType) {
    if (FilterType == "PeakSquareRoot")
        return new PeakSquareRoot();
    else if (FilterType == "PeakToProb")
        return new NormalizedToProb();
    else if (FilterType == "AboveMean") {
        return new RemovePeakLessThanMean();
    } else if (FilterType == "AboveMedian") {
        return new RemovePeakLessThanMedian();
    } else if (FilterType.length() > 3 && FilterType.substr(0, 3) == "Top") {
        KeepTopNPeaks *ret = new KeepTopNPeaks();
        int N = atoi(FilterType.substr(3).c_str());
        cout << "[Info] Top-" << N << " filter is created" << endl;
        ret->setTopN(N);
        return ret;
    } else {
        cout << "[Error] Invalid Filter type: " << FilterType << endl;
        cout << "[Alert] Using base class" << endl;
        return new VectorFilter();
    }
}

PeakListFilter *CreatePeakListFilterFactoryMethod(string param) {
    if (param == "Dummy") {
        return new PeakListFilter();
    } else if (param == "NormalizeByMax") {
        return new PeakListNormalizedByMaximal();
    } else if (param.length() > 3 && param.substr(0, 3) == "Top") {
        param = param.substr(3);
        int N = atoi(param.c_str());
        PeakListKeepTopN *ret = new PeakListKeepTopN();
        ret->setTopN(N);
        return ret;
    } else {
        cout << "[Info] Invalid MS2 PeakList filter: " << param << endl;
        throw "[Error] Invalid PeakList filter for MS2";
    }
}

void BinningPeakList::setNorm(double x) {
    norm = x;
}

vector<double> BinningPeakList::GetIntensityList() const {
    return m_intensityList;
}

void BinningPeakList::setIntensitylist(vector<double> intensitylist) {
    for (double & i : intensitylist) {
        cout << i;
        m_intensityList.push_back(i);
    }
}

BinningPeakList *BinningPeakList::filterPeaks(VectorFilter *vf) {
    vf->run(this);
    vector<int> tmp;

    for (int & nonzero : nonzeros) {
        if (fabs(m_intensityList[nonzero]) > EPSILON)
            tmp.push_back(nonzero);
    }
    nonzeros.swap(tmp);

    return this;
}

void BinningPeakList::setIntensity(int i, double intensity) {
    m_intensityList[i] = intensity;
}

void PeakSquareRoot::run(BinningPeakList *x) {
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], sqrt(intensity));
    }
}

PeakSquareRoot::PeakSquareRoot() { setType("PeakSquareRoot"); }

void VectorFilter::run(BinningPeakList *x) {
    //cout << "[Info] Running base class" << endl;
}

VectorFilter::VectorFilter() { m_Type = "Dummy"; }

void VectorFilter::setType(string type) { m_Type = std::move(type); }

string VectorFilter::getType() { return m_Type; }

VectorFilter::~VectorFilter() = default;

void NormalizedToProb::run(BinningPeakList *x) {
    double sum = 0;
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        sum += intensity;
    }

    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity) / sum);
    }
}

NormalizedToProb::NormalizedToProb() { setType("PeakToProb"); }

void RemovePeakLessThanMean::run(BinningPeakList *x) {
    double mean = 0;
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        mean += intensity;
    }
    mean = mean / x->GetBinNum();

    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity - mean > 0) ? intensity - mean : 0);
    }
}

RemovePeakLessThanMean::RemovePeakLessThanMean() { setType("AboveMean"); }

void RemovePeakLessThanMedian::run(BinningPeakList *x) {
    // how to ge the median?
    if (x->nonzeros.size() < 40) {
        return;
    }

    vector<double> p;
    for (int i = 0; i < x->nonzeros.size(); ++i) {
        p.push_back(x->GetIntensity(x->nonzeros[i]));
    }
    sort(p.begin(), p.end());
    double median = 0;
    int psize = p.size();
    if (psize < 3) {
        cout << "[Info] too small peaks" << endl;
    } else if (psize % 2 == 0) {
        median = (p[psize / 2 - 1] + p[psize / 2]) / 2;
    } else {
        median = p[psize / 2];
    }
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity - median > 0) ? intensity - median : 0);
    }
}

RemovePeakLessThanMedian::RemovePeakLessThanMedian() { setType("AboveMedian"); }

void KeepTopNPeaks::run(BinningPeakList *x) {
    if (m_topN > x->nonzeros.size())
        return;
    vector<double> p;
    for (int i = 0; i < x->nonzeros.size(); ++i) {
        p.push_back(x->GetIntensity(x->nonzeros[i]));
    }
    sort(p.begin(), p.end(), std::greater<double>());

    double threshold = p[m_topN - 1];
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity > threshold) ? intensity : 0);
    }
}

KeepTopNPeaks::KeepTopNPeaks() {
    m_topN = 50;
    setType("TopN");
}

void KeepTopNPeaks::setTopN(int topN) {
    m_topN = topN;
}

const string &PeakListFilter::getM_FilterType() const {
    return m_FilterType;
}

void PeakListFilter::setM_FilterType(const string &FilterType) {
    PeakListFilter::m_FilterType = FilterType;
}

PeakListFilter::PeakListFilter() {
    m_FilterType = "Dummy";
}

PeakListFilter::~PeakListFilter() = default;

void PeakListFilter::run(PeakList *x) {}

PeakListKeepTopN::PeakListKeepTopN() {
    m_topN = 20;
    setM_FilterType("KeepTopN");
}

void PeakListKeepTopN::setTopN(int topN) {
    m_topN = topN;
}

void PeakListKeepTopN::run(PeakList *x) {
    x->KeepTopN(m_topN);
}

PeakListNormalizedByMaximal::PeakListNormalizedByMaximal() {
    setM_FilterType("NormalizeByMax");
}

void PeakListNormalizedByMaximal::run(PeakList *x) {
    x->NormalizedToSum();
}

double PCC::calcDistance(vector<double> &a, vector<double> &b) {
    double XX = 0, YY = 0, XY = 0;
    double meana = statistic::calcmean<double>(a);
    double meanb = statistic::calcmean<double>(b);
    for (int i = 0; i < a.size(); i++) {
        XX += (a[i] - meana) * (a[i] - meana);
        XY += (a[i] - meana) * (b[i] - meanb);
        YY += (b[i] - meanb) * (b[i] - meanb);
    }
    double sim = 0;

    if (fabs(XX) > EPSILON && fabs(YY) > EPSILON)
        sim = XY / sqrt(XX * YY);

    return (sim + 1) / 2;

}

double PCC::calcMean(BinningPeakList &a) {
    double ret = 0;
    for (int i = 0; i < a.nonzeros.size(); ++i) {
        ret += a.GetIntensity(a.nonzeros[i]);
    }
    return ret / a.GetBinNum();
}

double PCC::calcNorm(BinningPeakList &a) {

    double ret = 0;
    double meana = calcMean(a);
    for (int i = 0; i < a.nonzeros.size(); i++) {
        double x = a.GetIntensity(a.nonzeros[i]);
        ret += ((x - meana) * (x - meana));
    }
    ret = sqrt(ret);
    a.setNorm(ret);
    return ret;
}

PCC::PCC() = default;

double PCC::calc(BinningPeakList &a, BinningPeakList &b) {
    double dot = 0;
    double meana = calcMean(a);
    double meanb = calcMean(b);

    for (int i = 0; i < a.GetBinNum(); i++) {

        dot += ((b.GetIntensity(i) - meanb) * (a.GetIntensity(i) - meana));
    }
    double normx = a.GetNorm();
    double normy = b.GetNorm();

    double ret = 0;
    if (normy > EPSILON && normx > EPSILON) {
        ret = dot / (normx * normy);
    }

    ret = (ret + 1) / 2;

    if (ret < 0) {
        cout << "x, y ,dot " << normy << " " << normx << " " << dot << endl;
    }
    return ret;
}

double PNorm::calcDistance(vector<double> &a, vector<double> &b) {
    vector<double> distance(a.size(), 0), abs_a(a.size(), 0), abs_b(b.size(), 0);
    for (int i = 0; i < a.size(); ++i) {
        distance[i] = fabs(a[i] - b[i]);
        abs_a[i] = fabs(a[i]);
        abs_b[i] = fabs(b[i]);

    }
    double ret = 0, norm_a = 0, norm_b = 0;
    if (m_p > 0) {
        ret = statistic::calcGeneralizedMean(distance, m_p);
        // fixed: I would like to introduce the new metric by m_p, however, there are some problem here
        // How could we bound the distance,
        // should we do the normalization first?
        // We bound this by Minkowski inequality.
        // from this site, http://math.stackexchange.com/questions/581257/equality-in-minkowskis-theorem
        // we got the following conclusion.
        // The equality of pNorm(f+g) == pNorm(f) + pNorm(g) will stand when f = cg almost everywhere.
        norm_a = statistic::calcGeneralizedMean(abs_a, m_p);
        norm_b = statistic::calcGeneralizedMean(abs_b, m_p);
//            cout <<"[Info] p-norm(a-b)="  << ret << " p-norm(a)=" << norm_a << " p-norm(b)=" << norm_b << endl;
        ret = ret / (norm_a + norm_b);
    }

    else if (m_p == -1) {

        double min = *min_element(distance.begin(), distance.end());
        double max = *max_element(distance.begin(), distance.end());
        statistic::checkAndFixVector(distance, min, max);
        if (max > 0) {
            ret = statistic::calcGeneralizedMean(distance, m_p);
        }
    }
    else if (m_p == 0) {
        double min = *min_element(distance.begin(), distance.end());
        double max = *max_element(distance.begin(), distance.end());
        statistic::checkAndFixVector(distance, min, max);
        if (max > 0) {
            ret = statistic::calcGeometricMean(distance);
        }
    }
    return 1 - ret;
}

double PNorm::calc(BinningPeakList &a, BinningPeakList &b) {
    double ret = 0;

    vector<double> v;
    vector<int> nonzeros_x;
    int m = 0, n = 0;
    while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
        double x = 0, y = 0;
        if (a.nonzeros[m] < b.nonzeros[n]) {
            x = a.GetIntensity(a.nonzeros[m]);
            m++;
        }
        else if (a.nonzeros[m] > b.nonzeros[n]) {
            y = b.GetIntensity(b.nonzeros[n]);
            n++;
        }
        else {
            x = a.GetIntensity(a.nonzeros[m]);
            m++;
            y = b.GetIntensity(b.nonzeros[n]);
            n++;
        }
        ret += pow(fabs(x - y), m_p);
    }
    while (m < a.nonzeros.size()) {
        ret += pow(fabs(a.GetIntensity(a.nonzeros[m])), m_p);
        m++;
    }

    while (n < b.nonzeros.size()) {
        ret += pow(fabs(b.GetIntensity(b.nonzeros[n])), m_p);
        n++;
    }

    ret /= a.GetBinNum();
    ret = pow(ret, 1 / m_p);

    double normx = a.GetNorm();
    double normy = b.GetNorm();
    if (normx == -1 || normy == -1)
        throw "Invalid normx and normy!";
    if (normy + normx > EPSILON) { ret /= (normx + normy + EPSILON); }
    // if two of them are zero
    if (ret > 1) {
        cout << "[Info] Error ret > 1" << endl;
        cout << ret << " " << normx << " " << normy << endl;
        a.Print();
        b.Print();
        throw "Error, ret > 1";
    }
    return 1 - ret;
}

PNorm::PNorm(double p) : m_p(p) {
    if ((p < 1 && p > 0) || (p < 0 && p > -1) || (p < -1)) {
        cout << "[Error] Invalid p=" << p << " for calculating p-norm. Valid range [1, +Inf) & {0,-1}" << endl;
        throw "Error in calculate p-norm! Program will exit.";
    }

    if (p == -1 || p == 0) {
        cout << "[Alert] L" << p << "-Norm is not normalized to [0,1] range" << endl;
    }
}

double PNorm::calcNorm(BinningPeakList &a) {
    double ret = 0;
    for (int i = 0; i < a.nonzeros.size(); ++i) {
        ret += pow(fabs(a.GetIntensity(a.nonzeros[i])), m_p);
    }
    ret /= a.GetBinNum();
    ret = pow(ret, 1 / m_p);
    a.setNorm(ret);
    return ret;
}

SimlarityMetric::SimlarityMetric() = default;

SimlarityMetric::~SimlarityMetric() = default;

double SimlarityMetric::calcDistance(vector<double> &a, vector<double> &b) {
    cout << "Base function called" << endl;
    return 0;
}

double SimlarityMetric::calc(BinningPeakList &a, BinningPeakList &b) {
    cout << "Base function called" << endl;
    return 0;
}

double SimlarityMetric::calcNorm(BinningPeakList &a) {
    cout << "Base function called" << endl;
    return 1;
}

double SimlarityMetric::calcDistance(PeakList &a, PeakList &b) {
    BinningPeakList *x = a.CreateBinningPeaks();
    BinningPeakList *y = b.CreateBinningPeaks();

    double dot = this->calc(*x, *y);
    return dot;
}

DotProduct::DotProduct() = default;

double DotProduct::calcDistance(vector<double> &a, vector<double> &b) {
    double XX = 0, YY = 0, XY = 0;
    for (int i = 0; i < a.size(); i++) {
        XX += a[i] * a[i];
        XY += a[i] * b[i];
        YY += b[i] * b[i];
    }
    double sim = 0;

    if (fabs(XX) > EPSILON && fabs(YY) > EPSILON)
        sim = XY / sqrt(XX * YY);

    return sim;

}

double DotProduct::calcNorm(BinningPeakList &a) {
    double ret = 0;
    for (int i = 0; i < a.nonzeros.size(); i++) {
        double x = a.GetIntensity(a.nonzeros[i]);
        ret += (x * x);
    }
    ret = sqrt(ret);
    a.setNorm(ret);
    return ret;
}

double DotProduct::calc(BinningPeakList &a, BinningPeakList &b) {
    double dot = 0;
    vector<int> v;
    int m = 0, n = 0;
    while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
        if (a.nonzeros[m] == b.nonzeros[n]) {
            v.push_back(a.nonzeros[m]);
            m++;
            n++;
        }
        else if (a.nonzeros[m] < b.nonzeros[n]) m++;
        else n++;
    }

    for (int i : v) {
        dot += (b.GetIntensity(i) * a.GetIntensity(i));
    }
    double normx = a.GetNorm();
    double normy = b.GetNorm();

    double ret = 0;
    if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
        ret = dot / (normx * normy);
    }

    if (ret < 0) {
        cout << "x, y ,dot " << normy << " " << normx << " " << dot << endl;
    }
    return ret;
}

double SPC::calcDistance(vector<double> &a, vector<double> &b) {
    double ret = 0, sharedcount = 0;
    double acount = 0, bcount = 0;
    for (int i = 0; i < a.size(); ++i) {
        if (fabs(a[i]) > EPSILON && fabs(b[i]) > EPSILON) {
            sharedcount += 1;
        }
        if (fabs(a[i]) > EPSILON) {
            acount += 1;
        }
        if (fabs(b[i]) > EPSILON) {
            bcount += 1;
        }

    }
    ret = sharedcount / (acount + bcount - sharedcount);

    return ret;
}

SPC::SPC() = default;

double SPC::calc(BinningPeakList &a, BinningPeakList &b) {
    double ret = 0, sharedcount = 0;
    double acount = a.nonzeros.size(), bcount = b.nonzeros.size();

    int m = 0, n = 0;
    while (m < acount && n < bcount) {
        if (a.nonzeros[m] == b.nonzeros[n]) {
            sharedcount += 1;
            m++;
            n++;
        }
        else if (a.nonzeros[m] < b.nonzeros[n]) {
            m++;
        }
        else { n++; }
    }
    ret = sharedcount / (acount + bcount - sharedcount);
    return ret;
}

double SPC::calcNorm(BinningPeakList &a) {
    return a.nonzeros.size();

}

JSD::JSD() = default;

double JSD::calcDistance(vector<double> &a, vector<double> &b) {
    return 1;
}

double JSD::calcNorm(BinningPeakList &a) {
    return 1;
}

double JSD::calc(BinningPeakList &a, BinningPeakList &b) {
    //cout << "JSD" << endl;
    double div = 0;
    double sum_p = 0;
    double sum_q = 0;
    int m = 0;
    int n = 0;
    vector<double> p;
    vector<double> q;
    for (int i = 0; i < a.nonzeros.size(); i++) {
        double x = a.GetIntensity(a.nonzeros[i]);
        sum_p += x;
        p.push_back(x);
    }
    for (int i = 0; i < b.nonzeros.size(); i++) {
        double x = b.GetIntensity(b.nonzeros[i]);
        sum_q += x;
        q.push_back(x);
    }
    while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
        if (a.nonzeros[m] < b.nonzeros[n]) {
            double p_m = p[m] / sum_p;
            div += p_m * log(2) / 2;
            m++;
        }
        else if (a.nonzeros[m] > b.nonzeros[n]) {
            double q_n = q[n] / sum_q;
            div += q_n * log(2) / 2;
            n++;
        }
        else {
            double p_m = p[m] / sum_p;
            double q_n = q[n] / sum_q;
            double mean = (q_n + p_m) / 2;
            div += (p_m * log(p_m / mean) + q_n * log(q_n / mean)) / 2;
            m++;
            n++;
        }
    }
    while (m < a.nonzeros.size()) {
        double p_m = p[m] / sum_p;
        div += p_m * log(2) / 2;
        m++;
    }

    while (n < b.nonzeros.size()) {
        double q_n = q[n] / sum_q;
        div += q_n * log(2) / 2;
        n++;
    }
    return 1 - div;

}
