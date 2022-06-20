#ifndef PEAKLIST_H_INCLUDED
#define PEAKLIST_H_INCLUDED

#include <vector>
//#include <string>
//#include <iostream>
//#include <thread>
//#include <algorithm>
//#include <chrono>
//#include <cfloat>
//#include <cstring>
//#include <fstream>
//#include <iomanip>
//#include <sstream>
//#include <cmath>
//#include <map>
//#include <numeric>
//#include <iterator>     // std::ostream_iterator


using namespace std;



// A new class for peaks transform/Filter
//
class BinningPeakList;

class VectorFilter {
    string m_Type;
public:
    VectorFilter();
    void setType(string type);
    string getType();
    virtual ~VectorFilter();
    virtual void run(BinningPeakList *x);
};

class PeakSquareRoot : public VectorFilter {
public:
    PeakSquareRoot();
    void run(BinningPeakList *x);
};

class NormalizedToProb : public VectorFilter {
public:
    NormalizedToProb();
    void run(BinningPeakList *x);
};

class RemovePeakLessThanMean : public VectorFilter {
public:
    RemovePeakLessThanMean();
    void run(BinningPeakList *x);
};

class RemovePeakLessThanMedian : public VectorFilter {
public:
    RemovePeakLessThanMedian();
    void run(BinningPeakList *x);
};

class KeepTopNPeaks : public VectorFilter {
private:
    int m_topN;
public:
    KeepTopNPeaks();
    void setTopN(int topN);
    void run(BinningPeakList *x) override;
};

// FactoryMethod
VectorFilter *VectorFilterFactoryMethod(const string& FilterType);

class BinningPeakList {
    vector<double> m_intensityList;
    int BinSize;
    double norm;
    int peakcountOutOfRange;
public:
    BinningPeakList();
    BinningPeakList(const vector<double> &x, const vector<double> &y, bool useFlankingBins, bool verbose=false);
    ~BinningPeakList();
    void setNorm(double x);
    double CalcNorm();// error here we will stop using this one
    double CalcDotProductByNonZeros(BinningPeakList &other);
    double CalcSimilarByPnorm(BinningPeakList &other, double p);
    int GetBinNum();
    double GetIntensity(int k) const;
    vector<double> GetIntensityList() const;
    void setIntensity(int i, double intensity);
    void setIntensitylist(vector<double> intensitylist);
    BinningPeakList *filterPeaks(VectorFilter *vf);
    void Print();
    double GetNorm();
public:
    vector<int> nonzeros;
};

class PeakList {
private:
    double m_RTinSeconds;
    vector<double> m_mzList;
    vector<double> m_intensityList;
    BinningPeakList *m_binningPeakList; // release the resources in the destruction function
public:
    PeakList();
    PeakList(const PeakList &pl);
    ~PeakList();
    void print();
    void removePeaksWithin(double mz_left, double mz_right);
    double getRTinSeconds() const;
    void setRTinSeconds(double RTinSeconds);
    void setM_mzList(const vector<double> &mzList);
    void setM_intensityList(const vector<double> &intensityList);
    BinningPeakList *getM_binningPeakList() const;
    vector<double> getM_mzList() const;
    vector<double> getM_intensityList() const;
    void InsertPeak(double mz, double intensity);
    void setBinningPeakList(BinningPeakList *binlist);
    PeakList getPeaksWithin(double leftmz, double rightmz);
    void addPeakList(PeakList &pl, int ms2_number = 0);
    double CalcDotProduct(PeakList &other);
//    double CalcDotProduct(PeakList &other, SimlarityMetric * simcalculator)
    BinningPeakList *CreateBinningPeaks(bool useFlankingBins=false, bool verbose=false);
    void NormalizedToSum();
    void KeepTopN(int N);
    void rankingAsIntensity(int maxRanking=50);
};

class AnalysisMassDiff {
public:
    static double minDiff(double mz, const vector<double>& mzlist);
    std::vector<double> CalculateMassDifference(std::vector<PeakList *> &vpl, double tolerance);
};

void releaseVectorPeakListPtr(vector<PeakList *> &vpl);
class PeakListFilter {
    string m_FilterType;
public:
    const string &getM_FilterType() const;
    void setM_FilterType(const string &FilterType);
public:
    PeakListFilter();
    virtual ~PeakListFilter();
    virtual void run(PeakList *x);
};


class PeakListKeepTopN : public PeakListFilter {
private:
    int m_topN;
public:
    PeakListKeepTopN();
    void setTopN(int topN);
    void run(PeakList *x);
};

class PeakListNormalizedByMaximal : public PeakListFilter {
public:
    PeakListNormalizedByMaximal();
    void run(PeakList *x);
};

class SimlarityMetric {
public:
    SimlarityMetric();
    virtual ~SimlarityMetric();
    virtual double calcDistance(vector<double> &a, vector<double> &b);
    virtual double calc(BinningPeakList &a, BinningPeakList &b);

    virtual double calcNorm(BinningPeakList &a);

    double calcDistance(PeakList &a, PeakList &b);
};

// Remember: This is cosine of two vector, not simply dot-product.
// Tested by gtest from Google
class DotProduct : public SimlarityMetric {
public:
    DotProduct();

    double calcDistance(vector<double> &a, vector<double> &b);

    double calcNorm(BinningPeakList &a);

    double calc(BinningPeakList &a, BinningPeakList &b);
};

//
// This is the pearson's correlation coefficient
// ToDo;  to be tested.
class PCC : public SimlarityMetric {
public:
    PCC();

    double calcDistance(vector<double> &a, vector<double> &b) override;

    double calcMean(BinningPeakList &a);

    double calcNorm(BinningPeakList &a);

    double calc(BinningPeakList &a, BinningPeakList &b);
};



// Todo: to be tested by Gtest.
// ToDo: new scores to be added
// 1. Shared peak counts:
// 2. K-L divergence.
// 3. Pearson Correlaiton Coefficient: r

// Todo: This is not real Minkowski distance,
// We add 1/n in the formula
// (sum(xi^p)/n)^(1/p)
// This term is removed when we normalize the right-hand side, divided by sum of p-norm(a) and pnorm(b)
class PNorm : public SimlarityMetric {
    // If each feature/peak contains some information, then combining all those informations by p-norm might lead to a very powerful discriminator.
    // Here we use Lp-distance
    // The distance is bounded by Minkowski inequality
    // pNorm(f+g) <= pNorm(f) + pNorm(g)
    // A problem is: if f and g are of very different magnitude, then, f could be the dominate factor in pNorm(f+g)
    //
private:
    double m_p;
    // valid for
    // m_p >= 1; Regular p-Norm
    // m_p = 0;  Geometric distance
    // m_p = -1;  Harmonic distance
public:
    PNorm(double p = 2);

    // Todo: check if a.size == b.size
    // if not, exit
    double calcDistance(vector<double> &a, vector<double> &b);

    double calcNorm(BinningPeakList &a);

    double calc(BinningPeakList &a, BinningPeakList &b);
};

class SPC : public SimlarityMetric {
public:
    SPC();

    double calcDistance(vector<double> &a, vector<double> &b);

    double calc(BinningPeakList &a, BinningPeakList &b);

    double calcNorm(BinningPeakList &a);

};

class JSD : public SimlarityMetric {
public:
    JSD();

    double calcDistance(vector<double> &a, vector<double> &b);

    double calcNorm(BinningPeakList &a);

    double calc(BinningPeakList &a, BinningPeakList &b);
};

PeakListFilter *CreatePeakListFilterFactoryMethod(string param);

#endif // PEAKLIST_H_INCLUDED
