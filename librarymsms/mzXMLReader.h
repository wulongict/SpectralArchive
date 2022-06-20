//
// Created by wulong on 7/12/15.
//

#ifndef PROJECT_MZXMLREADER_H
#define PROJECT_MZXMLREADER_H
//

#include <glob.h>
#include <ctime>
#include <zlib.h>
#include <random>

typedef std::string mzXMLFilename;
using namespace std;

class VectorFilter;
class PeakListFilter;
class SimlarityMetric;
class PeakList;

class mzXMLReader {
    int MaxScanDiff;
    string SLASH; // or
    string MSLEVEL;
    SimlarityMetric *m_SimCalculator;
    vector<VectorFilter *> m_PeakFilters;
    vector<PeakListFilter *> m_ms2PeakListFilters;
public:
    mzXMLReader(string mslevel, int max_scan_th = 10000);
    void setPeakFilters(string filters);

    void setMs2PeakListFilters(string ms2PeakFilters);

    void setSimlarityMetric(SimlarityMetric *simMetric);


    virtual ~mzXMLReader();

    vector<PeakList *> ReadmzXMLToPeakLists(const mzXMLFilename& f);

    void
    CalculateDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, int start, int end, int thread_ID, double *res) const;

    void SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputfile);

    void SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, double *res);

    //overload with vector of vector
    void SingleThreadDotProduct(vector<vector<double> > a, vector<vector<double> > b, double *res);


    void RandomVectorMultiplication(vector<PeakList *> &a, vector<vector<double> > &product, vector<double> random_vector,
                               int num_RT, int cycle);

    void CreateBinningList(vector<PeakList *> &a);

    void MultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile, int threadNum,
                               int cycle = 1);

    //
    void RandMultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile, int threadNum,
                                   int cycle);

    void ApplyMemoryforMatrix(const long aSize, const long bSize, long long int &nSize, double *&res);

    void CalcDotProduct(vector<PeakList *> &a, const vector<PeakList *> &b, double *res, int threadNum);


    void setMaxScanDiff(int _MaxScanDiff);

    void FilterMS2PeakList(PeakList *x);

    void exportToCSV(vector<PeakList *> &spec, const string& csvFileName);
};

vector<double> normalize(vector<vector<double> > a);

#endif //PROJECT_MZXMLREADER_H
