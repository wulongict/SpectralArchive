//
// Created by wulong on 10/3/19.
//

#ifndef MYTOOL_PVALUECALC_H
#define MYTOOL_PVALUECALC_H


#include <vector>
#include <iostream>
#include <memory>
using namespace std;
class CFragIndex;
class ICMzFile;
class ICQuery;
class CArxivSearchResult;
class SLinearRegressionModel;


class IndexListBuilder{
    float m_fraction;
    long m_MaxNum ;
public:
    IndexListBuilder(float fraction, long maxNum):m_fraction(fraction), m_MaxNum(maxNum){

    }
    IndexListBuilder() : IndexListBuilder(0.01, 10000){
    }
    void build(vector<long>& m_indexlist, long totalSpecNum, int  seed);
};

class CPValueCalculator {
    long m_total_spec_num;
    int m_threadnum;
    vector<long> m_indexlist;
    shared_ptr<CFragIndex> m_FragIdxPtr;
    int m_tol;
    bool m_saveBackgroundScore;


public:
    CPValueCalculator(long totalSpecNum, int threadnum, int seed, int tol, bool saveBackgroundScore,
                      bool verbose);
    ~CPValueCalculator();
    void buildFragIndex(shared_ptr<ICMzFile> &csa, bool verbose);

    // an different run?
    void run(vector<long> &idxListForQuery, ICQuery &query, CArxivSearchResult &vqr, bool plot, bool verbose,
             vector<int> *ptrDP, SLinearRegressionModel *ptrLR, vector<vector<int> >*ptrDP_q, vector<SLinearRegressionModel> *ptrLR_q);
    vector<long> getIndexList(){return m_indexlist;}
    void calcDotProducts(uint16_t *spec, vector<int> &score, bool normalize);
};

class CPvalueMultiCalculator
{
    vector<shared_ptr<CPValueCalculator>> m_pvalueCalc;
    int m_tol;
public:
    CPvalueMultiCalculator(int numCalc, long totalnum, int threadnum,
                           shared_ptr<ICMzFile> &pScorer, int startseed, int tol);
    // change it
    void run(vector<long> &idxListForQuery, ICQuery &query, CArxivSearchResult &vqr, bool plot, bool verbose ,vector<vector<int>> *ptrDPs,
             vector<SLinearRegressionModel> *ptrLRs,
             vector<vector<vector<int>>> *ptrDP_all_pv_query,
             vector<vector<SLinearRegressionModel>> *slRMs_pv_query);
    //vector<int> *ptrDP, SLinearRegressionModel *ptrLR);
};

#endif //MYTOOL_PVALUECALC_H
