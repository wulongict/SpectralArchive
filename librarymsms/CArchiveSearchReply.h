//
// Created by wulong on 4/13/19.
//

#ifndef MYTOOL_CARCHIVESEARCHREPLY_H
#define MYTOOL_CARCHIVESEARCHREPLY_H

#include <set>
#include <queue>
#include <memory>
#include <ostream>

using namespace std;
class SPsmAnnotation;
class CAnnotationDB;
class ICQuery;

struct SAnnSpectrum {
    long idx;
    float dist;
    int dotprod;
    float appdist;
    float pvalueAll;
    vector<float> pvalues;
    float pvaluePartial;
    bool neighborPvalueAll; // whether using neighbor calculating pvlaue;
    vector<float> spec;
    void print();
    void setpvaluePartial(float pvaluepartial);
    float getStdOfPvalue() const;
    SAnnSpectrum();
    SAnnSpectrum(long idx_, float dist_, float appdist_);
    SAnnSpectrum(long idx_, float dist_, float appdist_, int dotprod_);
    SAnnSpectrum(long idx_, float dist_, float appdist_,  vector<float>::iterator start,
            vector<float>::iterator end );
};

struct CompareEntryDist{
    bool operator()(const SAnnSpectrum &a, const SAnnSpectrum &b);
};

class CAnnSpectra {
    bool m_has_spec;
    int m_dim;
    long m_queryidx;
    long m_tophitIdx; // groundtruth

public:
    vector<SAnnSpectrum> m_anns;
    priority_queue<SAnnSpectrum, vector<SAnnSpectrum>,CompareEntryDist> m_queueneighbors;
public:
    CAnnSpectra(int dim, vector<long> &ind, vector<float> &dist, vector<float> &specs,
                long queryidx);
    CAnnSpectra(vector<long> &neighborIdx, vector<float> &dist,
                long queryidx,vector<int> &dpscore);
    CAnnSpectra(vector<long> &ind, vector<float> &dist, vector<float>& appdist,
                long queryidx);
    vector<SAnnSpectrum>::iterator find(long query_index);

    bool isExist(long query_index);
    bool empty() const;
    long getQueryIdx(){return m_queryidx;}
    void setTopHitIdx(long tophit){
        m_tophitIdx = tophit;
    }
    CAnnSpectra &keepTopN(int n, bool verbose);
    CAnnSpectra &keepWithMinDP(double mindp);
    void print();
    vector<float> ConcatnateSpec();

    vector<long> ConcatenateIndex();

    void toLinkStr(string &links);
    int size() const;

    string createUpdateNeighborSQL();
private:
    CAnnSpectra & print(int i);
    void tolink(ostringstream &oss, vector<SAnnSpectrum>::iterator &it) const;
};

class CArxivSearchResult {
    vector<CAnnSpectra *> m_queryresults;
public:
    CArxivSearchResult(){}
    int size();
    void copyto(vector<CAnnSpectra*> &vqrs, int offset);
    // copy pointer to an new vector and change the original pointer to nullptr.
    void moveTo(vector<CAnnSpectra*> &vqrs, int offset);
    void release(int i);
    vector<long> concatenateIndex(int i);
    void keepTopN(int topN, bool verbose);
    void keepWithMinDP(double mindp);
    void keepWithMinDP(int i, double mindp);
    void keepTopN(int i, int topN, bool verbose);
    void push_back(CAnnSpectra *p);
    void set(int i, CAnnSpectra *p, bool releasefirst=true);
    ~CArxivSearchResult();
    CAnnSpectra * get(int i );
    void init(ICQuery &query, bool verbose, vector<vector<long>> &allRetIdx, vector<vector<float>> &accDist,vector<vector<int>> &dpscores);
};

class CNodeGtInfo {
    shared_ptr<SPsmAnnotation> m_gtinfo;
    int m_group;
public:
    CNodeGtInfo(const SPsmAnnotation& gtinfo);
    string toNodeStr( );
    string getPepSeq();
    void setGroup(int group);
};


class CArchiveSearchReply{
    vector<CAnnSpectra*> &m_vqr;
    long m_query_index;
public:
    CArchiveSearchReply(vector<CAnnSpectra*> &vqr, long queryIndex);
    void tojsonstring(string &jsonstring, CAnnotationDB &m_AnnotationDB);
private:
    void toNodeStr(string &myString, CAnnotationDB &m_AnnotationDB);
    void toLinkStr(string &links);
    void getUniqueId(set<long> &nodeIdx, bool verbose);
    void createNodeList(set<long> &nodeIdx,CAnnotationDB &m_AnnotationDB, vector<CNodeGtInfo> &vNodes);
    void getNodeStr(vector<CNodeGtInfo> &vNodes, string &myNodes);
};

#endif //MYTOOL_CARCHIVESEARCHREPLY_H
