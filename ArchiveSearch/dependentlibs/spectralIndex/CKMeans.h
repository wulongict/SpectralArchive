//
// Created by wulong on 4/16/19.
// Using library Eigen!

#ifndef MYTOOL_CKMEANS_H
#define MYTOOL_CKMEANS_H

#include <memory>
#include <random>
#include "Core"
#include "Dense"
#include <sstream>
using namespace std;

// TODO: we need a factory method to create the objects.

class ICKMeansVec{
public:
    virtual ~ICKMeansVec(){}
    virtual double get(int i) = 0;
    virtual int size() = 0;
    virtual void set(int i, double val) = 0;

    void initAsZero();
    void display();
    double distance(ICKMeansVec &other);
    ICKMeansVec& operator+=(ICKMeansVec &other);
    double operator[](int i);
    ICKMeansVec & operator/=(double x);
    friend ostream & operator << (ostream &os, ICKMeansVec &other)    {
        for(int i = 0; i < other.size(); i ++)        {
            if(i!=0) os << "\t";
            os << other.get(i) ;
        }
        return os;
    }
};


class ICData
{
public:
    virtual ~ICData(){}
    virtual int size() = 0;
    virtual ICKMeansVec & getVec(int i ) = 0;

    ICKMeansVec& operator[](int i);
    void average(vector<int> &idxlist, ICKMeansVec &out);
    int getClusterID(ICKMeansVec &sample);
    void display();
};

class CVec: public ICKMeansVec {
    vector<double> &m_vec;
public:
    CVec(vector<double> &v);
    ~CVec(){}
    double get(int i ) override;
    void set(int i, double val) override;
    int size() override;
};

class CVecSelected: public ICKMeansVec
{
    vector<double> &m_vec;
    vector<int> &m_idx_selected;
public:
    CVecSelected(vector<double> &v, vector<int> &idx_selected);
    ~CVecSelected(){}
    double get(int i ) override;
    void set(int i, double val) override;
    int size() override;
};

class CVecFloatPtr:public ICKMeansVec
{
    float *m_ptr;
    int m_dim;
public:
    CVecFloatPtr(float *ptr, int dim): m_ptr(ptr), m_dim(dim){}
    double get(int i) override{return m_ptr[i];}
    void set(int i, double val) override{m_ptr[i]=val;}
    int size() override{return m_dim;}
};

class CData: public ICData
{
    vector<vector<double>> &m_data;
    vector<CVec> m_Vecs;
public:
    CData(vector<vector<double>> &data);
    ICKMeansVec & getVec(int i) override;
    int size() override ;
};


class CDataSelected: public ICData
{
    vector<int> m_idx_selected;
    vector<vector<double>> m_data;
    vector<CVecSelected> m_Vecs;
public:
    CDataSelected(vector<vector<double>> &data, vector<int> idx_selected);
    ICKMeansVec & getVec(int i) override;
    int size() override;
};

class CDataFloatPtr: public ICData
{
    float *m_ptr;
    int m_num;
    int m_dim;
    vector<CVecFloatPtr> m_Vecs;
public:
    CDataFloatPtr(float *p, int num, int dim):m_ptr(p), m_num(num),m_dim(dim)
    {
        for(int i = 0; i < m_num; i ++)        {
            m_Vecs.emplace_back(CVecFloatPtr(m_ptr + i*m_dim, m_dim));
        }
    }
    ICKMeansVec & getVec(int i) override {
        return m_Vecs[i];
    }
    int size(){return m_Vecs.size();}
};


void kmeans(ICData *data, ICData *centroids);

class CKMeans {
    ICData *_data;
    ICData *_centroids;
    vector<int> m_membership;
    int m_k;
    int m_num;
    bool m_release;
    bool m_verbose;
public:
    CKMeans(vector<vector<double>> &data, vector<vector<double>> &centroids, int k, bool verbose);
    CKMeans(ICData *data, ICData *centroids, bool verbose);
    CKMeans(shared_ptr<ICData> data, shared_ptr<ICData> centroids, bool verbose);
    ~CKMeans();
    void run();
    void getResidual(vector<vector<double>> &res);
    bool update_membership(int i, double &dist);
    void getNearestId(vector<int> &ids);

private:
    void update_centroid(int j);
    void update_centroids();
    bool update_membership(double &totaldist);
    void display();
};




struct CPQResultEntry{
    int m_idx;
    double m_dist;
    CPQResultEntry(int idx, double dist);
    CPQResultEntry();
    friend ostream & operator << (ostream &os, const CPQResultEntry & cpqrEntry);
    void set(int idx, double dist);
};
class CPQResult{
    vector<vector<CPQResultEntry>> m_res;
    int m_topN;
public:
    CPQResult(int topN){m_topN = topN;}
    ~CPQResult(){}
    void add(vector<int> &idx, vector<double> &dist);
    void display();
    void add(vector<CPQResultEntry> &x);
    const vector<CPQResultEntry> & get(int i) const;
    int size() const;
    bool operator==(const CPQResult &other) const;
    bool equal(const CPQResult & other);
    void flattern( vector<float>& dist, vector<long>& ind, int d);
};


class CPQParam{
public:
    enum PQ_OPTIONS {INITIAL_ONES, RANDOM_SAMPLES, RANDOM, RANDOM_PROJECTION_KMEANS};
    PQ_OPTIONS _option;
    CPQParam(){
        _option = INITIAL_ONES;
    }
    CPQParam(PQ_OPTIONS option){
        _option = option;
    }
    CPQParam(const CPQParam &other)    {
        _option = other._option;
    }
    CPQParam & operator=(const CPQParam &other)    {
        _option = other._option;
        return *this;
    }
    void display();
};

class CIndexCodeBook
{
    vector<int> m_indexCodeBook;
    int m_dim;
public:
    CIndexCodeBook(int dim):m_dim(dim){
    }
    CIndexCodeBook():CIndexCodeBook(-1){}
    void setDim(int dim)    {
        m_dim = dim;
    }
    long size(){return m_indexCodeBook.size()/m_dim;}
    void extend(int num){
        long len = size();
        m_indexCodeBook.resize((len + num)*m_dim);
    }
    inline long get(long i, long j)    {
        return m_indexCodeBook.at(i*m_dim + j);
    }
    inline void set(long i, long j, int val)    {
        m_indexCodeBook[i*m_dim +j] = val;
    }
    inline int getDim(){return m_dim;}
};

class CPreComputeDist
{
    vector<vector<double>> &m_initCenMatrix;
    vector<vector<double>> &m_resCenMatrix;
    vector<Eigen::MatrixXf> m_preDist;
    int m_dim;
    int m_clusterNum;
    Eigen::MatrixXf m_resM;
    Eigen::MatrixXf m_resCenM;
    int m_subspaceNum;
public:
    CPreComputeDist(vector<vector<double>> & initCen,
            vector<vector<double>> & resCen,int subspaceNum);
    float get(int iniCenId, int resCenId, int spaceId)    {
        return m_preDist[spaceId](iniCenId, resCenId);
    }

    void set(int iniCenId,int resCenId, int spaceId, float s)    {
    }

    void create(shared_ptr<ICKMeansVec> &a);
};

class CProductQuantization
{
    shared_ptr<ICData> m_initCenPtr; // in proper result
    ICData *m_data;
    int m_subSpaceNum;
    int m_clusterNum;
    int m_vectorSize;
    vector<vector<double>> m_initCenMatrix;
    vector<vector<double>> m_residualCenMatrix;

    bool istrained;
    int m_nprobe;
    CIndexCodeBook m_indexCodeBookObj;
    vector<float> m_precomputeDist;
    vector<vector<long>> m_idsInEachBuckets;
    shared_ptr<CPreComputeDist> m_preCompDistPtr;
public:
    void initPreComputeDist();
    inline float getPreCompDist(int iniCenId, int resCenId, int spaceId );
    CProductQuantization(int PQN,int IVF, int DIM);
    CProductQuantization();
    void createIdsInEachBuckets();
    void createPreComputeDist(shared_ptr<ICKMeansVec> &a);
    ~CProductQuantization();
    void train(ICData *data, CPQParam &cpqParam);
    void display(bool all);

    void add(int numspec, float *ptr, int d);
    void load(string filename);
    void save(string filename);

    int getCodeBookDim() {
        return m_indexCodeBookObj.getDim();
    }

    void resizeCodeIndex(int newnum)    {
        m_indexCodeBookObj.extend(newnum);
    }

    void setIndexCodeBook(long i, long j, int val)    {
        m_indexCodeBookObj.set(i,j,val);
    }

    int getFromIndexCodeBook(long i,long j)    {
        return m_indexCodeBookObj.get(i,j);
    }
    
    void search(float *queries, int querynum, CPQResult &res);

    // we need to go to lab now. as library will close now
    void calculateDistWithCodeBook(int idx, vector<int> &idx_list, vector<double> &dist_list);

    vector<double> reconstruct(vector<int> &code);
    long size(){
        return m_indexCodeBookObj.size();
    }
    bool isTrained(){return istrained;}
    int getDim(){return m_vectorSize;}
    void setnProbe(int nprobe){m_nprobe = nprobe;}
    int getClusterNum(){return m_clusterNum;}
    void chooseProperInitialVec(ICData *data, CPQParam &cpqParam);
private:
    vector<int> getIdx(int i );
    void useFirstFewAsCentroid(ICData *data);
    void useRandomVecAsCentroid(ICData *data);
    void useGaussianRandomVecAsCentroid(ICData *data);
    void useRandomProjection(ICData *data);
};


#endif //MYTOOL_CKMEANS_H
