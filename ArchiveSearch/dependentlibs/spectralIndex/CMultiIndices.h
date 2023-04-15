//
// Created by wulong on 1/25/19.
//

#ifndef MYTOOL_CMULTIINDICES_H
#define MYTOOL_CMULTIINDICES_H

#include <iostream>
#include <vector>
#include <memory>
using namespace std;

class ICIndexWrapper;
class CPQParam;
class DataFile;
class ICQuery;
class CMzFileReader;
class Progress;
class PeptideProphetParser;

// create a class to randomly choose vectors from many vectors
class CShuffledVector
{
	vector<long> m_shuffledIdx;
public:
	CShuffledVector(long specnum);
	void getVector(float* src, int offset, int d, int n, float * dest);
};

class CRecallStat {
public:
    long m_total;
    long m_correct;

    CRecallStat();
    void verify();
    void correct();
    double recall();
    void show();
};

class MultipleIndicesImpl {
	bool m_useCpu;
	vector<shared_ptr<ICIndexWrapper>> m_ShuffleIndex;
    int m_numIndexUsed;
	vector<int> m_seeds;
    int m_nprobe;
public:
    bool usingCPU();
	MultipleIndicesImpl();
	void display();
	long ntotal(int i);
	void read();
    void removeIds(vector<long> &idx);
	void write();
	bool ifAllFileExist();
	int getNum() const;

    void setNum(int indexNum);
	int getNumOfIndexBuilt() const{return m_ShuffleIndex.size();}
	void setNprobe(int nprobe);
    int getNprobe(){return m_nprobe;}
	void train(long specnum, int dim,float *vec);
	void createEmptyIndices(vector<string> &indexstrs, int dim);
	void toGPU();
	void toMultipleGPUs(vector<int> &gpu_idx);
	void toCPU();
	void initialize(int numIndex, string indexpath, string basename, bool useMyOwn, shared_ptr<CPQParam> option, string indexshuffleseeds);

	void add(float *vec, long specnum);
    void getANNByShuffleQueries(int numQuery, int ret_num,
                                const float *vquery, vector<vector<long>> &multipleInd,vector<vector<double>> &results_dist,
                                int i, vector<float> &dist, vector<long> &ind);

private:
	int getDim(int i);
	void search(int i, int ret_num, vector<float> &dist, vector<long>& ind, float *shuffledVquery, int numQuery) ;
	void shuffleVectors(int i, int numVec, float * p);
	bool istrained(int i);
	string getIndexFileName(int i);
	void add(int i, int batchsize, float *shuffleVec);
    ICIndexWrapper * getShuffledIndex(int i);
	void trainIndex(int i, int batchSize, float * vecBatch);
	void readIndex(int i);
	void readIndex(int i, string filename);
	void addshuffle(int i, float *vec, long specnum);

	void writeIndex(int i);
	void writeIndex(int i, string filename);
	bool ifFileExist(int i);
	void display(int i);
	void createEmptyIndex(int i, string indexstr, int dim);
	void toGPU(int i);
	void toMultipleGPUs(int i, vector<int> &gpu_idx);
	void shufflevector(int seed, float * p, int dim);
};

class CMultiIndices {
    vector<string> m_multiIndicesStr;
    bool m_doShuffle;
	MultipleIndicesImpl m_impl;
	double m_tolerance;
	bool m_removeprecursor;
	bool m_useFlankingBins;
	shared_ptr<CPQParam> option;
	int m_indexNum;
	int m_topPeakNum;
	int m_dim;
	const long TRAINING_SPEC_NUM;
    int m_nprobe_default;
public:
	bool m_debug;

public:
    void removeIds(vector<long> &idx){

        m_impl.removeIds(idx);
    }
    // get parameter ready, but not really create the index object.
    CMultiIndices(string indexstring, string indexshuffleseeds,string indexpath, string libname, bool doShuffle,
                  double tolerance, bool useMyOwn, shared_ptr<CPQParam> cpqParam, const int topPeakNum, bool removePrecursor, bool useFlankingBins, int dim);
    ~CMultiIndices();

    int getNum();
	int getNumOfindexInUse(){return m_impl.getNum();}
    void create(long &specnum, DataFile &splib, bool removeprecursor, bool useFlankingBins);
	void trainOnSingle(int dim, long &specnum, DataFile &splib, bool removeprecursor, bool useFlankingBins);
	bool isExist(){return m_impl.ifAllFileExist();} // to be called : todo
	void trainOnFileList(string datafilelist);
	void appendList(string datafilelist);
	void loadIndices(){
	    m_impl.read();
	}

    void append(DataFile &df);
	void setNprobe(int nprobe, bool verbose);
    void toGpu();
	void toMultipleGPUs(vector<int> &gpu_idx){
		m_impl.toMultipleGPUs(gpu_idx);
	}
    void toCpu();
    void display();
    void getAnns(ICQuery& q, int ret_num, vector<vector<long>> &results, vector<vector<double>> &results_dist, int indexNum);
    void recallOfAnn(DataFile &df, string ipropepxmlfilename, CMzFileReader &compactRawData, bool useflankingbin, int dim,
                     DataFile &splib, bool use_gpu, int indexChoice, int batchsize, string outputpath, int ret_num,
                     long spec_start, long spec_end);

    // to be optimized.
    void collectANNs(int indexChoice, int ret_num, const vector<vector<long>> &results,  const vector<vector<double>> &results_dist,int idx,
                     vector<long> &int_ind,  vector<double> &dist, bool verbose) const;


	long total();
	void write();
private:
    string getOutFileName(const string &ipropepxmlfilename, const string &outputpath) const;

    void verifyEachQuery(DataFile &df, DataFile &splib, int indexChoice, int ret_num, PeptideProphetParser &ppp,
                         double probThreshold, Progress &ps, const vector<vector<long>> &results,const vector<vector<double>> &result_dists,
                         const vector<long> &queries, int indexNum, ofstream &fout, CRecallStat &recallANNs) const;
    int findPeptide(DataFile &splib, const vector<long> &int_ind, const string &peptide) const;
    void printIndicesStr();
    void getAnnsForQueries(int numQuery, int ret_num, float *vquery, vector<vector<long>> &multipleInd,vector<vector<double>> &results_dist,
                           bool verbose);
    void train(int dim, long &specnum, float *vec) ;

    float * collectTrainingSpectra(const vector<string> &files,long & total) const;
};


#endif //MYTOOL_CMULTIINDICES_H
