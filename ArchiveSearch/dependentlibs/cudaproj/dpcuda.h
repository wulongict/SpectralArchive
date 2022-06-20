//
// Created by wulong on 7/15/18.
//

#ifndef MYTOOL_DPCUDA_H
#define MYTOOL_DPCUDA_H

#ifdef __CUDA__
#include <cuda_runtime_api.h>
#endif

#include <iostream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../../../librarymsms/ICMzFile.h"
#include <stdio.h>

#define cudaCheckErrors(msg) \
do {\
    cudaError_t __err = cudaGetLastError(); \
    if(__err != cudaSuccess){\
        fprintf(stderr, "Fetal error: %s (%s at %s: %d)\n", \
        msg, cudaGetErrorString(__err), \
        __FILE__, __LINE__);\
        fprintf(stderr, "*** FAILED - ABORTING\n");\
        exit(1);\
    }\
} while(0) 


using namespace std;
class ICQuery;
class ICMzFile;
class ICMzFactory;

class CUDAScore : public ICMzFile{
    string m_mzfile;
    const int m_PeakNumPerSpec;
    long m_filebytes;
    long m_peaknum;
    long m_specnum;
    uint16_t  *m_all; // spectral archive CPU/GPU shared
    uint16_t  *m_s; // score              CPU/GPU shared
    uint16_t *vecform;// vector form of query   To speedup
    long *gpuindexlist;
    int m_device;
    int m_max_query_index_size;
public:
    CUDAScore(string mzfile);
    ~CUDAScore();

	void calcDotProduct(int TopNPeak, int tol, uint16_t * queryspec,
                        int blockSize, vector<long> &indexlist, vector<int>&scores) override;
    long getSpecNum() override;
    string getClassName() override {return "CUDAScore";}
    int getPeakNumPerSpec() const override;
    uint16_t *getSpecBy(long queryindex) override;

    void Query(int topn, int tol, int queryindex, int blockSize);
    void scoreAllVecForm(int mzTopN, int tol, long queryindex, int blockSize, bool normalize);
    void queryFastOnIndex(int TopNPeak, int tol, long queryindex, int blockSize, vector<long> &indexlist);
    void QueryfastOnIndexWithUINT16(int TopNPeak, int tol, uint16_t * queryspec, int blockSize, vector<long> &indexlist);
    uint16_t * getscore();
//    vector<float> dist(long query_index, vector<long> &ind, int tol);
    vector<int> distributionAll(int tol, long queryindex, bool normalize);
    void normalizationAllScore(long queryindex);
    int initGpuIndexList(vector<long> &indexlist);
    void toVectorForm(uint16_t* queryspec,int tol,uint16_t* vecform, int len);
    void copyScore(vector<long> &indexlist, vector<int> &scores);
    void dpscore(double tolerance, vector<vector<long>> &allRetIdx, int threadnum, vector<vector<float>> &accDist,
                 ICQuery &query,vector<vector<int>> &dpscores);

};


class CMzCUDAFactory: public ICMzFactory{
public:
    virtual shared_ptr<ICMzFile> create() override {
        cout << "[Info] Creating GPU scorer Object" << endl;
        return make_shared<CUDAScore>(getmzfilename());
    }
};


#endif //MYTOOL_DPCUDA_H
