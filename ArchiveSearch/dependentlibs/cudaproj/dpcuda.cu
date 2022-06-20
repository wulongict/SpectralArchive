//
// Created by wulong on 7/15/18.
//

#include <iostream>
#include <vector>
#include <list>

#include "dpcuda.h"
#include "host_defines.h"
#include "device_launch_parameters.h"
#include "cuda_runtime_api.h"
#include "cuda_profiler_api.h"
#include "thrust/host_vector.h"

#include "../../../librarymsms/ICMzFile.h"
#include "../../../librarymsms/ICQuery.h"
#include "../../../librarymsms/ProteomicsDataTypes.h"
#include "../../../librarymsms/XMLFileParser.h"
#include "../../../librarymsms/ICGtInfoUpdate.h"
#include "../../../librarymsms/Util.h"

// CUDA kernel to add elements of two arrays
__global__
void add(int n, float *x, float *y) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
        y[i] = x[i] + y[i];
}

// CUDA kernel to add elements of two arrays
__global__
void adduint(int n, uint16_t *x, uint16_t *y) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
        y[i] = x[i] + y[i];
}

__global__
void dp_norec(int specnum, uint16_t *all, int PeakNum, uint16_t *s, int query, int mzTopN, int tol) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if (index == 0) {
//        cout << "stride: " << stride <<" block id="<<  blockIdx.x << " block size: " << blockDim.x << " gridDim/#blocks: " << gridDim.x << endl;
//        printf("stride: %d, blockid: %d, blocksize: %d, gridDim: %d\n", stride, blockIdx.x, blockDim.x, gridDim.x);
    }

    uint16_t *queryvector = all + PeakNum * query;

    s[index] = 0;
    for (int k = index; k < specnum; k += stride) {
        uint16_t *x = all + PeakNum * query, *y = all + index * PeakNum;

        // todo: the algorithm can be improved!
        for (int i = 0; i < mzTopN; i++) {
            if (x[i] == 0) break;
            for (int j = 0; j < mzTopN; j++) {
                if (y[j] == 0) { break; }

                if (x[i] > y[j]) // x[i] bigger
                {
                    if (x[i] - y[j] < tol) {
                        s[index] += (PeakNum - i) * (PeakNum - j);
                        break;
                    }
                } else {
                    if (y[j] - x[i] < tol) {

                        s[index] += (PeakNum - i) * (PeakNum - j);
                        break;
                    }
                }
            }

        }
    }


}

// todo: create another function to use GPU to calculate dot product

__global__
void dp_norec_vecform(int specnum, uint16_t *all, int PeakNum, uint16_t *s, int mzTopN, int tol,
                      uint16_t *vecform, bool debug=false, long debug_queryindex=-100) {
    long index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    s[index] = 0;
    for (int k = index; k < specnum; k += stride) {
        uint16_t *y = all + index * PeakNum;
        unsigned long rec = 0;
        // todo: the algorithm can be improved!
        uint16_t used[51]={0};
        for (int i = 0; i < mzTopN; i++) {
            if(used[vecform[y[i]]]){}
            else
            {

                if(vecform[y[i]]) used[vecform[y[i]]] = 1;
                s[index] += vecform[y[i]] * (PeakNum - i) ;
#ifdef _DEBUG_INDEX_WITH_PRINT_
                if(debug_queryindex == index and debug and vecform[y[i]]>0)
                {
                    printf("%d -th mz: %d" ,i+1,y[i]);
                    printf("score : %d x %d -> %d\n", vecform[y[i]], PeakNum-i, s[index]);

                }
#endif


            }

        }
    }


}
#define _DEBUG_INDEX_WITH_PRINT_
// potential bug is here!
__global__
void dp_norec_vecform_on_index(int specnum, uint16_t *all, int PeakNum, uint16_t *s, int mzTopN,
                               int tol, uint16_t *vecform, long *gpuindexlist, int indexlistsize, bool debug=false, long debug_queryindex=-100) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int k = index; k < indexlistsize; k += stride) {
        s[k] = 0;
//        if(index==0) printf("topn=%d \nk=%d \nstride=%d\n", mzTopN, k,stride);
        uint16_t *y = all + gpuindexlist[k] * PeakNum;

        // todo: the algorithm can be improved!
        uint16_t used[51]={0};
        for (int i = 0; i < mzTopN; i++) {
            if(used[vecform[y[i]]]){
#ifdef _DEBUG_INDEX_WITH_PRINT_
               if(debug and debug_queryindex == gpuindexlist[k]) printf("used peak: %d\n", vecform[y[i]]);
#endif

            }
            else{
                if(vecform[y[i]]) used[vecform[y[i]]] = 1;
                s[k] += vecform[y[i]] * (PeakNum - i);

#ifdef _DEBUG_INDEX_WITH_PRINT_
                if(debug_queryindex == gpuindexlist[k] and debug and vecform[y[i]]>0)
                {
                    printf("%d -th mz: %d" ,i+1,y[i]);
                    printf("score : %d x %d -> %d\n", vecform[y[i]], PeakNum-i, s[k]);

                }
#endif
            }

        }
    }


}

__global__
void dp(int specnum, uint16_t *all, int PeakNum, uint16_t *s, int query, int mzTopN, int tol) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    s[index] = 0;
    for (int k = index; k < specnum; k += stride) {
        uint16_t *x = all + PeakNum * query, *y = all + index * PeakNum;
        uint64_t rec = 0;

        // todo: the algorithm can be improved!
        for (int i = 0; i < mzTopN; i++) {
            if (x[i] == 0) break;
            for (int j = 0; j < mzTopN; j++) {
                if (y[j] == 0) { break; }

                if (rec & (1L << j)) {
                    continue;
                }
                if (x[i] > y[j])  {
                    if (x[i] - y[j] < tol) {
                        s[index] += (PeakNum - i) * (PeakNum - j);
                        rec = (rec | (1L << j));
                        break;
                    }
                } else {
                    if (y[j] - x[i] < tol) {
                        s[index] += (PeakNum - i) * (PeakNum - j);
                        rec = (rec | (1L << j));
                        break;
                    }
                }
            }

        }
    }

}

// CUDA kernel to add elements of two arrays
//__global__
//void dpsorted(int specnum, uint16_t *all, int PeakNum, uint16_t *s, int query, int mzTopN, int tol)
//{
//    int index = blockIdx.x * blockDim.x + threadIdx.x;
//    int stride = blockDim.x * gridDim.x;
//
//
//    // preprocessing of y
////        vector<int> index(mzTopN);
////        iota(index.begin(), index.end(),0);
////        sort(index.begin(), index.end(), [](const int &a, const int &b){return y[a] < y[b];});
//
//
////        vector<int> recordj(mzTopN,0);
//    s[index] = 0;
//    for(int k = index; k < specnum; k += stride)
//    {
//        uint16_t *x = all + PeakNum * query, *y=all + index*PeakNum;
////        uint64_t rec = 0;
//
//
//
//        //---------------------
//        vector<int> yi(mzTopN), xi(mzTopN);
//        iota(yi.begin(), yi.end(), 0);
//        sort(yi.begin(), yi.end(), [y](const int &a, const int &b) { return y[a] < y[b]; });
//        iota(xi.begin(), xi.end(), 0);
//        sort(xi.begin(), xi.end(), [x](const int &a, const int &b) { return x[a] < x[b]; });
//        int i = 0, j = 0;
//
//        uint64_t recX = 0, recY = 0;
//
//
//        while (i < mzTopN and j < mzTopN) {
//            uint16_t alpha = x[xi[i]];
//            uint16_t beta = y[yi[j]];
//            if (alpha == 0) {
//                i++;
//                continue;
//            }
//            if (beta == 0) {
//                j++;
//                continue;
//            }
//            // we have minimal of x and minimal of y
//            if (alpha < beta) //x is smaller
//            {
//                if (beta - alpha <= tol) {
//                    // beta in alpha tolerance
//                    // Match!!
//                    // ------alpha-tol----------alpha-----------alpha+tol---------
//                    //----------------|-----------*-----^--------|-----------------
//                    //----------------------------------|<-beta-------------------
//                    // Keep looking for bigger one! ?? NO;  we are not sure whetherr alpha is the largest one!
//                    s[index] += (PeakNum - xi[i])*(PeakNum-yi[j]);//prescore[xi[i]][yi[j]];
//                    j++;
//                    i++;
////                    cout << score << " a=" << alpha << " b=" << beta << " xi="
////                        << xi[i] << " yi=" << yi[j] << " i="  << i <<" j=" << j << endl;
//                } else {
//                    i++; // beta is out of windows of [alpha +/- tol]
//                }
//            } else// beta<=alpha
//            {
//                if (alpha - beta <= tol) {
//                    // Beta is in windows, but less than alpha
//                    // ------alpha-tol----------alpha-----------alpha+tol---------
//                    //----------------|---^--------*--------------|-----------------
//                    //--------------------|<-beta-------------------
//                    s[index] += (PeakNum - xi[i])*(PeakNum-yi[j]);//prescore[xi[i]][yi[j]];
//                    j++;
//                    i++;
////                    cout << score << " a=" << alpha << " b=" << beta << " xi="
////                         << xi[i] << " yi=" << yi[j] << " i="  << i <<" j=" << j << endl;
//                } else {
//                    j++; // alpha is too big, beta need to catch up;
//                }
//            }
//        }
//
//
//    }
//
//
//
//}


// CUDA kernel to add elements of two arrays
__global__
void init_s(int n, uint16_t *s) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
        s[i] = 0;
}

void run_cuda_dp_example() {
    // load data from file
    string mzfile = "/data/wulong/data/honeybee/all_mzXML.txt.mz";
    const int PeakNum = 50;
    long filebytes = 0;
    File::getfilesize(mzfile, filebytes);
    long peaknum = filebytes / 2;
    long specnum = peaknum / PeakNum;
    cout << "filebytes: " << filebytes << " peaknum " << peaknum << " specnum " << specnum << endl;

    // replace with new method
    uint16_t *all;
    uint16_t *s;
    cudaMallocManaged(&all, peaknum * sizeof(uint16_t));
    cudaMallocManaged(&s, specnum * sizeof(uint16_t));


    ifstream fin(mzfile.c_str(), ios::in | ios::binary);
    fin.read((char *) all, filebytes);
    fin.close();

    int mzTopN = PeakNum;
    int tol = 15;


    // copy data to gpu
    // let's move data to gpu after init
    // Prefetch the data to the GPU
    int device = -1;
    cudaGetDevice(&device);
    cudaMemPrefetchAsync(all, peaknum * sizeof(uint16_t), device, NULL);
    cudaMemPrefetchAsync(s, specnum * sizeof(uint16_t), device, NULL);
    int blockSize = 256;
    int numBlocks = (specnum + blockSize - 1) / blockSize;

    int query = 1312081;

    blockSize = 32;
    numBlocks = (specnum + blockSize - 1) / blockSize;
    for (int i = 1; i < 2; i++) {
        dp << < numBlocks, blockSize >> > (specnum, all, PeakNum, s, query, mzTopN, tol);
    }

    cudaDeviceSynchronize();
    for (int i = 0; i < specnum and i < 1; i++) {
        cerr << s[i] << endl;
    }
    // Check for errors (all values should be 3.0f)
    float maxError = 0.0f;
    for (int i = 0; i < specnum; i++)
        maxError = fmax(maxError, fabs(s[i] - 0));
    std::cout << "Max error: " << maxError << std::endl;

    cudaFree(all);
    cudaFree(s);
}

CUDAScore::CUDAScore(string mzfile) : m_PeakNumPerSpec(50), m_filebytes(0), m_all(nullptr), m_s(nullptr),
                                      m_mzfile(mzfile) {
    cout << "Loading mz file to GPU: " << mzfile << endl;
    File::getfilesize(m_mzfile, m_filebytes);
    m_peaknum = m_filebytes / 2;
    m_specnum = m_peaknum / m_PeakNumPerSpec;

    cudaMallocManaged(&m_all, m_peaknum * sizeof(uint16_t));
    cudaCheckErrors("Fail to allocate unified memory.");
    cout << "GPU RAM USED: MB: " << m_filebytes/1024/1024 << " #peak " << m_peaknum << " #spec " << m_specnum << endl;

    cudaMallocManaged(&m_s, m_specnum * sizeof(uint16_t));
    cudaCheckErrors("Fail to allocate unified memory.");
    cout << "GPU RAM USED for store score of all: " << sizeof(uint16_t) * m_specnum/1024/1024 << "MB" << endl;
    gpuindexlist = nullptr;
    m_max_query_index_size = 500000;
    cudaMallocManaged(&gpuindexlist, (m_max_query_index_size) * sizeof(long));
    cudaCheckErrors("Fail to allocate unified memory.");
    cout << "GPU RAM USED for index list: " << sizeof(long)*(m_max_query_index_size)/1024/1024 << "MB" << endl;

    int len = 1 + UINT16_MAX;
    cout << "Building vector of size  " << len << endl;
    cudaMallocManaged(&vecform, (len) * sizeof(uint16_t));
    cudaCheckErrors("Fail to allocate unified memory.");

    ifstream fin(m_mzfile.c_str(), ios::in | ios::binary);
    fin.read((char *) m_all, m_filebytes);
    fin.close();

    int mzTopN = m_PeakNumPerSpec;
    int tol = 15;


    m_device = -1;

    cudaGetDevice(&m_device);
    cudaCheckErrors("Fail to get CUDA device.");

    cudaMemPrefetchAsync(m_all, m_peaknum * sizeof(uint16_t), m_device, NULL);
    cudaCheckErrors("Fetch cudaMem fails");
    cudaMemPrefetchAsync(m_s, m_specnum * sizeof(uint16_t), m_device, NULL);
    cudaCheckErrors("Fetch cudaMem fails");
    cout << "--- GPU is ready! ----" << endl << endl;
}



CUDAScore::~CUDAScore() {
    cout << "Started to release GPU..." << endl;
    cudaFree(m_all);
    cudaFree(m_s);
    cudaFree(vecform);
    cudaFree(gpuindexlist);
    cout << "GPU released" << endl << endl;
}


void CUDAScore::Query(int topn, int tol, int queryindex, int blockSize) {

    cout << "Start query" << endl;
    int numBlocks = (m_specnum + blockSize - 1) / blockSize;

    dp_norec << < numBlocks, blockSize >> > (m_specnum, m_all, m_PeakNumPerSpec, m_s, queryindex, topn, tol);


    cudaDeviceSynchronize();
    cout << "End of query" << endl;
}

void CUDAScore::scoreAllVecForm(int mzTopN, int tol, long queryindex, int blockSize, bool normalize) {
    int len = 1 + UINT16_MAX;
    uint16_t *queryspec = getSpecBy(queryindex);

    for (int i = 0; i < len; i++) {
        vecform[i] = 0;
    }
    for (int i = 0; i < m_PeakNumPerSpec; i++) {
        if (queryspec[i] == 0) break;
        int j = queryspec[i] >= tol ? queryspec[i] - tol : 0;
        int max_j = queryspec[i] >= len - tol ? len - tol : queryspec[i] + tol;
        while (j < max_j) {
            if (vecform[j] == 0 and j != 0) vecform[j] = (m_PeakNumPerSpec - i);
            j++;
        }
    }
    cudaMemPrefetchAsync(vecform, len * sizeof(uint16_t), m_device, NULL);
    cudaCheckErrors("Fail in cudaMem");

    int numBlocks = (m_specnum + blockSize - 1) / blockSize;
    dp_norec_vecform << < numBlocks, blockSize >> >
                                         (m_specnum, m_all, m_PeakNumPerSpec, m_s, mzTopN, tol, vecform, true, 802792);

    cudaDeviceSynchronize();

    if(normalize)
    {
        normalizationAllScore(queryindex);
    }
}


void CUDAScore::calcDotProduct(int TopNPeak, int tol, uint16_t *queryspec, int blockSize,
                               vector<long> &indexlist, vector<int> & scores){

	QueryfastOnIndexWithUINT16(TopNPeak, tol, queryspec, blockSize, indexlist);
	copyScore(indexlist, scores);
}
#include <mutex>
std::mutex dp_lock_gpu;
void CUDAScore::QueryfastOnIndexWithUINT16(int TopNPeak, int tol, uint16_t *queryspec, int blockSize,
                                           vector<long> &indexlist) {
    std::lock_guard<std::mutex> guard(dp_lock_gpu);
    int len = 1 + UINT16_MAX;

    toVectorForm(queryspec,tol,vecform, len);
    int indexlistsize = initGpuIndexList(indexlist);
    cudaMemPrefetchAsync(vecform, len * sizeof(uint16_t), m_device, NULL);
    int numBlocks = (indexlistsize + blockSize - 1) / blockSize;
    long debug_index = 19266;
    debug_index = 658620;
    dp_norec_vecform_on_index << < numBlocks, blockSize >> >
    (m_specnum, m_all, m_PeakNumPerSpec, m_s, TopNPeak, tol, vecform, gpuindexlist, indexlistsize,false,debug_index);

    cudaDeviceSynchronize();
}

void CUDAScore::queryFastOnIndex(int TopNPeak, int tol, long queryindex, int blockSize, vector<long> &indexlist) {
    uint16_t *queryspec = getSpecBy(queryindex);
    QueryfastOnIndexWithUINT16(TopNPeak, tol, queryspec, blockSize, indexlist);

}

uint16_t *CUDAScore::getscore() {
    return m_s;
}

long CUDAScore::getSpecNum() {
    return m_specnum;
}

int CUDAScore::getPeakNumPerSpec() const {
    return m_PeakNumPerSpec;
}

uint16_t *CUDAScore::getSpecBy(long queryindex) {
    return m_all + queryindex * getPeakNumPerSpec();
}

// todo: make this availabel on CPU also: Nov 16 2019
// Why we are still using this function?
vector<int> CUDAScore::distributionAll(int tol, long queryindex, bool normalize){
    int blockSize = 32;
    scoreAllVecForm(50, tol, queryindex, blockSize, normalize);
    const int MAX_TOP50_COS = 42925;
    const int MAX_SCORE = MAX_TOP50_COS + 1;
    vector<int> histogram(MAX_SCORE,0);
    for(long i = 0; i < m_specnum; i ++)    {
        int key = m_s[i];
        if(key < 0 or key > MAX_SCORE)        {
            cout << "Building histogram: invalid score: s[" << i <<  "]" << key << endl;
        } else{
            histogram[key] ++;
        }
    }
    return histogram;
}

void CUDAScore::normalizationAllScore(long queryindex) {
    const int MAX_SCORE=42925;
    double querynorm = getSquaredNorm(queryindex);

    for(long i = 0; i < m_specnum; i ++)    {
        double s = m_s[i];
        if(s>EPSILON)        {
            s/=sqrt(getSquaredNorm(i));
            s/=sqrt(querynorm);
        }

        s *= MAX_SCORE;
        if(s<0 or s>MAX_SCORE or i==802792)        {
            cout << "Normalize All scores: Invalid score s: " << s << " on index: " << i << endl;
            cout << "origin: " << m_s[i] << " norm: " << querynorm << " & " << getSquaredNorm(i) << endl;
        }
        m_s[i] = uint16_t(s)>MAX_SCORE? MAX_SCORE: uint16_t(s);
    }
}


int CUDAScore::initGpuIndexList(vector<long> &indexlist) {
    int indexlistsize = indexlist.size();
    if (indexlistsize >= m_max_query_index_size) {
        indexlistsize = m_max_query_index_size;
        cout << "[Warning] Index list exceed largest number allowed: " << m_max_query_index_size << ". Program will resize to max value" << endl;
    }
    for (int i = 0; i < m_max_query_index_size; i++) {
        if (i >= indexlistsize) {
            gpuindexlist[i] = -1;
        } else {

            gpuindexlist[i] = indexlist[i];  // todo: to be fixed long to int;
        }
    }
    cudaMemPrefetchAsync(gpuindexlist, m_max_query_index_size * sizeof(long), m_device, NULL);  // wrong!!! fixed

    return indexlistsize;
}

void CUDAScore::toVectorForm(uint16_t *queryspec, int tol, uint16_t *vecform, int len) {
    const int PeakPerSpec = getPeakNumPerSpec();
    for (int i = 0; i < len; i++) {
        vecform[i] = 0;
    }
    for (int i = 0; i < PeakPerSpec; i++) {
        if (queryspec[i] == 0) break;
        int j = queryspec[i] >= tol ? queryspec[i] - tol : 0;
        int max_j = queryspec[i] >= len - tol ? len - tol : queryspec[i] + tol;
        while (j <= max_j) {  // fix this on GPU, Now CPU and GPU are the same
            if (vecform[j] == 0 and j != 0) vecform[j] = (PeakPerSpec - i);
            j++;
        }
    }
}

__global__
void pair_dp(long *index_pairs,uint16_t * query_spec_unimem, uint16_t * all, int index_pairs_len,int query_num,
        int PeakNum, int mzTopN, int tol, uint16_t *s, bool debug, long first_idx, long second_idx){
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x; // todo: what is this?

    for (int k = index; k < index_pairs_len; k += stride) {
        s[k] = 0;
        long Idx = index_pairs[2*k], Idy = index_pairs[2*k+1];
        if(debug and Idx == first_idx and Idy == second_idx){
            printf("Score---\n");
        }
        uint16_t *y = query_spec_unimem + PeakNum * Idx, *x = all + Idy * PeakNum;
        uint64_t rec = 0;

        for (int i = 0; i < mzTopN; i++) {
            if (x[i] == 0) break;
            for (int j = 0; j < mzTopN; j++) {
                if (y[j] == 0) { break; }

                if (rec & (1L << j)) {
                    continue;
                }
                if (x[i] > y[j])
                {
                    if (x[i] - y[j] <= tol) {

                        s[k] += (PeakNum - i) * (PeakNum - j);

                        if(debug and Idx == first_idx and Idy == second_idx){
                            printf("%ld  %ld  %d <- %d x %d mz: %d vs %d\n",Idx, Idy, s[k],PeakNum-i,PeakNum-j,x[i],y[j]);
                        }

                        rec = (rec | (1L << j));

                        break;
                    }
                } else {
                    if (y[j] - x[i] <= tol) {

                        s[k] += (PeakNum - i) * (PeakNum - j);

                        if(debug and Idx == first_idx and Idy == second_idx){
                            printf("%ld  %ld  %d <- %d x %d mz: %d vs %d\n",Idx, Idy, s[k],PeakNum-i,PeakNum-j,x[i],y[j]);
                        }

                        rec = (rec | (1L << j));

                        break;
                    }
                }
            }

        }
    }
}

// Case 1: one query
// Case 2: many queries
void CUDAScore::dpscore(double tolerance, vector<vector<long>> &allRetIdx, int threadnum, vector<vector<float>> &accDist,
                   ICQuery &query,vector<vector<int>> &dpscores) {
    bool debug = false;
    long first = 5;
    long debug_index = 984171;

    // step 1. memory for queries
    uint16_t  *querySpecGpu;
    cudaMallocManaged(&querySpecGpu, query.size()*getPeakNumPerSpec() * sizeof(uint16_t));
    long * index_pairs;
    int num = 0;
    for(auto x: allRetIdx) {
//        cout << "number of ret idx: " << x.size() << endl;
        num += x.size();
    }

    cudaMallocManaged(&index_pairs, 2*num*sizeof(long));
    uint16_t *dps;
    cudaMallocManaged(&dps,num*sizeof(uint16_t));
    //printf("allocated on GPU: querySpecGpu %ld; index_pairs: %ld; dps: %ld byts.\n",query.size()*getPeakNumPerSpec() * sizeof(uint16_t),2*num*sizeof(long),num*sizeof(uint16_t));

    // step 2. copy queryies to gpu
    for(int i = 0; i < query.size(); i ++)   {
        uint16_t  *p  = query.getPtrUint16(i);
        copy(p, p+getPeakNumPerSpec(),querySpecGpu + i*getPeakNumPerSpec() );
    }

    int k = 0;
    for(int i = 0; i < allRetIdx.size(); i ++)  {
        long first = i;
        for(int j = 0; j < allRetIdx[i].size();  j++)
        {
            long second = allRetIdx[i][j];
            index_pairs[2*k] = first;
            index_pairs[2*k+1] = second;
            k++;
        }
    }

    // step 3. calculate dp score for each pair of them
    int useTopN=50;
    int tol = 1;
    pair_dp<<<32,256>>>(index_pairs,querySpecGpu,m_all,num,query.size(),m_PeakNumPerSpec,useTopN,tol,dps, debug,first,debug_index);

    // step 4. copy score back
    cudaDeviceSynchronize();

    // step 5. return
    cudaFree(querySpecGpu);


    k = 0;
    for(int i = 0; i <query.size(); i ++)  {
        uint16_t  * queryspec = query.getPtrUint16(i);
        vector<long> &indexlist=allRetIdx[i];
        float querySquaredNorm = getSquaredNorm(queryspec);
        vector<float> accDist1st(indexlist.size(), 0);
        dpscores[i].assign(indexlist.size(),0);
        for (int j = 0; j < indexlist.size(); j++) {
            float candSquaredNorm = getSquaredNorm(indexlist[j]);
            accDist1st[j] = dps[k];
            dpscores[i][j] = dps[k];
            k++;
            if (candSquaredNorm > EPSILON and querySquaredNorm > EPSILON) {
                accDist1st[j] /= sqrt(candSquaredNorm * querySquaredNorm);
            }
            accDist1st[j] = sqrt(2.0) * sqrt(1.0 - accDist1st[j] < EPSILON ? 0 : 1.0 - accDist1st[j]);
        }
        accDist[i].swap(accDist1st);
    }

    cudaFree(dps);
    cudaFree(index_pairs);
}

void CUDAScore::copyScore(vector<long> &indexlist, vector<int> &scores) {
    int scoreLen = indexlist.size() > m_max_query_index_size? m_max_query_index_size: indexlist.size();
    scores.assign(scoreLen,0);
    for (int i = 0; i < scores.size(); i ++)    {
        scores[i] = getscore()[i];
    }
}

