//
// Created by wulong on 6/19/19.
//

#ifndef MYTOOL_CFAISSINDEX_H
#define MYTOOL_CFAISSINDEX_H

#include "ICIndexWrapper.h"
#include <faiss/Index.h>
#include <faiss/impl/AuxIndexStructures.h>

class CFaissIndexWrapper: public ICIndexWrapper {
    faiss::Index *m_index;
    // faiss::Index * m_tmp;
    bool m_isCPU;
    vector<int> m_gpu_idx;
public:
    CFaissIndexWrapper();
    void removeIds(vector<long> &idx) override;
    virtual ~CFaissIndexWrapper();
    void read(string filename) override ;
    void add(long newspecnum, float *vec) override;
    void display() override;
    void write(string filename) override;
    long total() override;
    void search(int num, float * queries, int ret_num, vector<float>& dist, vector<long>& ind);
    void setnProb(int nprobe) override;
    bool empty() override;
    int dim() override;
    bool istrained() override;
    void toCPU() override;
    void createEmptyIndex(int dim, string indexstr) override;
    void toGPU() override; // todo
    void toMultipleGPUs(vector<int> &gpu_idx) override; // todo
    void train(int batchSize, float * vecBatch) override;
private:
    void setPtr(faiss::Index *newptr);
    // void backup();
    faiss::Index * getPtr();
};


#endif //MYTOOL_CFAISSINDEX_H
