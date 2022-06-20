//
// Created by wulong on 6/19/19.
//

#ifndef MYTOOL_CMYINDEX_H
#define MYTOOL_CMYINDEX_H

#include "ICIndexWrapper.h"
class CProductQuantization;



class CMyIndex: public ICIndexWrapper {
    /*
     * This index is implemented by myself currently only work on CPU
     * */
    shared_ptr<CProductQuantization> cpq;
    shared_ptr<CPQParam> m_option;
public:
//    CMyIndex(CPQParam option);
    CMyIndex(shared_ptr<CPQParam> optionptr);
//    CMyIndex();
    virtual ~CMyIndex();
    void read(string filename);
    void add(long newspecnum, float *vec) override;
    void display() override;
    void write(string filename);
    long total() override;
    void search(int num, float * queries, int ret_num, vector<float>& dist, vector<long>& ind);
    void setnProb(int nprobe) override;
    bool empty() override;
    int dim() override;
    bool istrained() override;
    void createEmptyIndex(int dim, string indexstr) override;
    void train(int batchSize, float * vecBatch) override;
    void toGPU() override; // todo
    virtual void toCPU() override;
};


#endif //MYTOOL_CMYINDEX_H
