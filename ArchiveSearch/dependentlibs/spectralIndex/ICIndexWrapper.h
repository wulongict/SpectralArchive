//
// Created by wulong on 2/18/19.
//

#ifndef MYTOOL_ICINDEXWRAPPER_H
#define MYTOOL_ICINDEXWRAPPER_H

#include <string>
#include <vector>
#include <memory>
using namespace std;
class CPQParam;

class ICIndexWrapper
{
    string m_filename;
public:
    virtual ~ICIndexWrapper();
    void read();
    void write();
    void setfilename(string filename);
    string getfilename();

    virtual void toGPU()=0;
    virtual void toCPU()=0;

    // API: Pure Virtual functions
    virtual void read(string filename)=0;
    virtual void write(string filename)=0;
    virtual long total()=0;
    virtual bool istrained()=0;
    virtual bool empty()=0;
    virtual int dim()=0;
    virtual void createEmptyIndex(int dim, string indexstr) = 0;
    virtual void display()=0;
    virtual void setnProb(int nprobe)=0;
    virtual void train(int batchSize, float * vecBatch) = 0; // todo: after training step, we do not have the centroid object anymore!!!
    virtual void add(long newspecnum, float *vec) = 0;
    virtual void search(int num, float * queries, int ret_num, vector<float>& dist, vector<long>& ind)=0;
};


shared_ptr<ICIndexWrapper> IndexFactory(bool myIndex, shared_ptr<CPQParam> option);

#endif //MYTOOL_ICINDEXWRAPPER_H
