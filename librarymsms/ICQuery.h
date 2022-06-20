//
// Created by wulong on 2/22/19.
//

#ifndef MYTOOL_ICQUERY_H
#define MYTOOL_ICQUERY_H


#include <memory>
#include <vector>
using namespace std;
class CMzSpec;

class ICQuery  {
public:
    virtual ~ICQuery();
    virtual void get(vector<float> &queries)=0;
    virtual float* get()=0;
    virtual long dim()=0;
    virtual int size()=0;
    virtual uint16_t * getPtrUint16(int i) = 0;
    virtual long getQueryIndex(int i) = 0;
    ICQuery& L2Normalization();
	void print(int i);
	long getSpecNum() {return size();}
    CMzSpec getMzSpec(long queryindex);
    int getPeakNumPerSpec() const {return 50;}
    uint16_t * getSpecBy(long queryindex) {return getPtrUint16(queryindex);}
    void QueryfastOnIndexWithUINT16(int TopNPeak, int tol, uint16_t *queryspec, int blockSize,
                                    vector<long> &indexlist, vector<int> &scores);
	void conv(int i, int j, int tol = 1);
    void conv(uint16_t *a, uint16_t * b, int tol=1) const;
};

#endif //MYTOOL_ICQUERY_H
