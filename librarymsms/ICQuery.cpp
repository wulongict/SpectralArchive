//
// Created by wulong on 2/22/19.
//
#include <iostream>
#include <memory>
#include <algorithm>
#include <vector>
#include "ICQuery.h"
#include "ProteomicsDataTypes.h"
#include <cmath>
//#include "PeakList.h"
#include "Util.h"

ICQuery::~ICQuery(){
}

ICQuery& ICQuery::L2Normalization(){
	float *p = get();
	long d = dim();
	for (int j = 0; j < size(); j++)	{
		float *v = p + j * d;

		// normalize number!
		double sum = 0;
		for (int i = 0; i < dim(); i++){
			sum += (v[i] * v[i]);
		}
        sum = sqrt(sum);
		if (sum > EPSILON){
            for (int i = 0; i < dim(); i++) {
                v[i] /= sum;
            }
		}
	}
	return *this;
}

void ICQuery::print(int j){
	cout << "query: " << j << endl;
	float *p = get();
	long d = dim();
	
	float *v = p + j * d;
	for (int i = 0; i < dim(); i++)	{
		if(v[i]>EPSILON)
		cout << i*0.5 << ":" << v[i] << "\t";
	}
	cout << endl;
	
}

void ICQuery::conv(uint16_t *a, uint16_t *b, int tol) const {
    const int PeakNum = getPeakNumPerSpec();
    vector<int> a_matched(PeakNum,0), b_matched(PeakNum,0);
    for(int i = 0; i < PeakNum; i ++) {
        if(a[i] == 0) break;
        for(int j = 0; j < PeakNum; j ++) {
            if(b[j] ==0) break;
            if((a[i] >= b[j] and a[i]-b[j] < tol )
               or (a[i]<=b[j] and b[j] -a[i] < tol)){
                a_matched[i] = 1;
                b_matched[j] = 1;
            }
        }
    }

    const int MAX_PEAK_MZ = UINT16_MAX;  // 65535
    const double max_mz = 2000.0;
    const double max_score = 42925.0;
    vector<float> conv_vec(MAX_PEAK_MZ*2-1);
    for(int k = 0; k < conv_vec.size(); k ++)   {
        for(int i = 0; i < PeakNum; i ++)  {
            if(a[i] == 0) break;
            if(a_matched[i]) continue;
            for(int j = 0; j < PeakNum; j ++)  {
                if(b_matched[j]) continue;
                if(b[i] == 0) break;
                if(abs(MAX_PEAK_MZ+a[i]-b[j] -k) < tol) conv_vec[k] += (PeakNum-i) * (PeakNum-j);
//                    if(MAX_PEAK_MZ+a[i]-b[j] == k) conv_vec[k] += (PeakNum-i) * (PeakNum-j);
            }
        }
    }
    cout << "threshold: " << 3000/max_score << endl;
    for(int k = 0; k < conv_vec.size(); k ++)  {
        if(conv_vec[k] > 3000) cout << (k-MAX_PEAK_MZ)*2000.0/MAX_PEAK_MZ << " :" << conv_vec[k]/max_score << endl;
    }
}

CMzSpec ICQuery::getMzSpec(long queryindex) {
    return CMzSpec(getSpecBy(queryindex), getPeakNumPerSpec());
}

void ICQuery::conv(int i, int j, int tol) {
    conv(getPtrUint16(i), getPtrUint16(j), tol);
}

void
ICQuery::QueryfastOnIndexWithUINT16(int TopNPeak, int tol, uint16_t *queryspec, int blockSize, vector<long> &indexlist,
                                    vector<int> &scores) {
    cout << "do nothing!" << endl;
}
