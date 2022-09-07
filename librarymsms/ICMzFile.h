//
// Created by wulong on 2/28/19.
//

#ifndef MYTOOL_ICMZFILE_H
#define MYTOOL_ICMZFILE_H

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include "ICHistdata.h"
using namespace std;
class ICQuery;
class CMzSpec;
class ICHistdata;

class ICMzFile
{
public:
    virtual string getClassName();
    static vector<int> score_to_histogram(vector<int> &scores);

    // Create Histogram with list of score and bin_num
    // bin_num: 10,000
    //
    template<typename T>
    static shared_ptr<ICHistdata> score2histogram(vector<T> &scores, int bin_num){
        // find max and min value from scores
        auto minmax_score = minmax_element(scores.begin(), scores.end());
        double minval = *minmax_score.first, maxval = *minmax_score.second;
        bool verbose = false;
        if(verbose)cout << "max, min " << maxval << "\t" << minval << endl;

        // divide the histogram into 10,001 buckets. and create its boundaries
        vector<int> histogram(bin_num+1,0);
        vector<double> left_boundary(bin_num+1, 0);

        // the step is based on the max and minimum value...
        double step = (maxval - minval) / bin_num;
        if(fabs(step)<1e-10){
            if(verbose){
                cout << "max, min " << maxval << "\t" << minval << endl;
                cout << "no valid dot product found will return a zero histogram (nullptr)" << endl;
            }
            return nullptr;
        }
        // the last left boundary is the max value
        for(int i = 0; i <= bin_num; i ++)  {
            left_boundary[i] = minval + i * step;
        }

        // create the histogram...
        for (int i = 0; i < scores.size(); i++) {
            int k = floor((scores[i]-minval)/step);
            histogram[k] ++;
        }
        return make_shared<CHistInt>(histogram,left_boundary);// 20200806
    }
//    static vector<int> score_to_histogram(vector<double> &scores, int bin_num);
    static float getSquaredNorm(const uint16_t *p, const int PEAKNUM_PER_SEPC) ;
    virtual ~ICMzFile(){}
    virtual long getSpecNum()=0;
    virtual uint16_t * getSpecBy(long queryindex)=0;
    virtual int getPeakNumPerSpec() const = 0;

    // calculating dot product between queryspec and spectra listed in indexlist
    // implementd in subclasses.
    // The GPU version has an upper limit on the size of indexlist, 500,000
    // The CPU version has no limit.
    virtual void calcDotProduct(int TopNPeak, int tol, uint16_t *queryspec, int blockSize,
                                vector<long> &indexlist, vector<int> &scores) = 0;
    virtual void dpscore(double tolerance, vector<vector<long>> &allRetIdx, int threadnum, vector<vector<float>> &accDist,
                          ICQuery &query,vector<vector<int>> &dpscores)=0;

    virtual int getPeakNum(long queryindex);

    virtual float getSquaredNorm(long queryindex);
    virtual float getSquaredNorm(uint16_t *p);
    void get_vector_form(long idx, int tol, vector<int> &vecform);  // only used once
	void get_vector_form(uint16_t * x, int tol, vector<int> &vecform) const;
    CMzSpec getMzSpec(long queryindex);


    float *buildNormalizedQueries(int dim, const vector<long> &candidates, bool useFlankingBins, int numQuery);

    void getQueryVecByMZArrayIndex(int queryindex, vector<float> & v, bool useFlankingBins=false) ;

    // Several different implementation of dot product calcualtion.
    long calculate_dot_product(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore);

    long calculate_dot_product_with_vecfrom(long queryX, vector<int> &vecformY, int mzTopN, bool debug,
                                            long debug_index);
    long calculate_dot_product_with_hist(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore, vector<double> &diff);
    long calculate_dot_product_with_hist(int tol, int mzTopN, vector<vector<int>> &prescore, vector<double> &diff,
                                         uint16_t *x, uint16_t *y, bool print_matched_pks) const;
    long calculate_dot_product_norec(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore);
    long calculate_dot_product_norec_short_if_else(long queryX, long queryY, int tol, int mzTopN,
                                                  vector<vector<int>> &prescore);
    long calculate_dot_product_sorted(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore);
    void get_prescoreMatrix(vector<vector<int>> &prescorematrix) const;
    void distOnSpec(uint16_t *queryspec, vector<long> &ind, int tol, vector<float> &dist,vector<int> &dpscore);
    void dist(long query_index, vector<long> &ind, int tol, vector<float> &distances,vector<int> &dpscore);

    static void displayspec(vector<uint16_t> &spec);

	virtual void scorePartiallyWithVecForm(int mzTopN, int tol, int blockSize, bool normalize, vector<long> &indexlist,
                                                      uint16_t  *queryspec, vector<int> &scores);

    virtual shared_ptr<ICHistdata> distributionPartialGeneral(int tol, bool normalize, vector<long> &indexlist, uint16_t *queryspec);
	virtual vector<int> distributionPartial(int tol, bool normalize, vector<long> &indexlist, uint16_t *queryspec);
    virtual vector<int> distributionAll(int tol, long queryindex, bool normalize);
    void getCompactForm(const double *mz, const double *intensity, vector<uint16_t> &newspec) const;
protected:

    void getsortindex(vector<int> &idx, const double *intensity) const;
};

shared_ptr<ICMzFile> createScorer(bool use_gpu,string mzXMLListFileName);

class ICMzFactory{
    string mzfilename;
    shared_ptr<ICMzFile> m_csa;
public:
    void setfilename(string mzfile){
        mzfilename = mzfile;
    }
    void setMzFileObject(shared_ptr<ICMzFile> csa){
        m_csa = csa;
    }
    shared_ptr<ICMzFile> getMzFileObject();
    string getmzfilename();
    virtual ~ICMzFactory(){}
    virtual shared_ptr<ICMzFile> create()=0;
};



#endif //MYTOOL_ICMZFILE_H
