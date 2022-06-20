//
// Created by wulong on 12/9/19.
//

#ifndef MYTOOL_ICHISTDATA_H
#define MYTOOL_ICHISTDATA_H


#include <vector>
#include <iostream>

using namespace std;

class ICHistdata {
    // Abstract class of all the histogram
public:
    virtual int size() = 0;
    bool empty(){return 0==size();}
    virtual double at(int i) = 0;
    virtual double accumulate(int shift)=0;
//    virtual double nonzerosum(double leftLim)=0;
    virtual double gettotal()=0;
    virtual vector<double> getX()=0;
    virtual vector<double> getSurvivalFn(bool nonZero)=0;
};

// can be improved!
class CHistInt: public ICHistdata{
    // This class is designed for score of dot product with 50 peaks.
    // too limited!!!
    vector<int> m_counts;           // The frequency
    vector<double> m_left_boundary; // The left boundary of the bin/bucket
    vector<double> m_cdf;
    double m_total;
public:
    explicit CHistInt(vector<int> allscores);
    CHistInt(vector<int> allscores, vector<double> left_boundary);
    vector<double> getSurvivalFn(bool nonZero) override;

    vector<double> getX() override;
    int size() override {return m_counts.size();}
    double at(int i) override{return m_counts.at(i);}
    double accumulate(int shift) override;
    double gettotal() override{return m_total;}
private:
    // add up frequency of bins with positive left bounary.
    double nonzerosum(double leftLim) ;


};


#endif //MYTOOL_ICHISTDATA_H
