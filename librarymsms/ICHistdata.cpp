//
// Created by wulong on 12/9/19.
//

#include "ICHistdata.h"
#include <memory>
#include <numeric>
#include <algorithm>
#include <utility>

double CHistInt::accumulate(int shift) {
    // empirical CDF...
    return std::accumulate(m_counts.begin() + shift, m_counts.end(), 0);
}

double CHistInt::nonzerosum(double leftLim) {
    // CDF to left lim
    double sum = 0;
    for(int i = 0; i < m_left_boundary.size(); i ++)   {
        if(m_left_boundary[i]> leftLim)    {
            sum += m_counts.at(i);
        }
    }
    return sum;
}

CHistInt::CHistInt(vector<int> allscores) : m_counts(std::move(allscores)){
    cout << "---------------------------------"
         << endl << " The code should be not use as it has assumption of score between 0 to 42925 "
         << endl <<"-----------------------------------" << endl;
    if(size()==0){throw logic_error("score should not be zero!");}
    m_left_boundary.assign(size(),0);
    const int MAX_TOP50_COS = 42925;
    double step = 1.0 / MAX_TOP50_COS;
    m_left_boundary[0] = 0;
    for (int i = 1; i < size(); i++) {
        m_left_boundary[i] = m_left_boundary[i-1] + step;
    }
    m_cdf.assign(size(),0);
}

CHistInt::CHistInt(vector<int> allscores, vector<double> left_boundary):m_counts(std::move(allscores)) {
    double left=0.0, right=1.0;
    auto it = std::find_if(left_boundary.begin(), left_boundary.end(),[&](const double x)->bool{return x < left or x > right;});
    if( it==left_boundary.end())   { // not found
        m_left_boundary = left_boundary;
    } else{
        m_left_boundary.assign(size(),0);
        const int MAX_TOP50_COS = 42925;
        double step = 1.0 / MAX_TOP50_COS;
        for (int i = 0; i < size(); i++) {
            m_left_boundary[i] = left_boundary[i] * step;
        }
    }
    m_cdf.assign(size(),0);
}

vector<double> CHistInt::getSurvivalFn(bool nonZero) {
    // this is survival function, not CDF
    double leftlim = 1e-3;
    m_cdf[size()-1]=at(size()-1);
    for (int i = size() - 2; i >= 0 and m_left_boundary.at(i)>leftlim; i--) {
        m_cdf[i] = at(i) + m_cdf.at(i+1);
        m_total=m_cdf[i];
    }
    for(auto & x: m_cdf){x/=m_total;}
    return m_cdf;
}

vector<double> CHistInt::getX() {
    return m_left_boundary;
}
