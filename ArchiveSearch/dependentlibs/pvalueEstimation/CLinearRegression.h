//
// Created by wulong on 4/3/19.
//

#ifndef MYTOOL_CLINEARREGRESSION_H
#define MYTOOL_CLINEARREGRESSION_H

#include "Core"
#include "Dense"
#include <iostream>
#include <vector>
#include <iterator>
#include <random>
using namespace Eigen;
using namespace std;

class CLinearRegression
{
    vector<vector<double>> m_A;
    vector<double> m_b;
    vector<double> m_x; // solution
    // model A x = b
    // solution (no weight):
    // A' * A x = A' b
    // x = inv(A'*A) A' b

    vector<double> m_w;
    double m_R_square;
    // solution (with weight)
    // x= inv(A'WA) A' b
public:
    CLinearRegression();
    CLinearRegression(const vector<vector<double>> &A, const vector<double> &b) :m_R_square{0}{
        m_A = A;
        m_b = b;
    }
    void solve(bool verbose);

    vector<double> getCoefficients(){return m_x;}
    double getRsquare(){return m_R_square;}
};

#endif //MYTOOL_CLINEARREGRESSION_H
