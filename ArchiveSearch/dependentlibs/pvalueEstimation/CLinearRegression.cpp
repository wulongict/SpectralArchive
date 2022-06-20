//
// Created by wulong on 4/3/19.
//

#include "CLinearRegression.h"
#include <iostream>

using namespace std;

CLinearRegression::CLinearRegression() : m_R_square{0} {
    m_A = vector<vector<double>>{vector<double>{1, 2}, vector<double>{2, 4}, vector<double>{3, 6}};
    m_b = vector<double>{5, 9.0001, 13.01};
}

// Eigen module is using multiple threads
void CLinearRegression::solve(bool verbose) {
    if (verbose)cout << m_A.size() << "----" << endl;
    if (verbose and m_A.size() > 0) { cout << m_A[0].size() << "0---size---" << endl; }
    const int rows = m_A.size();

    const int cols = m_A[0].size();

    MatrixXd A(rows, cols + 1);
    VectorXd b(rows);
    MatrixXd w(rows, rows);
    w.setIdentity();

    VectorXd x(rows);
    for (int i = 0; i < rows; i++) {
        b(i) = m_b[i];
        for (int j = 0; j < cols; j++) {
            A(i, j) = m_A[i][j];
        }
        A(i, cols) = 1;
    }

    VectorXd Ax = A.transpose() * b;
    MatrixXd AwA = A.transpose() * w * A;
    x = AwA.jacobiSvd(ComputeThinU | ComputeThinV).solve(Ax);
    if (verbose)cout << "solution is: \n" << x << endl;
    m_x.assign(cols + 1, 0);
    for (int i = 0; i < cols + 1; i++) {
        m_x[i] = x(i);
    }

    VectorXd deltaResidual = b - A * x;
    VectorXd b_bar = VectorXd::Constant(b.size(), b.mean());
    double RSS = deltaResidual.transpose() * deltaResidual;
    double TSS = (b - b_bar).transpose() * (b - b_bar);
    m_R_square = 1 - RSS / TSS;
    if (verbose)cout << "R^2 is: " << m_R_square << endl;
    if (verbose)cout << "--------------------------------" << endl;

}


