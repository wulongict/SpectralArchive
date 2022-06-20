//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_GNUPLOT_FUNCTIONS_H
#define MYTOOL_GNUPLOT_FUNCTIONS_H



#include <string>
#include <vector>
using namespace std;

void gnuplot_histogram(const string &title, const vector<double> &score_plot);
void scatterplot(const string title, const vector<vector<double>> &fc_score);
//void gnuplot_histogram(const string &testingfile, const vector<double> &score_plot);

#endif //MYTOOL_GNUPLOT_FUNCTIONS_H
