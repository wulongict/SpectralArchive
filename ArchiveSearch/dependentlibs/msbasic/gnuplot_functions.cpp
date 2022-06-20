//
// Created by wulong on 11/25/17.
//

#include "gnuplot_functions.h"
#include "../../../External/gnuplot-iostream/gnuplot-iostream.h"

void scatterplot(const string title, const vector<vector<double>> &fc_score) {
    Gnuplot gnuplot;
    gnuplot << "unset key" << endl;
    gnuplot << "set title \'" << title << "\'" << endl;
    gnuplot << "set xlabel ' intensity fold ' " << endl;
    gnuplot << "set ylabel ' score' " << endl;
    gnuplot << "plot '-' with points " << endl;
    gnuplot.send1d(fc_score);
}

void gnuplot_histogram(const string &title, const vector<double> &score_plot) {
    Gnuplot gnuplot;
    gnuplot << "unset key" << endl;
    gnuplot << "set title \'" << title << "\'" << endl;
    gnuplot << "set xlabel ' score ' " << endl;
    gnuplot << "set ylabel ' frequency ' " << endl;
    gnuplot << "set style histogram" << endl;
    gnuplot << "    set boxwidth 0.009 absolute " << endl;
    gnuplot << "set style fill solid 1.0 noborder " << endl;

    gnuplot << "bin_width = 0.01; " << endl;
    gnuplot << "bin_number(x) = floor(x/bin_width)" << endl;
    gnuplot << "rounded(x) = bin_width * ( bin_number(x) + 0.5 )" << endl;
    gnuplot << "plot " << gnuplot.file1d(score_plot) << " using (rounded($1)):(1) smooth frequency with boxes " << endl;
}