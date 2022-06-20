//
// Created by wulong on 1/31/19.
//

#ifndef MYTOOL_VISUAL_H
#define MYTOOL_VISUAL_H

#include <string>
#include <vector>
#include <algorithm>
#include <vector>
#include <memory>
using namespace std;
namespace gnuplotio
{
    class Gnuplot;
}

namespace CVisual
{
    struct SImageSize    {
        int m_dpi ;
        int m_width ; // inch
        int m_height ; // inch
        int m_wpx;
        int m_hpx;

        SImageSize(int width, int height, int dpi);
        SImageSize();
        void calcPixel();
        void init();
    };

    class gnuplotWrapper{
        string xlabel;
        string ylabel;
        string outputfilename;
        string terminal_type;
        int width;
        int height;

        shared_ptr<gnuplotio::Gnuplot> m_gnuplot;
        double minval;
        double maxval;

        double yminval;
        double ymaxval;
    public:
        gnuplotWrapper();
        gnuplotWrapper & set_xlabel(string xlabel_);
        gnuplotWrapper & set_ylabel(string ylabel_);
        gnuplotWrapper & set_filename(string filename_);
        gnuplotWrapper & set_terminaltype(string terminaltype_);
        gnuplotWrapper & set_height(int height_);
        gnuplotWrapper & set_width(int width_);
        gnuplotWrapper & set_minmax(double minval_, double maxval_);
        gnuplotWrapper & set_yminmax(double yminval_, double ymaxval_);
        void setLabel(gnuplotio::Gnuplot &gnuplot);
        void setTerminal(gnuplotio::Gnuplot &gnuplot);
        void setMinMaxRange(gnuplotio::Gnuplot &gnuplot);
        void setYMinMaxRange(gnuplotio::Gnuplot &gnuplot);
        void setOutput(gnuplotio::Gnuplot &gnuplot);
        // This function is never used!
        void config(gnuplotio::Gnuplot &gnuplot);
        void config();
        void pngparams();
        void plotLinearRegression(vector<double> &coefs, double rSquared, const vector<vector<double>> &rx_logry) const;
    };
    namespace impl{
        void gnuplot_histogram_impl(const vector<double> &s, int binNum, gnuplotio::Gnuplot &gnuplot, bool replot=false,  double minVal=0, double maxVal=-1);

        void gnuplot_curves_impl(const vector<vector<tuple<double, double>>> &curves, gnuplotio::Gnuplot &gnuplot, bool replot=false);

    }
    void gnuplot_histogram(const vector<double> &s, int binNum=-1, double minval=0, double maxval=-1, int width=100,int heigth=40);
    void gnuplot_histogram_topng(const vector<double> &s, int binNum, double minVal, double maxVal,
                                 gnuplotWrapper &info);
    void gnuplot_curve(const vector<double> &x, const vector<double> &y);
    void gnuplot_curve_topng(const vector<double> &x, const vector<double> &y, string type, gnuplotWrapper &info);
    void gnuplot_curves(const vector<vector<tuple<double, double>>> &curves);
    void gnuplot_curves(const vector<vector<tuple<double, double>>> &curves, string pngFilename);
    void gnuplotPeptideProphetFigure(const vector<vector<tuple<double, double>>> &curves, const vector<double> &hist, string pngFilename);




};

#endif //MYTOOL_VISUAL_H
