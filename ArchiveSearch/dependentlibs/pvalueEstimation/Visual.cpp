//
// Created by wulong on 1/31/19.
//

#include "Visual.h"
#include "../../../External/gnuplot-iostream/gnuplot-iostream.h"

void CVisual::gnuplot_histogram(const vector<double> &s, int binNum, double minval, double maxval,int width, int height) {
    cout << "[gnuplot] histogram: xrange: " << minval << " " << maxval << endl;
    cout << "[gnuplot] histogram: figure size: width  " << width << ", height "<< height << endl;
    Gnuplot gnuplot;

    gnuplot << "set term dumb size " << width << " , " <<  height << endl;
//    if(minval )
//    gnuplot << "set xrange [" << minval << ":" << maxval <<"]" << endl;
    CVisual::impl::gnuplot_histogram_impl(s, binNum*2, gnuplot, false, minval, maxval);

}

void CVisual::impl::gnuplot_histogram_impl(const vector<double> &s, int binNum, Gnuplot &gnuplot, bool replot, double minVal, double maxVal) {
    double minval = minVal, maxval = maxVal;
    cout << "[gnuplot] histogram: xrange: " << minval << "\t" << maxval  << endl;
    if(minVal>=maxVal)    {
        maxval = *max_element(s.begin(), s.end());
        minval = *min_element(s.begin(), s.end());
    }

    cout << "[gnuplot] histogram xrange (updated): " << minval << "\t" << maxval << endl;

    int N = s.size();
    if(binNum<=0)    {
        binNum =  N/20;
    }
//    int binNum = N/10;
    double binWidth = (maxval-minval)/binNum;
    gnuplot << "n=" << binNum << endl;
    gnuplot << "max=" << maxval << endl;
    gnuplot << "min=" << minval << endl;
    gnuplot << "width=(max-min)/n" << endl;
    gnuplot << "hist(x, width)=width*floor(x/width)+width/2.0" << endl;
    gnuplot << "set boxwidth width*0.9" << endl;
    gnuplot << "set style fill solid 0.5" << endl;
    gnuplot << "set xrange [" << minval << ":" << maxval << "]" << endl;

    gnuplot << "stat " << gnuplot.file1d(s) << " u 1 nooutput" << endl;
    gnuplot << "sum = STATS_records" << endl;
    if(replot)    {
        gnuplot << ", " ;
    } else{
        gnuplot << "plot " ;
    }
    gnuplot << gnuplot.file1d(s)<<"  u (hist($1,width)):(1.0/width/sum) smooth freq w boxes lc rgb\"green\" notitle" << endl;
}

void CVisual::gnuplot_curve(const vector<double> &x, const vector<double> &y) {
    Gnuplot gnuplot;
    gnuplot << "set term dumb size 100, 40" << endl;
    vector<tuple<double, double>> data;
    assert(x.size()==y.size());
    data.reserve(x.size());
    for(int i = 0; i < x.size(); i ++)    {
        data.emplace_back(make_tuple(x[i],y[i]));
    }
    gnuplot << "plot " << gnuplot.file1d(data) << " u 1:2 with lines notitle" << endl;
}

void CVisual::gnuplot_curve_topng(const vector<double> &x, const vector<double> &y, string type, gnuplotWrapper &info) {
    Gnuplot gnuplot;
    info.config(gnuplot);
    vector<tuple<double, double>> data;

    data.reserve(x.size());
    for(int i = 0; i < x.size(); i ++)    {
        data.emplace_back(make_tuple(x[i],y[i]));
    }
    gnuplot << "plot " << gnuplot.file1d(data) << " u 1:2 with " << type <<  " notitle" << endl;
}

void CVisual::gnuplot_curves(const vector<vector<tuple<double, double>>> &curves) {
    Gnuplot gnuplot;
    gnuplot << "set term dumb size 100, 40" << endl;
    CVisual::impl::gnuplot_curves_impl(curves, gnuplot);
}


void CVisual::gnuplot_curves(const vector<vector<tuple<double, double>>> &curves, string pngFilename) {
    Gnuplot gnuplot;
    gnuplot << "set term pngcairo size 1200, 900" << endl;
    gnuplot << "set output \""<< pngFilename<< "\"" << endl;
    CVisual::impl::gnuplot_curves_impl(curves, gnuplot);
}

void
CVisual::gnuplotPeptideProphetFigure(const vector<vector<tuple<double, double>>> &curves, const vector<double> &hist,
                                     string pngFilename) {
    Gnuplot gnuplot;
    gnuplot << "set term pngcairo size 1200, 900" << endl;
    gnuplot << "set output \""<< pngFilename<< "\"" << endl;
    CVisual::impl::gnuplot_histogram_impl(hist,100, gnuplot);
    cout << "Plot Second half, the curves" << endl;
    CVisual::impl::gnuplot_curves_impl(curves, gnuplot,true);
    cout << "Curves done" << endl;
}


// the distribution will be either split into two group or just merged into a large flat distribution
// this is because, when the nodes are clustered into several clusters, the distance inside cluster
// and between cluster are quite different, so the value goes to two extreme region.
// when the nodes are not that close, not form nice clusters, the nodes will have loose edges connected.
// Then we see the distribution is more neutral and flat.
void
CVisual::gnuplot_histogram_topng(const vector<double> &s, int binNum, double minVal, double maxVal,
                                 gnuplotWrapper &info) {
    Gnuplot gnuplot;
    info.setTerminal(gnuplot);
    info.setOutput(gnuplot);
    info.setLabel(gnuplot);
    info.setMinMaxRange(gnuplot);

    cout << "[gnuplot] histogram: xrange: " << minVal << "\t" << maxVal << endl;
    CVisual::impl::gnuplot_histogram_impl(s, binNum, gnuplot,false, minVal, maxVal);

}

void CVisual::impl::gnuplot_curves_impl(const vector<vector<tuple<double, double>>> &curves, Gnuplot &gnuplot, bool replot) {
    if(replot){
        gnuplot << ", ";
    }else{
        gnuplot << "plot ";
    }

    for(int i = 0; i < curves.size(); i ++)    {
        gnuplot << gnuplot.file1d(curves[i]) << " u 1:2 with lines notitle";
        if(i<curves.size() -1)  {
            gnuplot << ", ";
        }
    }
}

CVisual::gnuplotWrapper::gnuplotWrapper() {
    pngparams();
    xlabel = "xlabel";
    ylabel = "ylabel";
    minval=0;
    maxval=-1;
    yminval = 0;
    ymaxval = -1;
    m_gnuplot = make_shared<gnuplotio::Gnuplot>();
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_xlabel(string xlabel_) {
    xlabel = xlabel_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_ylabel(string ylabel_) {
    ylabel = ylabel_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_filename(string filename_) {
    outputfilename = filename_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_terminaltype(string terminaltype_) {
    terminal_type = terminaltype_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_height(int height_) {
    height = height_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_width(int width_) {
    width = width_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_minmax(double minval_, double maxval_) {
    minval = minval_;
    maxval = maxval_;
    return *this;
}

CVisual::gnuplotWrapper &CVisual::gnuplotWrapper::set_yminmax(double yminval_, double ymaxval_) {
    yminval = yminval_;
    ymaxval = ymaxval_;
    return *this;
}

void CVisual::gnuplotWrapper::setLabel(gnuplotio::Gnuplot &gnuplot) {
    gnuplot << "set xlabel \"" << xlabel << "\""<< endl;
    gnuplot << "set ylabel \"" << ylabel << "\""<< endl;
}

void CVisual::gnuplotWrapper::setTerminal(gnuplotio::Gnuplot &gnuplot) {
    gnuplot << "set term " << terminal_type << " size " << width << ",  " << height << endl;

}

void CVisual::gnuplotWrapper::setMinMaxRange(gnuplotio::Gnuplot &gnuplot) {
    if(minval<maxval)   gnuplot << "set xrange [" << minval << ":" << maxval << "]" << endl;
}

void CVisual::gnuplotWrapper::setYMinMaxRange(gnuplotio::Gnuplot &gnuplot) {
    if(yminval<ymaxval)   gnuplot << "set yrange [" << yminval << ":" << ymaxval << "]" << endl;
}

void CVisual::gnuplotWrapper::setOutput(gnuplotio::Gnuplot &gnuplot) {
    if(terminal_type!="dumb")    gnuplot << "set output '" << outputfilename << "'" << endl;

}

void CVisual::gnuplotWrapper::config(gnuplotio::Gnuplot &gnuplot) {
    setTerminal(gnuplot);
    setOutput(gnuplot);
    setLabel(gnuplot);
    setMinMaxRange(gnuplot);
}

void CVisual::gnuplotWrapper::config() {
    setTerminal(*m_gnuplot);
    setOutput(*m_gnuplot);
    setLabel(*m_gnuplot);
    setMinMaxRange(*m_gnuplot);
}

void CVisual::gnuplotWrapper::pngparams() {
    terminal_type = "pngcairo enhanced ";
    outputfilename = "test.png";
    width = 1200;
    height = 900;
}

void CVisual::gnuplotWrapper::plotLinearRegression(vector<double> &coefs, double rSquared,
                                                   const vector<vector<double>> &rx_logry) const {
    ostringstream osstitle;
    osstitle << setprecision(4) << "log(p)=" << coefs.at(0) << "x+" << coefs.at(1) << " R^2=" << rSquared;
    *m_gnuplot << "f(x)=" << coefs.at(0) << "*x+" << coefs.at(1) << endl;
    *m_gnuplot << "plot " << m_gnuplot->file1d(rx_logry)
              << "  u 1:2 w points lc rgb'green' title 'p-value in log scale' ,"
                 " f(x) w line lc rgb'red' linetype '-' title '" << osstitle.str() << "'" << endl;
}
void CVisual::SImageSize::calcPixel() {
    m_wpx=floor(m_width*m_dpi);
    m_hpx=floor(m_height*m_dpi);
}

CVisual::SImageSize::SImageSize() {
    init();
}

CVisual::SImageSize::SImageSize(int width, int height, int dpi) {
    m_dpi = dpi;
    m_width = width;
    m_height = height;
    calcPixel();
}

void CVisual::SImageSize::init() {
    m_dpi = 72;
    m_width = 7;
    m_height = 5;
    calcPixel();
}