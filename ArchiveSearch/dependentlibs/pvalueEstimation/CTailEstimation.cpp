//
// Created by wulong on 4/13/19.
//

#include <cmath>
#include <ostream>
#include <iostream>
#include "CTailEstimation.h"
#include "CLinearRegression.h"
#include "Visual.h"
#include "../../../librarymsms/Util.h"
#include "../../../External/gnuplot-iostream/gnuplot-iostream.h"
#include "../../../librarymsms/ICHistdata.h"

using namespace std;

void prepare_linear_regression(const vector<vector<double>> &rx_logry, double minCDF, double maxCDF, vector<double> &b,
                               vector<vector<double>> &A) {
    cout << "[Warning] use the following function instead! This one will be deleted later" << endl;
    int count = 0;
    for (const auto & i : rx_logry) {
        if (i[1] < minCDF or i[1] > maxCDF) {
            continue;
        } else {
            A[count][0] = i[0];
            b[count] = i[1];
            count++;
        }
    }
    b.resize(count);
    A.resize(count);
}

void prepare_linear_regression(vector<double> &rx, vector<double> &logry, double minCDF, double maxCDF,
                               vector<double> &b, vector<vector<double>> &A) {
    int count = 0;
    for (int i = 0; i < logry.size(); i++) {
        if (logry[i] < minCDF or logry[i] > maxCDF) {
            continue;
        } else {
            A[count][0] = rx[i];
            b[count] = logry[i];
            count++;
        }
    }
    b.resize(count);
    A.resize(count);
}

//
//       S_Y S_{XX}-S_X S_{XY}
// b   = ---------------------------------
//       n  S_{XX} -  S_X S_X
//
//        n S_{XY} - S_Y S_X
// k   =   ---------------------------
//        n  S_{XX} -  S_X S_X

// 1D linear regression problem.
// y = kx +b
bool simple_1D_linear_regression(std::vector<double> x, std::vector<double> y, double &k, double &b, double &r2){
    cout << "running " << __FUNCTION__  << endl;
    bool ret = false;
    int n = x.size();
    double sx = 0, sy=0, sxx=0, sxy=0;
    for(int i = 0; i < x.size(); i ++){
        sx += x[i];
        sy += y[i];
        sxy += x[i]*y[i];
        sxx += x[i] * x[i];
    }

    double denominator = n * sxx - sx * sx;
    double EPSILON = 1e-8;
    if(fabs(denominator)<EPSILON){
        ret = false;
    }else{
        b = (sy * sxx - sx * sxy) / denominator;
        k = (n * sxy - sy * sx)/ denominator;
        // calculate r2
        double y_bar = sy / n;
        double SSR = 0, SST = 0;
        for(int i = 0; i < x.size(); i ++){
            double t = y[i]-(k * x[i] + b), tt = y[i]-y_bar;
            SSR += t*t;
            SST += tt * tt;
        }
//        std::cout << SSR << "\t" << SST << std::endl;
        r2 = 1- SSR/SST; // is it possible that SST = 0?
        ret = true;
    }

    return ret;

}


// the output of linear regression
struct CLinRegOut{
    double k;
    double b;
    double r2;
};

// diff of the two function
// linear system
// Ax=b
bool linear_regression_on_logCDF(vector<double> &x, vector<double> &y, vector<double> &coefs, double &r2,
                                 bool verbose, double &min_y, double &max_y, int &sampleSize) {
    bool ret = false;
    double MIN_R_SQUARE_ALLOWED = 0.95;
    r2 = 0.001;

    while (min_y < max_y ) {
        vector<double> b(x.size(), 0);
        vector<vector<double>> A(x.size(), vector<double>{0});
        prepare_linear_regression(x, y, min_y, max_y, b, A);
        sampleSize = b.size();
        CLinRegOut lro;
        {

            vector<double> x;
            for(int i = 0; i < A.size(); i ++){
                x.push_back(A[i][0]);
            }
            simple_1D_linear_regression(x,b,lro.k,lro.b,lro.r2);
        }

        if (sampleSize < 50 ) {
//            cout << "Warining: Too few sample: " << sampleSize << endl;
            break;
        }
        CLinearRegression lr(A, b);
        lr.solve(verbose);
        coefs = lr.getCoefficients();
        r2 = lr.getRsquare();


        if(fabs(coefs[0]-lro.k) + fabs(coefs[1]-lro.b) > 1e-3 or fabs(r2-lro.r2)>1e-3  ){
            cout << "[Error] inconsistence lr resutls " << endl;
            cout << setprecision(5) << "#sample:" << sampleSize << " on y axis interval [" << min_y << ", " << max_y << "] Model: logCDF = " << coefs.at(0) << "x+" << coefs.at(1) << "\tR^2=" << r2
                 << endl;
            cout << "Mine: " << lro.k << "x + " << lro.b << " R2 = " << lro.r2 << endl;
        }else{
            cout << fabs(coefs[0]-lro.k) + fabs(coefs[1]-lro.b) << "\t"  << fabs(r2-lro.r2) << endl;

            cout << "[Good] good lr resutls " << endl;
            cout << setprecision(5) << "#sample:" << sampleSize << " on y axis interval [" << min_y << ", " << max_y << "] Model: logCDF = " << coefs.at(0) << "x+" << coefs.at(1) << "\tR^2=" << r2
                 << endl;
            cout << "Mine: " << lro.k << "x + " << lro.b << " R2 = " << lro.r2 << endl;

        }

        if(r2>=MIN_R_SQUARE_ALLOWED){
            ret=true;
            break;
        }
        if (verbose){
            cout << setprecision(5) << "#sample:" << sampleSize << " on y axis interval [" << min_y << ", " << max_y << "] Model: logCDF = " << coefs.at(0) << "x+" << coefs.at(1) << "\tR^2=" << r2
                 << endl;

        }
        min_y += 0.1;
    }
    return ret;
}

bool linear_regression_on_logCDF(vector<vector<double>> &rx_logry, vector<double> &coefs, double &r2) {
    SimpleTimer st("Linear Regression");
    double MIN_R_SQUARE_ALLOWED = 0.99;
    r2 = 0.001;
    double minCDF = -2.5, maxCDF = -1;
    bool success = true;
    while (r2 < MIN_R_SQUARE_ALLOWED and minCDF < maxCDF - 1) {
        cout << "log CDF range: " << minCDF << " --> " << maxCDF << endl;
        vector<double> b(rx_logry.size(), 0);
        vector<vector<double>> A(rx_logry.size(), vector<double>{0});
        prepare_linear_regression(rx_logry, minCDF, maxCDF, b, A);
        cout << "Info: sample number = " << b.size() << endl;
        if (b.size() < 10) {
            success = false;
            break;
        }
        CLinearRegression lr(A, b);
        lr.solve(true);
        coefs = lr.getCoefficients();
        r2 = lr.getRsquare();
        cout << setprecision(5) << "Model: logCDF = " << coefs[0] << "x+" << coefs[1] << "\tR^2=" << r2 << endl;
        minCDF += 0.1;
    }
    return success;
}


void CTailEstimation::terminalVisualization(vector<double> coefs, double rSquared, vector<vector<double>> &rx_logry) {
    if (coefs.size() < 2) {
        return;
    }
    CVisual::gnuplotWrapper info;
    info.set_xlabel("dot product")
            .set_ylabel("log10(1-cdf) or log_{10}(p-value)")
            .set_height(40)
            .set_minmax(0, 1)
            .set_yminmax(-20, 0)
            .set_width(100)
            .set_terminaltype("dumb");
    info.config();
    info.plotLinearRegression(coefs, rSquared, rx_logry);
}

void CTailEstimation::pngVisualization(const string& outputname, double r2, vector<double> coefs,
                                       vector<vector<double>> &rx_logry) {
    if (coefs.size() < 2) {
        return;
    }
    CVisual::gnuplotWrapper info;
    CVisual::SImageSize imsize;
    info.set_xlabel("dot product")
            .set_ylabel("log_{10}(1-CDF(x))")
            .set_minmax(0, 1)
            .set_yminmax(-20, 0)
            .set_width(imsize.m_wpx)
            .set_height(imsize.m_hpx)
            .set_filename(outputname + "_rx_logcdf_allscores_adaptive_lr_model.png");
    info.config();
    info.plotLinearRegression(coefs, r2, rx_logry);
}

void CTailEstimation::getSurvivalFn(bool noZero, string outputbasename, bool outputToFile) {
    if (m_histogram->empty()) {
        cout << "Error: empty score vector " << m_histogram->size() << endl;
        throw logic_error("Score vector for tail estimation is empty!");
    }

    m_x = m_histogram->getX();
    m_y = m_histogram->getSurvivalFn(noZero);
    m_total = m_histogram->gettotal();

    if (outputToFile) {
        CVisual::gnuplotWrapper info;
        CVisual::SImageSize imsize;
        string pngfilename = outputbasename + to_string("", "_", "nozero", noZero) + "_cdf_scores.png";
        info.set_xlabel("normalized dot product")
                .set_ylabel("1-CDF(x)")
                .set_height(imsize.m_hpx)
                .set_width(imsize.m_wpx)
                .set_filename(pngfilename);
        CVisual::gnuplot_curve_topng(m_x, m_y, "lines", info);
        exportTable(m_y, outputbasename + "_cdf_allscores.txt", '\t');
    }
}

void CTailEstimation::getLog10SurvivalFn(double tailFraction) {
    for (int i = m_histogram->size() - 1; i >= 0; i--) {
        if (m_y[i] > 0) {
            if (px.empty() or m_y[i] - m_y[i + 1] > EPSILON) {
                px.push_back(m_x[i]);
                logy.push_back(log10(m_y[i]));
                px_logy.push_back(vector<double>{m_x[i], log10(m_y[i])});
            }
        }
        if (m_y[i] > tailFraction) {
            break;
        }
    }
}

bool CTailEstimation::buildLinearRegressionModel(bool noZero, const string& outputbasename, bool exportToFile, bool plot, SLinearRegressionModel *lrPtr) {
    getSurvivalFn(noZero, outputbasename, exportToFile);
    getLog10SurvivalFn();

    succeed = linear_regression_on_logCDF(px, logy, coefs, r2, m_verbose, m_min_y, m_max_y, m_sampleSize);
    if (succeed and plot) {
        string pngfilename = outputbasename + to_string("", "_", "nozero", noZero) + "_logpValue_scores.png";
        CVisual::gnuplotWrapper info;
        CVisual::SImageSize imsize;
        info.set_xlabel("normalized dot product")
                .set_ylabel("log10(1-CDF(x))")
                .set_height(imsize.m_hpx).set_width(imsize.m_wpx).set_minmax(0,1)
                .set_filename(pngfilename);

        CVisual::gnuplot_curve_topng(px, logy, "points pt 7 ", info);
        cout << "terminal visualization and png visualization" << endl;
        plotModel(coefs, r2, px_logy, outputbasename, "dumb");
        plotModel(coefs, r2, px_logy, outputbasename, "png");
        exportTable(px_logy, outputbasename + "_rx_logcdf_allscores.txt", '\t');
    }
    if(lrPtr!=nullptr and succeed)
    {
        lrPtr->sampleSize = m_sampleSize;
        lrPtr->coefs = coefs;
        lrPtr->rSquare = r2;
        lrPtr->maxY = m_max_y;
        lrPtr->minY = m_min_y;
    }
    return succeed;
}

double CTailEstimation::getlogpvalue(double x) {
    if (coefs.empty() or not succeed) {
        if (m_verbose)cout << "Linear regression Model empty:" << endl;
        return 0;
    }
    double ret = coefs.at(0) * x + coefs.at(1);
    //cout << "log(pvlaue) = " << coefs.at(0) << "* x + " << coefs.at(1) <<   "\tscore:  " << x << " log(pvalue): " << ret << endl;
    return ret;
}

double CTailEstimation::getTopHitpvalue(double x, float fraction, bool verbose) {
    double constant = 1;
    // we are having N/fraction
    double log10_pvalue = getlogpvalue(x);

    double pvalue = exp10(log10_pvalue);
    if (std::isnan(log10_pvalue) or std::isnan(pvalue)) {
//        cout << "log10 pvalue and pvalue contain nan: "<<  log10_pvalue << " " << pvalue << endl;
    }
    double N = m_total / fraction / constant;
    // We look into N candidates from archive, and find the top hit score: s
    // we know by random chance, the probability of getting a candidates of score s is
    // Pr ( score >= s) = pvalue
    // to random observe such a score, on N random candiates, is
    // Pr(s_max >= s ) = 1-Pr( s_1 < s)Pr(s_2 < s)...Pr(s_N < s) = 1- (1-pvalue)^N
    // that's why we have the following correction for pvalue of top hit!
    double adj_pvalue = 1 - pow(1 - pvalue, N);  // How to choose N?

    if (std::isnan(adj_pvalue) or pvalue>0.1) {
        // Aug 19.. This piece of code might change a lot.
//        if(verbose) cout << "adj_pvalue nan: p: " << pvalue << " n: " << N << " 1-(1-p)^n = " << adj_pvalue << endl;
        adj_pvalue = 1;
    }

    if (fabs(pvalue - 1.0) < 1e-16L) {
        adj_pvalue = 1;
    }

    if (not succeed) {
        adj_pvalue = 2;
    }
    return adj_pvalue;
}

CTailEstimation::CTailEstimation(shared_ptr<ICHistdata> ptrHistogram, bool verbose) {
    m_histogram = ptrHistogram;
    m_verbose = verbose;
    succeed = false;
    m_min_y=-2.5, m_max_y=-1;
    m_sampleSize = 0;
}

void CTailEstimation::plotModel(vector<double> coefs, double r2, vector<vector<double>> &rx_logry, string outputname,
                                string termtype) {
    if (termtype == "dumb") {
        terminalVisualization(coefs, r2, rx_logry);
    } else if (termtype == "png") {
        pngVisualization(std::move(outputname), r2, coefs, rx_logry);
    } else {
        cout << "[Error] Invalid terminal type " << termtype << "Please try dumb or png" << endl;
    }
}

string SLinearRegressionModel::toJsonString() {
    ostringstream  oss;
    oss << R"("lrModel": { "coef_a": )" << this->coefs.at(0)
        << ", \"coef_intercept\": " << this->coefs.at(1) << ", \"rSquare\": " << this->rSquare
        << ", \"sampleSize\": " << this->sampleSize
        << ", \"maxLog10Pvalue\": " << this->maxY
        << ", \"minLog10Pvalue\": " << this->minY
        << " }"
        << endl;
    return oss.str();
}
