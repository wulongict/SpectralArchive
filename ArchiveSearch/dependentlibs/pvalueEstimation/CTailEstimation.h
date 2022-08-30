//
// Created by wulong on 4/13/19.
//

#ifndef MYTOOL_CTAILESTIMATION_H
#define MYTOOL_CTAILESTIMATION_H

#include <vector>
#include <string>
#include <memory>
using namespace std;
class ICHistdata;
namespace gnuplotio
{
    class Gnuplot;
}


bool linear_regression_on_logCDF(vector<vector<double>> &rx_logry, vector<double> &coefs, double &r2);


struct SLinearRegressionModel{
    vector<double> coefs;
    int sampleSize;
    double rSquare;
    double minY;
    double maxY;
    SLinearRegressionModel(){
        coefs.assign(2,0);
        sampleSize=0;
        rSquare=0;
        minY=0;
        maxY=0;
    }

    string toJsonString();
};

class CTailEstimation {
    shared_ptr<ICHistdata> m_histogram;
    vector<double> m_x;
    vector<double> m_y;
    bool succeed;
    vector<double> px;
    vector<double> logy;
    vector<vector<double>> px_logy;
    vector<double> coefs;
    double r2;
    bool m_verbose;
    double m_total;
    int m_sampleSize;
    double m_min_y;
    double m_max_y;
public:
    CTailEstimation(shared_ptr<ICHistdata> ptrHistogram, bool verbose);
    double getTopHitpvalue(double x, float fraction, bool verbose);
    double getlogpvalue(double x);
    bool buildLinearRegressionModel(bool noZero, const string& outputbasename, bool exportToFile, bool plot = false, SLinearRegressionModel *lrPtr=nullptr);

private:
    void plotModel(vector<double> coefs, double r2, vector<vector<double>> &rx_logry, string outputname, string termtype);
    void getLog10SurvivalFn(double tailFraction=0.1);
    void getSurvivalFn(bool noZero, string outputbasename, bool outputToFile);
    void terminalVisualization(vector<double> coefs, double rSquared, vector<vector<double>> &rx_logry);
    void pngVisualization(const string& outputname, double r2, vector<double> coefs, vector<vector<double>> &rx_logry);
};


#endif //MYTOOL_CTAILESTIMATION_H
