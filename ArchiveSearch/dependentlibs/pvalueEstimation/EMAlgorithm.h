//
// Created by wulong on 1/30/19.
//

#ifndef MYTOOL_EMALGORITHM_H
#define MYTOOL_EMALGORITHM_H

#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <iomanip>
#include <map>
#include <algorithm>
#include "../../../librarymsms/Util.h"
#include "Visual.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions.hpp>
using namespace std;


double calcSum(const vector<double> &x, bool stable = true);

void visualizeHistogram(const vector<double> &s);

class Distribution
{
public:
    virtual ~Distribution(){}
    virtual double getCdf(double x)=0;
    virtual double getPdf(double x)=0;
    virtual void update(const vector<double> &sample, const vector<double> &R)=0;
    virtual void display()=0;
    virtual void visualize()=0;
    virtual string getname()=0;
};

class Gaussian: public Distribution
{
    double m_mu;
    double m_sigma;
public:
    Gaussian(double mu, double sigma);
    double getCdf(double x) override;
    double getPdf(double x) override;
    void update(const vector<double> &sample, const vector<double> &R) override;
    void display() override;
    string getname() override{
        ostringstream oss;
        oss << "Gaussian_" << m_mu << "_" << m_sigma << "_";
        return oss.str();
    }
    double generate()
    {
        random_device rd;
        mt19937 gen(rd());
        std::normal_distribution<> distribution{m_mu, m_sigma};
        return distribution(gen);
    }
    void visualize()
    {
        cout << "---------Gaussian-----------" << endl;
        cout << "mu\tsigma" << endl;
        cout << m_mu << "\t" << m_sigma << endl;
        int N = 10000;

        vector<double> s(N,0);
        for(int i =  0; i < N; i ++)
        {
            s[i]=generate();
        }

        double mean = statistic::calcmean(s);
        double sigma = statistic::calcstd(s);
        cout << "mean=" << mean << "\tsigma=" << sigma << endl;
        visualizeHistogram(s);
        CVisual::gnuplot_histogram(s,100);
        cout << "----------Gaussian Done--------------------" << endl;
    }
};

class Gumbel: public Distribution
{
    double m_mu;
    double m_sigma;
    // We could also present the model using alpha and beta
    double m_alpha;
    double m_beta;
public:

    Gumbel(double mu, double sigma);
    Gumbel(double alpha, double beta, bool natural):m_alpha(alpha),m_beta(beta){
        const double gamma = 0.57721566490153286060651209008240243;
        const double pi = 3.141592653589793238462643383279;
        m_mu = alpha+beta*gamma;
        m_sigma = pi/sqrt(6)*beta;
    }
    void updateAlphaBeta();
    double getCdf(double x) {
        double a = (m_alpha-x)/m_beta;
        return exp(-exp(a));
    }
    double getPdf(double x) override;
    void update(const vector<double> &sample, const vector<double> &R) override;
    void display() override;
    string getname() override{
        ostringstream oss;
        oss << "Gumbel_" << m_mu << "_" << m_sigma << "_";
        return oss.str();
    }
    double generate()
    {
        random_device rd;
        mt19937 gen(rd());
        std::extreme_value_distribution<> distribution{m_alpha, m_beta};
        return distribution(gen);

    }
    void visualize()
    {
        cout << "----------Gumbel--------------------" << endl;
        cout << "mu\tsigma\talpha\tbeta" << endl;
        cout << m_mu << "\t" << m_sigma << "\t" << m_alpha << "\t" << m_beta << endl;
        int N = 10000;
        vector<double> s(N,0);
        for(int i =  0; i < N; i ++)
        {
            s[i]=generate();
        }

        double mean = statistic::calcmean<double>(s);
        double sigma = statistic::calcstd<double>(s);
        cout << "mean=" << mean << "\tsigma=" << sigma << endl;


        visualizeHistogram(s);
        CVisual::gnuplot_histogram(s,100);

        cout << "----------Gumbel-----Done---------------" << endl;


    }



};


class Gamma: public Distribution
{
    double m_shape;
    double m_scale;
    boost::math::gamma_distribution<> m_gamma;
    double m_mu;
    double m_std;
public:
    Gamma(double shape, double scale):m_shape(shape),m_scale(scale),m_gamma(shape, scale){
        m_mu = boost::math::mean(m_gamma);
        m_std = boost::math::standard_deviation(m_gamma);
    }
    double getCdf(double x) override{return boost::math::cdf(m_gamma,x);}
    double getPdf(double x){return boost::math::pdf(m_gamma, x) ;}
    void display()override
    {
        cout << "Gamma(" << m_mu << " " << m_std << ") ";
    }
    string getname() override{
        ostringstream oss;
        oss << "_Gamma_" << m_mu << "_" << m_std << "_";
        return oss.str();
    }

    void visualize()override
    {
        cout << "----------Gamma--------------------" << endl;
        cout << "mu\tsigma\tshape\tscale" << endl;
        cout << m_mu << "\t" << m_std << "\t" << m_shape << "\t" << m_scale << endl;
        int N = 10000;
        vector<double> s(N,0);
        for(int i =  0; i < N; i ++)
        {
            s[i]=generate();
        }

        double mean = statistic::calcmean<double>(s);
        double sigma = statistic::calcstd<double>(s);
        cout << "mean=" << mean << "\tsigma=" << sigma << endl;


        visualizeHistogram(s);
        CVisual::gnuplot_histogram(s,100);

        cout << "----------Gamma-----Done---------------" << endl;
    }
    void update(const vector<double> &sample, const vector<double> &R)override{
        double sum = 0;
        double sumr = 0;
        // calculate mean
        for(int j = 0; j < sample.size(); j ++)
        {
            sum += (sample[j]*R[j]);
            sumr += R[j];
        }
        //double tmu = m_mu;
        m_mu = sum/sumr;

        sum = 0;
        for(int j = 0; j < sample.size(); j ++)
        {
            sum += (sample[j]-m_mu)*(sample[j]-m_mu)*R[j];
        }
        sum /= sumr;
        m_std = sqrt(sum);

        updateShapeScale();
    }

    void updateShapeScale()
    {
        m_scale = m_std*m_std/m_mu;
        m_shape = m_mu*m_mu/m_std/m_std;
        m_gamma = boost::math::gamma_distribution<>(m_shape, m_scale);
    }
    double generate()
    {
        random_device rd;
        mt19937 gen(rd());
        std::gamma_distribution<> distribution{m_shape, m_scale};
        return distribution(gen);

    }
};



class EM
{
public:
    vector<vector<double>> r;
    vector<double> &Sample;
    vector<Distribution*> D;
    vector<double> W;
public:
    EM(vector<Distribution*> d, vector<double> &sample);

    double pdfMix(double x);
    void display();
    void visualizeSample(){
        cout << "----------------------sample----distribution-----------" << endl;
        CVisual::gnuplot_histogram(Sample,100);
        cout << "----------------------sample----distriubtion-----------" << endl;
    }
    void updateR();
    void updateModel();
    double loglikelihood();
    string getNames()
    {
        string ret = "";
        for(int i = 0; i < D.size(); i ++)
        {
            ret += D[i]->getname();
        }
        return ret;
    }

    void run(int n = 100);
    void getFdrPep(vector<double> &s, vector<double> &scoreFDR, vector<double> &scorePEP, vector<double> &scorePos, vector<double> &scoreNeg)
    {
        double maxdp = *max_element(Sample.begin(), Sample.end());
        double mindp = *min_element(Sample.begin(), Sample.end());

        int samplenum = 1000;
        s.assign(samplenum,0);
        for(int i = 0; i < samplenum ; i ++)
        {
            s[i] = mindp + i*(maxdp-mindp)/(samplenum-1);
        }

        vector<double> positiveArea(samplenum,0), negativeArea(samplenum,0);
        scoreFDR.assign(samplenum, 0);
        scorePEP.assign(samplenum,0);
        int pD = 1, nD= 0;
     //   cout << "score\tFDR\tPEP" << endl;
        for(int i = 0; i < samplenum; i ++)
        {
            positiveArea[i] = (1- D[pD]->getCdf(s[i])) * W[pD];
            negativeArea[i] = (1- D[nD]->getCdf(s[i])) * W[nD];

            scoreFDR[i] = negativeArea[i]/(positiveArea[i]+negativeArea[i]);
            double positiveValue = D[pD]->getPdf(s[i]) * W[pD];
            double negativeValue = D[nD]->getPdf(s[i]) * W[nD];
            scorePos.push_back(positiveValue);
            scoreNeg.push_back(negativeValue);
            scorePEP[i] = negativeValue/ (positiveValue+negativeValue);
           // cout << s[i] << "\t" << scoreFDR[i] << "\t" << scorePEP[i] << endl;
        }
    }
};

#endif //MYTOOL_EMALGORITHM_H
