#ifndef SPECTRASTDENOISER_HPP_
#define SPECTRASTDENOISER_HPP_

#include "SpectraSTPeakList.hpp"
#include <vector>

using namespace std;

class SpectraSTDenoiserFeature {

public:

    SpectraSTDenoiserFeature();

    void initialize(unsigned int nBins);

    void calculateLogOdds(double logOddsPrior, unsigned int minCount = 2);

    void smoothLogOdds();

    void printLogOdds(ofstream &fout);

    void incrementSignal(unsigned int bin);

    void incrementNoise(unsigned int bin);

    unsigned int getNumBins() { return (numBins); }

    double getLogOdds(unsigned int bin) { if (bin < numBins) return (logOdds[bin]); else return (0.0); }

    void setLogOdds(unsigned int bin, double value) { if (bin < numBins) logOdds[bin] = value; }

private:

    unsigned int numBins;
    vector<unsigned int> signalCount;
    vector<unsigned int> noiseCount;
    vector<double> logOdds;


};


class SpectraSTDenoiser {

public:

    SpectraSTDenoiser();

    virtual ~SpectraSTDenoiser();

    void initializeTraining();

    void addTrainingSpectrum(SpectraSTPeakList *pl);

    void generateBayesianModel();

    double predictPrior(vector<double> &priorData);

    void filter(SpectraSTPeakList *peakList, unsigned int maxNumPeaks, double minSignalProb = -1.0);

    void useDefault();

    void readModel(ifstream &fin);

    void printModel(ofstream &fout);

    bool isFilterReady() { return (m_ready); }

private:

    bool m_ready;

    unsigned int m_numSpectra;
    unsigned int m_numSignal;
    unsigned int m_numNoise;

    double m_defaultPrior;

    vector<double> m_priors;

    vector<double> m_sisters;

    SpectraSTDenoiserFeature *m_featureMzPos;
    SpectraSTDenoiserFeature *m_featureRank;
    // SpectraSTDenoiserFeature* m_featureDeltaInt;
    SpectraSTDenoiserFeature *m_featureWithSister;
    SpectraSTDenoiserFeature *m_featureWithoutSister;
    SpectraSTDenoiserFeature *m_featureWithComplement;
    SpectraSTDenoiserFeature *m_featureWithoutComplement;

    vector<double> m_priorCoefficients;

    vector<vector<double> > m_priorRegressionData;

    vector<int> m_numPeaksAtRank;
    vector<double> m_sumIntensityAtRank;
    vector<double> m_sumSqIntensityAtRank;

    static bool sortPeaksByOddsDesc(pair<double, Peak> a, pair<double, Peak> b);

};

#endif
