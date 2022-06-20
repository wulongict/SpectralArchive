#include "SpectraSTDenoiser.hpp"
#include "FileUtils.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>

#ifndef __LGPL__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>

#endif

using namespace std;

SpectraSTDenoiserFeature::SpectraSTDenoiserFeature() {}

void SpectraSTDenoiserFeature::initialize(unsigned int nBins) {

    numBins = nBins;
    signalCount.insert(signalCount.end(), nBins, 1);
    noiseCount.insert(noiseCount.end(), nBins, 1);
    logOdds.insert(logOdds.end(), nBins, 0.0);

}

void SpectraSTDenoiserFeature::incrementSignal(unsigned int bin) {
    if (bin < numBins) signalCount[bin]++;
}

void SpectraSTDenoiserFeature::incrementNoise(unsigned int bin) {
    if (bin < numBins) noiseCount[bin]++;
}

void SpectraSTDenoiserFeature::calculateLogOdds(double logOddsPrior, unsigned int minCount) {

    unsigned int cumSignalCount = 0;
    unsigned int cumNoiseCount = 0;
    int lastGoodIndex = numBins;

    for (int index = numBins - 1; index >= 0; index--) {

        cumSignalCount += signalCount[index];
        cumNoiseCount += noiseCount[index];

        if (cumSignalCount + cumNoiseCount >= minCount) {
            // enough
            for (int i = lastGoodIndex - 1; i >= index; i--) {
                if (cumSignalCount < 1 || cumNoiseCount < 1) {
                    logOdds[i] = 0.0; // should never happen, as the signal and noise vectors are initialized with one's
                } else {
                    logOdds[i] = log((double) cumSignalCount) - log((double) cumNoiseCount) - logOddsPrior;
                }
            }

            cumSignalCount = 0;
            cumNoiseCount = 0;
            lastGoodIndex = index;

        }
    }

    if (lastGoodIndex != 0) {
        for (int i = lastGoodIndex - 1; i >= 0; i--) {
            logOdds[i] = logOdds[lastGoodIndex];
        }
    }

}

void SpectraSTDenoiserFeature::smoothLogOdds() {

    if (numBins < 4) return; // avoid out-of-bound

    vector<double> smooth(logOdds); // copy

    logOdds[0] = (7 * smooth[0] + 6 * smooth[1] + 3 * smooth[2] - 2 * smooth[3]) / 14.0;
    logOdds[1] = (7 * smooth[1] + 6 * (smooth[2] + smooth[0]) + 3 * smooth[3] - 2 * smooth[4]) / 20.0;
    logOdds[2] = (7 * smooth[2] + 6 * (smooth[3] + smooth[1]) + 3 * (smooth[4] + smooth[0]) - 2 * smooth[5]) / 23.0;
    for (int index = 3; index < numBins - 3; index++) {
        logOdds[index] = (7 * smooth[index] + 6 * (smooth[index + 1] + smooth[index - 1]) +
                          3 * (smooth[index + 2] + smooth[index - 2]) - 2 * (smooth[index + 3] + smooth[index - 3])) /
                         21.0;
    }
    logOdds[numBins - 3] = (7 * smooth[numBins - 3] + 6 * (smooth[numBins - 2] + smooth[numBins - 4]) +
                            3 * (smooth[numBins - 1] + smooth[numBins - 5]) - 2 * smooth[numBins - 6]) / 23.0;
    logOdds[numBins - 2] =
            (7 * smooth[numBins - 2] + 6 * (smooth[numBins - 1] + smooth[numBins - 3]) + 3 * smooth[numBins - 4] -
             2 * smooth[numBins - 5]) / 20.0;
    logOdds[numBins - 1] =
            (7 * smooth[numBins - 1] + 6 * smooth[numBins - 2] + 3 * smooth[numBins - 3] - 2 * smooth[numBins - 4]) /
            14.0;


}

void SpectraSTDenoiserFeature::printLogOdds(ofstream &fout) {

    fout << numBins;
    for (int index = 0; index < numBins; index++) {
        fout << '\t' << logOdds[index];
    }

}


SpectraSTDenoiser::SpectraSTDenoiser() :
        m_numSpectra(0),
        m_numSignal(0),
        m_numNoise(0),
        m_ready(false),
        m_numPeaksAtRank(200, 0),
        m_sumIntensityAtRank(200, 0.0),
        m_sumSqIntensityAtRank(200, 0.0) {

    m_featureMzPos = new SpectraSTDenoiserFeature();
    m_featureRank = new SpectraSTDenoiserFeature();
//  m_featureDeltaInt = new SpectraSTDenoiserFeature();
    m_featureWithSister = new SpectraSTDenoiserFeature();
    m_featureWithoutSister = new SpectraSTDenoiserFeature();
    m_featureWithComplement = new SpectraSTDenoiserFeature();
    m_featureWithoutComplement = new SpectraSTDenoiserFeature();

}

SpectraSTDenoiser::~SpectraSTDenoiser() {

    if (m_featureRank) delete (m_featureRank);
//  if (m_featureDeltaInt) delete (m_featureDeltaInt);
    if (m_featureMzPos) delete (m_featureMzPos);
    if (m_featureWithSister) delete (m_featureWithSister);
    if (m_featureWithoutSister) delete (m_featureWithoutSister);
    if (m_featureWithComplement) delete (m_featureWithComplement);
    if (m_featureWithoutComplement) delete (m_featureWithoutComplement);

}

void SpectraSTDenoiser::initializeTraining() {

    m_sisters.push_back(-18.01057);
    m_sisters.push_back(-17.0027);
    m_sisters.push_back(-9.005285);
    m_sisters.push_back(-2.016);
    m_sisters.push_back(-1.008);
    m_sisters.push_back(1.008);
    m_sisters.push_back(2.016);
    m_sisters.push_back(9.005285);
    m_sisters.push_back(17.0027);
    m_sisters.push_back(18.01057);
    m_sisters.push_back(27.9949);
    m_sisters.push_back(57.02146);
    m_sisters.push_back(71.03711);
    m_sisters.push_back(79.966331);
    m_sisters.push_back(87.03202);
    m_sisters.push_back(97.05276);
    m_sisters.push_back(97.976896);
    m_sisters.push_back(99.06841);
    m_sisters.push_back(101.0476);
    m_sisters.push_back(103.0091);
    m_sisters.push_back(113.084);
    m_sisters.push_back(114.0429);
    m_sisters.push_back(115.0269);
    m_sisters.push_back(115.987461);
    m_sisters.push_back(128.0585);
    m_sisters.push_back(128.0949);
    m_sisters.push_back(129.0425);
    m_sisters.push_back(131.0404);
    m_sisters.push_back(137.0589);
    m_sisters.push_back(147.0684);
    m_sisters.push_back(156.1011);
    m_sisters.push_back(163.0633);
    m_sisters.push_back(186.0793);

    m_featureRank->initialize(150);
    m_featureMzPos->initialize(100);
//  m_featureDeltaInt->initialize(20);
    m_featureWithSister->initialize(m_sisters.size());
    m_featureWithoutSister->initialize(m_sisters.size());
    m_featureWithComplement->initialize(7);
    m_featureWithoutComplement->initialize(7);

    m_priorCoefficients.resize(17, 0.0);

}

void SpectraSTDenoiser::useDefault() {

    m_defaultPrior = 0.425597;
    m_featureRank->initialize(150);

    float rankLogOdds[] = {4.1824, 4.02584, 3.84031, 3.5457, 3.35526, 3.18706, 3.03382, 2.89556, 2.77159, 2.66182,
                           2.55555, 2.45662, 2.35724, 2.26668, 2.17767, 2.09578, 2.01949, 1.94732, 1.87646, 1.80855,
                           1.74037, 1.67251, 1.60289, 1.54318, 1.48164, 1.42276, 1.36782, 1.31481, 1.26229, 1.2139,
                           1.16323, 1.10972, 1.05745, 1.00635, 0.959576, 0.913063, 0.871409, 0.828742, 0.786954,
                           0.745091,
                           0.706025, 0.661849, 0.619583, 0.582152, 0.546317, 0.510238, 0.473625, 0.432658, 0.394975,
                           0.359724,
                           0.327071, 0.294045, 0.259905, 0.222172, 0.18875, 0.15242, 0.120852, 0.0919143, 0.0622389,
                           0.0316789,
                           0.00160561, -0.0315469, -0.0628154, -0.0899432, -0.118559, -0.145458, -0.173643, -0.205532,
                           -0.236153, -0.260781,
                           -0.285765, -0.310873, -0.335457, -0.363311, -0.390157, -0.414493, -0.436925, -0.462431,
                           -0.488903, -0.518141,
                           -0.543916, -0.564566, -0.583791, -0.606335, -0.631028, -0.659311, -0.685798, -0.707176,
                           -0.730214, -0.753894,
                           -0.776487, -0.798325, -0.817541, -0.837389, -0.86271, -0.886906, -0.908754, -0.927517,
                           -0.94514, -0.962288,
                           -0.98413, -1.00636, -1.02599, -1.0492, -1.06794, -1.08414, -1.10186, -1.12213, -1.1424,
                           -1.16392,
                           -1.18125, -1.20071, -1.22033, -1.24009, -1.25914, -1.27765, -1.29749, -1.31712, -1.33935,
                           -1.36006,
                           -1.37792, -1.39782, -1.41187, -1.42653, -1.44253, -1.46037, -1.47889, -1.50016, -1.51232,
                           -1.52873,
                           -1.54191, -1.55834, -1.57326, -1.59346, -1.61403, -1.63503, -1.65719, -1.67541, -1.69187,
                           -1.7045,
                           -1.7241, -1.7407, -1.76292, -1.78104, -1.7984, -1.81121, -1.83015, -1.85384, -1.8699,
                           -1.88268};

    for (int i = 0; i < 150; i++) {
        m_featureRank->setLogOdds(i, rankLogOdds[i]);
    }

    m_featureMzPos->initialize(100);
    float mzPosLogOdds[] = {-1.34657, -1.34657, -1.34657, -1.34657, -1.37168, -1.31375, -1.30542, -1.23415, -1.11003,
                            -0.972022,
                            -0.83933, -0.710855, -0.713758, -0.718696, -0.692271, -0.640897, -0.552938, -0.499806,
                            -0.459431, -0.443607,
                            -0.413077, -0.387957, -0.310813, -0.24045, -0.179349, -0.146892, -0.129763, -0.127471,
                            -0.0766487, -0.0902686,
                            -0.0081166, -0.010562, -0.0261707, -0.0484099, -0.0616517, -0.114581, -0.0458415,
                            -0.0315522, -0.00296058, 0.0309577,
                            0.0819645, 0.126689, 0.156309, 0.154758, 0.127764, 0.0202552, -0.0550962, 0.146643,
                            0.234947, 0.24339,
                            0.162565, -0.0596549, -0.232069, -0.126285, -0.106357, -0.0617923, 0.0420685, 0.0405285,
                            -0.0501975, -0.0450146,
                            -0.0446294, -0.0163725, 0.0113579, 0.0231656, 0.0449489, 0.0680706, 0.114268, 0.176331,
                            0.215261, 0.259835,
                            0.329065, 0.375858, 0.42079, 0.432524, 0.427477, 0.447402, 0.463656, 0.471419, 0.493349,
                            0.518875,
                            0.526353, 0.521248, 0.506394, 0.527454, 0.558109, 0.577996, 0.554414, 0.534352, 0.48976,
                            0.564687,
                            0.364586, 0.00121955, -0.542152, -0.777207, -0.923613, -0.666413, -0.502393, -0.616991,
                            -0.839414, -0.902675};

    for (int i = 0; i < 100; i++) {
        m_featureMzPos->setLogOdds(i, mzPosLogOdds[i]);
    }

    m_sisters.push_back(-18.01057);
    m_sisters.push_back(-17.0027);
    m_sisters.push_back(-9.005285);
    m_sisters.push_back(-2.016);
    m_sisters.push_back(-1.008);
    m_sisters.push_back(1.008);
    m_sisters.push_back(2.016);
    m_sisters.push_back(9.005285);
    m_sisters.push_back(17.0027);
    m_sisters.push_back(18.01057);
    m_sisters.push_back(27.9949);
    m_sisters.push_back(57.02146);
    m_sisters.push_back(71.03711);
    m_sisters.push_back(79.966331);
    m_sisters.push_back(87.03202);
    m_sisters.push_back(97.05276);
    m_sisters.push_back(97.976896);
    m_sisters.push_back(99.06841);
    m_sisters.push_back(101.0476);
    m_sisters.push_back(103.0091);
    m_sisters.push_back(113.084);
    m_sisters.push_back(114.0429);
    m_sisters.push_back(115.0269);
    m_sisters.push_back(115.987461);
    m_sisters.push_back(128.0585);
    m_sisters.push_back(128.0949);
    m_sisters.push_back(129.0425);
    m_sisters.push_back(131.0404);
    m_sisters.push_back(137.0589);
    m_sisters.push_back(147.0684);
    m_sisters.push_back(156.1011);
    m_sisters.push_back(163.0633);
    m_sisters.push_back(186.0793);

    m_featureWithSister->initialize(33);
    float withSisterLogOdds[] = {0.350278, 0.164399, 0.228062, 0.462959, 1.06885, 0.654958, 0.145021, 0.502616, 0.61074,
                                 0.802457,
                                 0.200676, 0.22245, 0.195438, 0.0630093, 0.1806, 0.182611, 0.222371, 0.311623, 0.179564,
                                 0.0406221,
                                 0.455735, 0.396997, 0.320707, 0.214306, 0.383839, 0.383838, 0.416872, 0.334572,
                                 0.084313, 0.382453,
                                 0.130522, 0.241397, 0.216694};

    for (int i = 0; i < 33; i++) {
        m_featureWithSister->setLogOdds(i, withSisterLogOdds[i]);
    }

    m_featureWithoutSister->initialize(33);
    float withoutSisterLogOdds[] = {-0.0873323, -0.0334474, -0.0434743, -0.094576, -0.440145, -0.270833, -0.0292009,
                                    -0.0967725, -0.126143, -0.201188,
                                    -0.0294654, -0.032661, -0.027162, -0.00764404, -0.023798, -0.0251726, -0.0302876,
                                    -0.0457539, -0.0238243, -0.00466176,
                                    -0.0751136, -0.062578, -0.0470924, -0.027647, -0.0554484, -0.0554921, -0.0639662,
                                    -0.0461552, -0.00916711, -0.0518711,
                                    -0.0141299, -0.0275938, -0.0232646};

    for (int i = 0; i < 33; i++) {
        m_featureWithoutSister->setLogOdds(i, withoutSisterLogOdds[i]);
    }


    m_featureWithComplement->initialize(7);
    float withComplementLogOdds[] = {0.299837, 0.299837, 0.838681, 0.352363, 0.0629686, -0.308204, -0.686169};

    for (int i = 0; i < 7; i++) {
        m_featureWithComplement->setLogOdds(i, withComplementLogOdds[i]);
    }

    m_featureWithoutComplement->initialize(7);
    float withoutComplementLogOdds[] = {0.299837, 0.299837, -0.0782653, -0.115022, -0.109954, -0.17707, -0.150952};

    for (int i = 0; i < 7; i++) {
        m_featureWithoutComplement->setLogOdds(i, withoutComplementLogOdds[i]);
    }

    m_priorCoefficients.push_back(0.103068);
    m_priorCoefficients.push_back(-0.00883144);
    m_priorCoefficients.push_back(0.00468882);
    m_priorCoefficients.push_back(0.20002);
    m_priorCoefficients.push_back(0.0116696);
    m_priorCoefficients.push_back(2.86969e-05);
    m_priorCoefficients.push_back(0.0108838);
    m_priorCoefficients.push_back(0.0116174);
    m_priorCoefficients.push_back(0.00205072);
    m_priorCoefficients.push_back(0.000605753);
    m_priorCoefficients.push_back(-0.00115826);
    m_priorCoefficients.push_back(-0.00178328);
    m_priorCoefficients.push_back(-0.00427138);
    m_priorCoefficients.push_back(-0.00246449);
    m_priorCoefficients.push_back(-0.00351068);
    m_priorCoefficients.push_back(-0.00531488);
    m_priorCoefficients.push_back(0.0526398);

    m_ready = true;

}


void SpectraSTDenoiser::addTrainingSpectrum(SpectraSTPeakList *pl) {

    if (m_ready) return;

    m_numSpectra++;

    float cumInten = 0.0;
    unsigned int numSignal = 0;
    unsigned int numNoise = 0;

    unsigned int rawNumPeaks = pl->getNumPeaks(); // before simplify

    int precursorCharge = pl->getParentCharge();
    double precursorMass = pl->getParentMz() * (double) precursorCharge;

    vector<double> priorData(17, 0.0);
    priorData[0] = 1.0;
    priorData[1] = log((double) rawNumPeaks); // raw num peaks before simplify
    priorData[2] = log(pl->getPrecursorIntensity());
    priorData[3] = pl->calcXrea();
    priorData[4] = log(pl->getParentMz() * (double) precursorCharge);
    priorData[5] = pl->getMzRange();

    unsigned int totalNumSistersFound = 0;

    //  pl->simplify(150, 99999, false, false, false);

    unsigned int maxRank = (rawNumPeaks > 150 ? 150 : rawNumPeaks);

    Peak thisPeak;
    pl->getNthLargestPeak(maxRank, thisPeak);

    double intensityAtMaxRank = thisPeak.intensity;

    Peak smallerPeak = thisPeak;

    Peak biggerPeak;

    for (unsigned int rank = maxRank; rank >= 1; rank--) {

        if (rank > 1) pl->getNthLargestPeak(rank - 1, biggerPeak);

        bool isSignal = false;
        if (thisPeak.info == "S") isSignal = true;

        // calculate CNI
        cumInten += thisPeak.intensity;

        m_numPeaksAtRank[rank - 1]++;
        m_sumIntensityAtRank[rank - 1] += thisPeak.intensity;
        m_sumSqIntensityAtRank[rank - 1] += (thisPeak.intensity * thisPeak.intensity);

        float cni = cumInten / pl->getTotalIonCurrent();
        int cniBin = (int) (cni * 100.0);

        // sum of top xx intense peaks' intensity and number of peaks the sum of whose CNI is within the certain range of TIC; e.g. 0-0.05
        if (rank <= 10) priorData[6] += cni;
        if (rank <= 20) priorData[7] += cni;
        if (rank <= 30) priorData[8] += cni;
        if (rank <= 40) priorData[9] += cni;
        if (rank <= 50) priorData[10] += cni;
        if (0 <= cni && cni < 0.05) priorData[11]++;
        if (0.05 <= cni && cni < 0.2) priorData[12]++;
        if (0.2 <= cni && cni < 0.35) priorData[13]++;
        if (0.35 <= cni && cni < 0.5) priorData[14]++;
        if (cni >= 0.5) priorData[15]++;


        // calculate rank
        int rankBin = rank - 1;
        if (rankBin > m_featureRank->getNumBins() - 1) rankBin = m_featureRank->getNumBins() - 1;


        /*
        // calculate deltaInt
        double meanIntAtRank = 10000.0;
        double stDevIntAtRank = 10000.0;
        if (rank > 1) {
          meanIntAtRank = 15000.0 * pow((double)rank, -0.85);
          stDevIntAtRank = meanIntAtRank * (0.14 * log((double)rank) + 0.22);
        }
        double deltaInt = (thisPeak.intensity - meanIntAtRank) / stDevIntAtRank;
        if (deltaInt > 2.999) deltaInt = 2.999;
        if (deltaInt < -2.999) deltaInt = -2.999;
        int deltaIntBin = (int)((deltaInt + 3.0) / 6.0 * (double)(m_featureDeltaInt->getNumBins()));
        */


        // calculate m/z position
        int mzPosBin = (int) (thisPeak.mz / precursorMass * (double) (m_featureMzPos->getNumBins()));

        if (isSignal) {
            numSignal++; // for this spectrum
            m_numSignal++; // sum of all spectra processed by this denoiser
            m_featureRank->incrementSignal(rankBin);
            //    m_featureDeltaInt->incrementSignal(deltaIntBin);
            m_featureMzPos->incrementSignal(mzPosBin);
        } else {
            numNoise++; // for this spectrum
            m_numNoise++; // sum of all spectra processed by this denoiser
            m_featureRank->incrementNoise(rankBin);
            //    m_featureDeltaInt->incrementNoise(deltaIntBin);
            m_featureMzPos->incrementNoise(mzPosBin);
        }

        // find complement
        // complementBin = thisCharge + compCharge (up to 6)
        for (int thisCharge = 1; thisCharge == 1 || thisCharge < precursorCharge; thisCharge++) {

            if (thisPeak.mz > (precursorMass / (double) thisCharge)) continue;

            for (int compCharge = 1; compCharge == 1 || compCharge < precursorCharge; compCharge++) {
                int numLostProtons = precursorCharge - thisCharge - compCharge;
                double compMz = (precursorMass - thisPeak.mz * (double) thisCharge -
                                 (double) numLostProtons * ((*Peptide::AAMonoisotopicMassTable)['+'])) /
                                (double) compCharge;
                if (compMz > (precursorMass / (double) compCharge)) continue;
                double complement = pl->findPeak(compMz, 0.25);

                unsigned int complementBin = thisCharge + compCharge;
                if (complementBin > 6) complementBin = 6;

                if (complement >= intensityAtMaxRank) {

                    if (isSignal) {
                        m_featureWithComplement->incrementSignal(complementBin);
                    } else {
                        m_featureWithComplement->incrementNoise(complementBin);
                    }

                } else {
                    // can't find a complement
                    if (isSignal) {
                        m_featureWithoutComplement->incrementSignal(complementBin);
                    } else {
                        m_featureWithoutComplement->incrementNoise(complementBin);
                    }
                }
            }
        }

        // find sisters
        for (unsigned int sisIndex = 0; sisIndex < (unsigned int) (m_sisters.size()); sisIndex++) {
            double mzDiff = m_sisters[sisIndex];
            if (pl->findPeak(thisPeak.mz - mzDiff, 0.25) >= intensityAtMaxRank) {

                totalNumSistersFound++;

                if (isSignal) {
                    m_featureWithSister->incrementSignal(sisIndex);
                } else {
                    m_featureWithSister->incrementNoise(sisIndex);
                }
            } else {
                if (isSignal) {
                    m_featureWithoutSister->incrementSignal(sisIndex);
                } else {
                    m_featureWithoutSister->incrementNoise(sisIndex);
                }
            }

        }

        smallerPeak = thisPeak;
        thisPeak = biggerPeak;

    } // END for all peaks

//  numNoise = rawNumPeaks - numSignal;
//  m_numNoise += numNoise;

//  priorData[16] = totalNumSistersFound / ((double)rawNumPeaks);
    priorData[16] = totalNumSistersFound / ((double) rawNumPeaks);

    for (vector<double>::iterator p = priorData.begin(); p != priorData.end(); p++) {
        cerr << *p << "  ";
    }

//  double prior = (double)numSignal / ((double)rawNumPeaks);
    double prior = (double) numSignal / ((double) rawNumPeaks);

    cerr << ": " << prior << endl;


    m_priors.push_back(prior);
    m_priorRegressionData.push_back(priorData);


}

void SpectraSTDenoiser::generateBayesianModel() {

    if (m_ready) return;

    if (m_numSignal < 10000 || m_numNoise < 10000) {
        // not enough data points
        m_ready = false;
        return;
    }

    // m_defaultPrior = (double)m_numSignal / m_numSpectra;
    m_defaultPrior = (double) m_numSignal / ((double) m_numSignal + (double) m_numNoise);

    double logOddsPrior = log((double) m_numSignal) - log((double) m_numNoise);

    // cout << "Rank\tMean\tStDev" << endl;
    /*
    for (unsigned int rank = 1; rank <= 200; rank++) {
      double mean = m_sumIntensityAtRank[rank - 1] / (double)(m_numPeaksAtRank[rank - 1]);
      double stdev = sqrt(m_sumSqIntensityAtRank[rank - 1] / (double)(m_numPeaksAtRank[rank - 1]) - mean * mean);
      //  cout << rank << "\t" << mean << "\t" << stdev << endl;
    }
    */

    cout << "Generating Bayesian Model -- S=" << m_numSignal << "; N=" << m_numNoise << "; logO=" << logOddsPrior
         << endl;

    m_featureRank->calculateLogOdds(logOddsPrior, 500);
//  m_featureDeltaInt->calculateLogOdds(logOddsPrior, 500);
    m_featureMzPos->calculateLogOdds(logOddsPrior, 500);
    m_featureWithSister->calculateLogOdds(logOddsPrior);
    m_featureWithoutSister->calculateLogOdds(logOddsPrior);
    m_featureWithComplement->calculateLogOdds(logOddsPrior);
    m_featureWithoutComplement->calculateLogOdds(logOddsPrior);

    m_featureRank->smoothLogOdds();
//  m_featureDeltaInt->smoothLogOdds();
    m_featureMzPos->smoothLogOdds();

#ifdef __LGPL__

    // LGPL -- no regression available. Just use the observed prior -- put it in the constant term and set all other
    // coefficients to zero
    m_priorCoefficients.assign(17, 0.0);
    m_priorCoefficients[0] = m_defaultPrior;

#endif

#ifndef __LGPL__
    gsl_multifit_linear_workspace *priorFittingWorkSpace = gsl_multifit_linear_alloc(m_priors.size(), 17);
    gsl_vector *gslPriors = gsl_vector_alloc(m_priors.size());
    gsl_matrix *gslDataMatrix = gsl_matrix_alloc(m_priors.size(), 17);

    for (int i = 0; i < m_priors.size(); i++) {
        gsl_vector_set(gslPriors, i, m_priors[i]);
        for (int j = 0; j < 17; j++) {
            gsl_matrix_set(gslDataMatrix, i, j, m_priorRegressionData[i][j]);
        }
    }

    gsl_vector *gslRegressionCoefficients = gsl_vector_alloc(17);
    gsl_vector_set_zero(gslRegressionCoefficients);
    gsl_matrix *gslRegressionCovMatrix = gsl_matrix_alloc(17, 17);
    gsl_matrix_set_zero(gslRegressionCovMatrix);

    double chisq;
    double tolerance = 0.01;
    size_t rank = 0;

    int errorCode = gsl_multifit_linear_svd(gslDataMatrix, gslPriors, tolerance, &rank, gslRegressionCoefficients,
                                            gslRegressionCovMatrix, &chisq, priorFittingWorkSpace);

    //cout << "chisq:" << chisq <<endl;
    double tss = gsl_stats_tss(m_priors.data(), 1, m_priors.size());

    double Rsq = 1 - chisq / tss;

    for (unsigned int i = 0; i < m_priors.size(); i++) {
        gsl_vector *gslPriorData = gsl_vector_alloc(17);
        for (unsigned int j = 0; j < 17; j++) {
            gsl_vector_set(gslPriorData, j, m_priorRegressionData[i][j]);
        }
        double regressedPrior = 0.0;
        double regressedPriorError = 0.0;
        errorCode = gsl_multifit_linear_est(gslPriorData, gslRegressionCoefficients, gslRegressionCovMatrix,
                                            &regressedPrior, &regressedPriorError);
//    cerr << m_priors[i] << '\t' << regressedPrior << '\t' << regressedPriorError << endl;
        gsl_vector_free(gslPriorData);
    }

    cout << "CHISQ= " << chisq << "; RSQ=" << Rsq << endl;

    for (unsigned int i = 0; i < 17; i++) {
        m_priorCoefficients[i] = gsl_vector_get(gslRegressionCoefficients, i);
    }

    gsl_vector_free(gslPriors);
    gsl_matrix_free(gslDataMatrix);
    gsl_matrix_free(gslRegressionCovMatrix);
    gsl_multifit_linear_free(priorFittingWorkSpace);
#endif

    m_ready = true;

    //  printModel();

}

void SpectraSTDenoiser::readModel(ifstream &fin) {

    string line("");
    while (nextLine(fin, line, "END CHARGE")) {
        string::size_type pos = 0;
        string feature = nextToken(line, pos, pos, " \t\r\n", " \t");
        string numBinStr = nextToken(line, pos, pos, " \t\r\n", " \t");
        int numBins = atoi(numBinStr.c_str());

        if (numBins < 1) {
            cerr << "ERROR: Bad spectrast.denoiser file" << endl;
            exit(1);
        }

        if (feature == "PRIOR") {
            string value = nextToken(line, pos, pos, " \t\r\n", " \t");
            m_defaultPrior = atof(value.c_str());
        }

        if (feature == "RANK") {
            m_featureRank->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureRank->setLogOdds(i, atof(value.c_str()));
            }

//    } else if (feature == "DELTA_INT") {
//      m_featureDeltaInt->initialize(numBins);
//
//      for (int i = 0; i < numBins; i++) {
//        string value = nextToken(line, pos, pos, " \t\r\n", " \t");
//	m_featureDeltaInt->setLogOdds(i, atof(value.c_str()));
//      }

        } else if (feature == "MZ_POS") {
            m_featureMzPos->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureMzPos->setLogOdds(i, atof(value.c_str()));
            }

        } else if (feature == "SISTER_DELTA") {

            m_sisters.resize(0);
            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_sisters.push_back(atof(value.c_str()));
            }

        } else if (feature == "WITH_SISTER") {
            m_featureWithSister->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureWithSister->setLogOdds(i, atof(value.c_str()));
            }

        } else if (feature == "WITHOUT_SISTER") {
            m_featureWithoutSister->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureWithoutSister->setLogOdds(i, atof(value.c_str()));
            }

        } else if (feature == "WITH_COMPLEMENT") {
            m_featureWithComplement->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureWithComplement->setLogOdds(i, atof(value.c_str()));
            }

        } else if (feature == "WITHOUT_COMPLEMENT") {
            m_featureWithoutComplement->initialize(numBins);

            for (int i = 0; i < numBins; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_featureWithoutComplement->setLogOdds(i, atof(value.c_str()));
            }

        } else if (feature == "PRIOR_REGRESSION") {
            string value("");
            m_priorCoefficients.resize(0, 0.0);
            for (int i = 0; i < 17; i++) {
                string value = nextToken(line, pos, pos, " \t\r\n", " \t");
                m_priorCoefficients.push_back(atof(value.c_str()));
            }

        } else {
            // unknown feature, do nothing
        }

    }

    if (m_sisters.size() != m_featureWithSister->getNumBins() ||
        m_sisters.size() != m_featureWithoutSister->getNumBins()) {
        cerr << "Inconsistent number of sisters considered. Quit.";
        exit(1);
    }
    if (m_featureWithComplement->getNumBins() != m_featureWithoutComplement->getNumBins()) {
        cerr << "Inconsistent number of complement bins. Quit.";
        exit(1);
    }

    m_ready = true;

    // check if this is a valid model
    if (m_featureRank->getLogOdds(0) < 0.01) m_ready = false;

}


void SpectraSTDenoiser::printModel(ofstream &fout) {

//  cout << "CNI" << endl;

    fout << "PRIOR\t1\t" << m_defaultPrior << endl;

    fout << "RANK\t";
    m_featureRank->printLogOdds(fout);
    fout << endl;

//  fout << "DELTA_INT\t";
//  m_featureDeltaInt->printLogOdds(fout);
//  fout << endl;

    fout << "MZ_POS\t";
    m_featureMzPos->printLogOdds(fout);
    fout << endl;

    fout << "SISTER_DELTA\t" << m_sisters.size();
    for (unsigned int sisIndex = 0; sisIndex < m_sisters.size(); sisIndex++) {
        fout << '\t' << m_sisters[sisIndex];
    }
    fout << endl;

    fout << "WITH_SISTER\t";
    m_featureWithSister->printLogOdds(fout);
    fout << endl;

    fout << "WITHOUT_SISTER\t";
    m_featureWithoutSister->printLogOdds(fout);
    fout << endl;

    fout << "WITH_COMPLEMENT\t";
    m_featureWithComplement->printLogOdds(fout);
    fout << endl;

    fout << "WITHOUT_COMPLEMENT\t";
    m_featureWithoutComplement->printLogOdds(fout);
    fout << endl;

    fout << "PRIOR_REGRESSION\t" << 17;
    for (int i = 0; i < 17; i++) {
        fout << "\t" << m_priorCoefficients[i];
    }
    fout << endl;

}

void SpectraSTDenoiser::filter(SpectraSTPeakList *pl, unsigned int maxNumPeaks, double minSignalProb) {

    if (!m_ready) return;
    // calculate signal prob for all peaks

    float cumInten = 0.0;

    unsigned int rawNumPeaks = pl->getNumPeaks();

    int precursorCharge = pl->getParentCharge();
    double precursorMass = pl->getParentMz() * (double) precursorCharge;

    vector<double> priorData(17, 0.0);
    priorData[0] = 1.0;
    priorData[1] = log((double) rawNumPeaks); // raw num peaks before simplify
    priorData[2] = log(pl->getPrecursorIntensity());
    priorData[3] = pl->calcXrea();
    priorData[4] = log(pl->getParentMz() * (double) precursorCharge);
    priorData[5] = pl->getMzRange();

    unsigned int totalNumSistersFound = 0;

    //  pl->simplify(150, 99999, false, false, false);

    vector<pair<double, Peak> > peakOdds; // vector of (logOdds, Peak) pairs

    unsigned int maxRank = (rawNumPeaks > 150 ? 150 : rawNumPeaks);

    Peak thisPeak;
    pl->getNthLargestPeak(maxRank, thisPeak);

    double intensityAtMaxRank = thisPeak.intensity;

    Peak smallerPeak = thisPeak;

    Peak biggerPeak;

    for (int rank = maxRank; rank >= 1; rank--) {

        if (rank > 1) pl->getNthLargestPeak(rank - 1, biggerPeak);

        double logOdds = 0.0;

        // calculate CNI
        cumInten += thisPeak.intensity;
        float cni = cumInten / pl->getTotalIonCurrent();
        int cniBin = (int) (cni * 100.0);

        // sum of top xx intense peaks' intensity and number of peaks the sum of whose CNI is within the certain range of TIC; e.g. 0-0.05
        if (rank <= 10) priorData[6] += cni;
        if (rank <= 20) priorData[7] += cni;
        if (rank <= 30) priorData[8] += cni;
        if (rank <= 40) priorData[9] += cni;
        if (rank <= 50) priorData[10] += cni;
        if (0 <= cni && cni < 0.05) priorData[11]++;
        if (0.05 <= cni && cni < 0.2) priorData[12]++;
        if (0.2 <= cni && cni < 0.35) priorData[13]++;
        if (0.35 <= cni && cni < 0.5) priorData[14]++;
        if (cni >= 0.5) priorData[15]++;

        // calculate rank
        int rankBin = rank - 1;
        if (rankBin > m_featureRank->getNumBins() - 1) rankBin = m_featureRank->getNumBins() - 1;
        logOdds += m_featureRank->getLogOdds(rankBin);

        /*
        // calculate deltaInt
        double meanIntAtRank = 10000.0;
        double stDevIntAtRank = 10000.0;
        if (rank > 1) {
          meanIntAtRank = 15000.0 * pow((double)rank, -0.85);
          stDevIntAtRank = meanIntAtRank * (0.14 * log((double)rank) + 0.22);
        }
        double deltaInt = (thisPeak.intensity - meanIntAtRank) / stDevIntAtRank;
        if (deltaInt > 2.999) deltaInt = 2.999;
        if (deltaInt < -2.999) deltaInt = -2.999;
        int deltaIntBin = (int)((deltaInt + 3.0) / 6.0 * (double)(m_featureDeltaInt->getNumBins()));

        // float deltaInt = (biggerPeak.intensity - thisPeak.intensity) / (biggerPeak.intensity - smallerPeak.intensity);
        //int deltaIntBin = (int)(deltaInt * (double)(m_featureDeltaInt->getNumBins()));
        logOdds += m_featureDeltaInt->getLogOdds(deltaIntBin);
        */

        int precursorCharge = pl->getParentCharge();
        double precursorMass = pl->getParentMz() * (double) precursorCharge;

        // calculate m/z position
        int mzPosBin = (int) (thisPeak.mz / precursorMass * (double) (m_featureMzPos->getNumBins()));
        logOdds += m_featureMzPos->getLogOdds(mzPosBin);

        // find complement
        // complementBin = thisCharge + compCharge (up to 6)
        for (int thisCharge = 1; thisCharge == 1 || thisCharge < precursorCharge; thisCharge++) {

            if (thisPeak.mz > (precursorMass / (double) thisCharge)) continue;

            for (int compCharge = 1; compCharge == 1 || compCharge < precursorCharge; compCharge++) {
                int numLostProtons = precursorCharge - thisCharge - compCharge;
                double compMz = (precursorMass - thisPeak.mz * (double) thisCharge -
                                 (double) numLostProtons * ((*Peptide::AAMonoisotopicMassTable)['+'])) /
                                (double) compCharge;
                if (compMz > (precursorMass / (double) compCharge)) continue;
                double complement = pl->findPeak(compMz, 0.25);

                unsigned int complementBin = thisCharge + compCharge;
                if (complementBin > 6) complementBin = 6;

                if (complement >= intensityAtMaxRank) {
                    logOdds += m_featureWithComplement->getLogOdds(complementBin);
                } else {
                    logOdds += m_featureWithoutComplement->getLogOdds(complementBin);
                }
            }
        }

        // find sisters
        for (unsigned int sisIndex = 0; sisIndex < (unsigned int) (m_sisters.size()); sisIndex++) {
            double mzDiff = m_sisters[sisIndex];
            if (pl->findPeak(thisPeak.mz - mzDiff, 0.25) >= intensityAtMaxRank) {
                totalNumSistersFound++;
                logOdds += m_featureWithSister->getLogOdds(sisIndex);
            } else {
                logOdds += m_featureWithoutSister->getLogOdds(sisIndex);
            }
        }

        pair<double, Peak> op;
        op.first = logOdds;
        op.second = thisPeak;

        peakOdds.push_back(op);

        smallerPeak = thisPeak;
        thisPeak = biggerPeak;

    } // END for all peaks

//  priorData[16] = totalNumSistersFound / (double)rawNumPeaks;
    priorData[16] = totalNumSistersFound / (double) rawNumPeaks;

    // sort all peaks by odds
    sort(peakOdds.begin(), peakOdds.end(), SpectraSTDenoiser::sortPeaksByOddsDesc);

    double prior = predictPrior(priorData);

    double logOddsPrior = log(prior / (1.0 - prior));

    unsigned int numPeaksKept = 0;
    double expectedNumSignalsKept = 0.0;

//  unsigned int predictedNumSignals = (unsigned int)(prior * (double)rawNumPeaks);
    unsigned int predictedNumSignals = (unsigned int) (prior * (double) rawNumPeaks);

    if (predictedNumSignals < 30) predictedNumSignals = 30;

    // insert peaks back into new peak list
    pl->m_peaks.clear();

    for (vector<pair<double, Peak> >::iterator i = peakOdds.begin(); i != peakOdds.end(); i++) {

        if (numPeaksKept >= maxNumPeaks || expectedNumSignalsKept >= predictedNumSignals) break;

        double probability = exp(i->first + logOddsPrior);
        probability = probability / (probability + 1.0);
        if (probability < minSignalProb) break;

        pl->insert(i->second.mz, i->second.intensity, i->second.annotation, i->second.info);
        numPeaksKept++;
        expectedNumSignalsKept += probability;
    }

    sort(pl->m_peaks.begin(), pl->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
    pl->m_isSortedByMz = true;


    if (pl->m_intensityRanked) {
        delete (pl->m_intensityRanked);
        pl->m_intensityRanked = NULL;
    }

    if (pl->m_peakMap) {
        delete (pl->m_peakMap);
        pl->m_peakMap = NULL;
    }

    if (pl->m_bins) {
        delete (pl->m_bins);
        pl->m_bins = NULL;
    }

}

double SpectraSTDenoiser::predictPrior(vector<double> &priorData) {

    double prior = m_defaultPrior;

    double testPrior = 0.0;

    for (unsigned int i = 0; i < 17; i++) {
        testPrior += m_priorCoefficients[i] * priorData[i];
    }

    if (testPrior <= 1.0 && testPrior >= 0.0) prior = testPrior;
    return (prior);

}

bool SpectraSTDenoiser::sortPeaksByOddsDesc(pair<double, Peak> a, pair<double, Peak> b) {

    return (a.first > b.first);

}
