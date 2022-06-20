#ifndef SPECTRASTPEAKLIST_HPP_
#define SPECTRASTPEAKLIST_HPP_

#include "SpectraSTSearchParams.hpp"
#include "SpectraSTCreateParams.hpp"
#include "SpectraSTSimScores.hpp"
#include "Peptide.hpp"

#include <string>
#include <vector>
#include <map>

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/

/* Class: SpectraSTPeakList
 * 
 * Implements a peak list for both a library spectrum and a query spectrum. This is where
 * all of the spectrum processing, filtering and comparing takes place. 
 */

using namespace std;

// a peak - note that we use float's for intensities to save memory
typedef struct _peak {
    double mz;
    float intensity;
    string annotation;
    string info;

} Peak;

class SpectraSTDenoiser;

class SpectraSTPeakList {

    friend class SpectraSTDenoiser;

public:

    // Constructor - general
    SpectraSTPeakList(double parentMz, int parentCharge, unsigned int numPeaks = 0, bool useBinIndex = false,
                      string fragType = "");

    // Constructor for creating consensus spectrum
    SpectraSTPeakList(vector<SpectraSTPeakList *> &pls, Peptide *pep, unsigned int totalNumRep, double quorum,
                      unsigned int maxNumPeaks, vector<SpectraSTDenoiser *> *denoisers = NULL,
                      bool keepRawIntensities = false);

    // Constructor for creating theoretical spectrum
    SpectraSTPeakList(Peptide *pep, string fragType = "");

    // Copy constructor
    SpectraSTPeakList(SpectraSTPeakList &other);

    SpectraSTPeakList &operator=(SpectraSTPeakList &other);

    // destructor
    virtual ~SpectraSTPeakList();

    // accessors
    unsigned int getNumPeaks();

    double getWeight() { return (m_weight); }

    double getMaxMz();

    double getMinMz();

    double getMzRange() { return (getMaxMz() - getMinMz()); }

    double getParentMz() { return (m_parentMz); }

    int getParentCharge() { return (m_parentCharge); }

    float getOrigMaxIntensity() { return (m_origMaxIntensity); }

    float getTotalIonCurrent() { return (m_totalIonCurrent); }

    float getPrecursorIntensity() { return (m_precursorIntensity); }

    unsigned int getNumAssignedPeaks() { return (m_numAssignedPeaks); }

    void getPeak(unsigned int index, Peak &peak) {
        peak = m_peaks[index];
        return;
    }

    void getNthLargestPeak(unsigned int n, Peak &p);

    string getFragType() { return (m_fragType); }

    bool isAnnotated() { return (m_isAnnotated); }

    string getFracUnassignedStr();

    // setters
    void insert(double mz, float intensity, string annotation, string info);

    void insertForSearch(double mz, float intensity, string annotation);

    void setWeight(double weight);

    void setNoiseFilterThreshold(double noiseFilterThreshold);

    void setPeptidePtr(Peptide *pep);

    void setPrecursorIntensity(float precursorIntensity) { m_precursorIntensity = precursorIntensity; }

    void setParentMz(double parentMz) { m_parentMz = parentMz; }

    void setParentCharge(int parentCharge, bool deleteBins = true);

    void setFragType(string &fragType);

    void setBinMagnitude(double value) { m_binMagnitude = value; }

    // comparing two peak lists by the dot product
    double compare(SpectraSTPeakList *other);

    double calcDot(SpectraSTPeakList *other);

    double calcDotAndDotBias(SpectraSTPeakList *other, double &dotBias);

    double calcDotNoBinning(SpectraSTPeakList *other, float mzTolerance);

    double calcDotTierwiseOpenModSearch(SpectraSTPeakList *other, float mzTolerance, pair<double, string> &openMod,
                                        int &numTiersUsed);

    // File output methods
    void writeToFile(ofstream &libFout);

    void writeToBinaryFile(ofstream &libFout);

    void writeToDtaFile(ofstream &dtaFout);

    // preprocessing methods
    bool passFilter(SpectraSTSearchParams &searchParams);

    bool passFilterUnidentified(SpectraSTCreateParams &params);

    void rankTransform(unsigned int maxRank, bool removePrecursor = false, double removeLightIonsMzCutoff = 0.0);

    void
    rankTransformWithQuota(unsigned int maxRank, bool removePrecursor = false, double removeLightIonsMzCutoff = 0.0,
                           int halfRange = 30, int quotaInRange = 4);

    void scalePeaks(double mzPower, double intensityPower, double unassignedFactor, bool rescale = false,
                    bool removePrecursor = false, double removeLightIonsMzCutoff = 0.0);

    void binPeaks(int numBinsPerMzUnit, double fractionToNeighbor, bool rebin = false);

    void binPeaksWithScaling(double mzPower, double intensityPower, double unassignedFactor,
                             int numBinsPerMzUnit, double fractionToNeighbor, bool rebin = false,
                             bool removePrecursor = false, double removeLightIonsMzCutoff = 0.0);

    void rankByIntensity(bool redo = false, unsigned int maxRank = 999999, bool removePrecursor = false,
                         double removeLightIonsMzCutoff = 0.0);

    void createPeakMap(bool redo = false);

    void normalizeTo(float basePeakValue, float maxDynamicRange = 0.0);

    void prepareForSearch(SpectraSTSearchParams &params, bool isLibrarySpectrum);

    void prepareForNoBinningDot();

    // peak list annotation method
    void annotate(bool redo = false, bool fixMz = false);

    // MRM table output
    void writeMRM(ofstream &fout, string pre, string post, string format);

    // peak list analysis methods
    double calcFractionUnassigned(unsigned int maxRank, unsigned int &numUnassigned, unsigned int &numAssigned,
                                  bool ignoreNearParentRegion = false, bool ignoreRareAnnotations = false);

    bool isSinglyCharged();

    double calcSignalToNoise();

    double calcXreaOld();

    double calcXrea(bool ignorePrecursorRegion = false);

    double calcXCorr();

    bool hasConsecutiveIonSeries();

    // peak list refining/cleaning methods
    void removeInquoratePeaks(unsigned int minNumRepWithPeak);

    void removeITRAQReporterPeaks(bool setToZeroOnly = false);

    void removeTMTReporterPeaks(bool setToZeroOnly = false);

    void quickSimplify(unsigned int maxNumPeaks, double maxDynamicRange, bool removePrecursor = false,
                       double removeLightIonsMzCutoff = 0.0);

    double
    simplify(unsigned int maxNumPeaks, double maxDynamicRange, bool deisotope = false, bool removePrecursor = false,
             double removeLightIonsMzCutoff = 0.0);

    double reduce(unsigned int maxNumPeaks, double minMz, double maxMz, unsigned int numRepsUsed = 1);

    double clean(bool deisotope, bool removePrecursor, double removeLightIonsMzCutoff);

    void centroid(string instrument);

    void removeNoncanonicalPeaks();

    void flattenAllPeaks();

    // peak-finding and peak-checking methods
    static bool isAssigned(string annotation, unsigned int NAA, bool ignoreRareAnnotations = false);

    static bool isBorY(string annotation);

    bool isNearPrecursor(double mz);

    double findPeak(double mz, double tolerance);

    Peak *findPeakPtr(double mz, double tolerance);

    Peak *findPeakPtr(double mz, double tolerance, map<double, double> *foundMZs);

    // convenient methods for debugging
    void printPeaks();

    void printBins();

    void plot(string name, string label);

    // decoy generation and semi-empirical prediction
    void repositionPeaks(bool keepEffectivePeakCountConstant = false);

    void shiftAllPeaks(double mzShift, double randomizeRange = 0.0);

    void useBinIndex();

    static void generatePoissonCutoffTable(int maxNumTrials, double alpha);

    static int findMinSignificantNumMatchedPeaks(int numBins, int numLibraryPeaks, int numQueryPeaks);

    static vector<double> *poissonCutoffTable;

private:


    double m_parentMz;

    int m_parentCharge;

    // m_peaks - the peaks. NOTE that there is no particular order maintained for it. The
    // peaks will just be in the order they are inserted into the vector
    vector<Peak> m_peaks;

    // m_bins - peaks will be put into bins for calculating the dot product.
    vector<float> *m_bins;

    // m_binIndex - when bin index is used, for storing the indices of the bins in m_bins
    vector<unsigned int> *m_binIndex;

    // m_intensityRanked - an index to m_peaks where the Peak pointers are sorted by decreasing intensity
    // for efficiency, this won't be instantiated at construction (since many operations on peak lists do not
    // require such a sorted list), but rather will only be created when rankByIntensity() is called.
    vector<Peak *> *m_intensityRanked;

    // m_peakMap - a map to m_peaks where the Peak pointers are keyed by m/z. This is used for quick binary search of
    // a peak of specified m/z.
    map<double, pair<Peak *, double> > *m_peakMap;

    // m_pep - points to the Peptide object that is the ID of this peak list (used only when the peak list is
    // a library spectrum, NULL otherwise). NOT the property of this class
    Peptide *m_pep;

    // m_fragType - the type of fragmentation (e.g. CID_IT, CID_Q, ETD)
    string m_fragType;

    // various fields to keep track of properties of the peak list
    unsigned int m_numBinsPerMzUnit;
    float m_binMagnitude;
    float m_peakMagnitude;
    double m_weight;
    double m_signalToNoise;
    float m_origMaxIntensity; // the "original" base peak intensity. this number is set when the peaks are inserted, and never changes
    // in any peak list manipulation.
    float m_totalIonCurrent; // the "original" total intensity. this number is set when the peaks are inserted, and never changes
    // in any peak list manipulation.
    double m_noiseFilterThreshold; // set this so that when inserting,
    // peaks smaller than threshold will not be inserted at all
    float m_precursorIntensity; // read off the mzXML file during pepXML input, not used typically

    // flags to make sure the peak list is not accidentally scaled or binned more than once.
    // NOTE: this is important, since the same library peak list object will be compared to multiple
    // query peak list objects in cached retrieval mode. This means that the first time the
    // library peak list object is compared, the preprocessing (scaling, binning, etc) is done, but
    // subsequently the same object will be compared to other query peak list objects without
    // the need to redo the preprocessing.

    bool m_isScaled;
    bool m_isAnnotated;
    bool m_isSortedByMz;
    bool m_isSearchOnly;

    unsigned int m_numAssignedPeaks;

    double m_scaleMzPower;
    double m_scaleIntensityPower;
    double m_scaleUnassignedFactor;

    unsigned int m_numBins;
    unsigned int m_numOccupiedBins;


    // helper methods
    unsigned int calcBinNumber(double mz);

    double calcBinMz(unsigned int binNum);

    float scale(Peak &p, double mzPower, double intensityPower, double unassignedFactor, bool removePrecursor,
                double removeLightIonsMzCutoff = 0.0);

    // annotation methods
    Peak *annotateIon(FragmentIon *fi, bool fixMz = false);

    void annotateIsotopicIons(Peak *assignedMonoisotopic, FragmentIon *fi, bool fixMz = false);

    void annotateInternalFragments(FragmentIon *fi, bool fixMz = false);

    map<double, pair<Peak *, double> >::iterator findBestPeakToAssign(double mzTh,
                                                                      map<double, pair<Peak *, double> >::iterator low,
                                                                      map<double, pair<Peak *, double> >::iterator high,
                                                                      double mzAccuracy);

    // consensus creation method
    void createConsensusSpectrum(vector<SpectraSTPeakList *> &pls, unsigned int totalNumRep, double quorum,
                                 unsigned int maxNumPeaks, vector<SpectraSTDenoiser *> *denoisers,
                                 bool keepRawIntensities);

    // theoretical spectrum creation method
    void generateTheoreticalSpectrum();

    // comparators for sorting
    static bool sortPeaksByIntensityDesc(Peak a, Peak b);

    static bool sortPeakPtrsByIntensityDesc(Peak *a, Peak *b);

    static bool sortPeaksByMzAsc(Peak a, Peak b);

    static bool sortFragmentIonsByProminence(FragmentIon a, FragmentIon b);

    static bool sortByMScoreDesc(pair<int, unsigned int> a, pair<int, unsigned int> b);
};

#endif /*SPECTRASTPEAKLIST_HPP_*/
