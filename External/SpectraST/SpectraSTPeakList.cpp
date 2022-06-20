
#include "SpectraSTPeakList.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTDenoiser.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>


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
 * all of the spectrum processing, filtering and comparing takes place. Can be subclassed
 * to use different spectrum processing, filtering and comparing methods.
 */



using namespace std;

extern SpectraSTLog *g_log;

vector<double> *SpectraSTPeakList::poissonCutoffTable = NULL;

// Constructor
SpectraSTPeakList::SpectraSTPeakList(double parentMz, int parentCharge, unsigned int numPeaks, bool useBinIndex,
                                     string fragType) :
        m_peaks(),
        m_parentMz(parentMz),
        m_parentCharge(parentCharge),
        m_pep(NULL),
        m_isScaled(false),
        m_isAnnotated(false),
        m_isSortedByMz(true),
        m_isSearchOnly(false),
        m_intensityRanked(NULL),
        m_peakMap(NULL),
        m_signalToNoise(-1.0),
        m_origMaxIntensity(-1.0),
        m_totalIonCurrent(0.0),
        m_precursorIntensity(0.0),
        m_noiseFilterThreshold(0.0),
        m_weight(1.0),
        m_binMagnitude(1.0),
        m_peakMagnitude(0.0),
        m_numAssignedPeaks(0),
        m_bins(NULL),
        m_binIndex(NULL),
        m_scaleMzPower(0.0),
        m_scaleIntensityPower(1.0),
        m_scaleUnassignedFactor(1.0),
        m_fragType("") {

    setFragType(fragType);

    if (numPeaks > 0) {
        m_peaks.reserve(numPeaks);
    }
    if (useBinIndex) {
        m_binIndex = new vector<unsigned int>;
    }

}

// constructor for consensus creation
SpectraSTPeakList::SpectraSTPeakList(vector<SpectraSTPeakList *> &pls, Peptide *pep, unsigned int totalNumRep,
                                     double quorum, unsigned int maxNumPeaks, vector<SpectraSTDenoiser *> *denoisers,
                                     bool keepRawIntensities) :
        m_peaks(),
        m_parentMz(0.0),
        m_parentCharge(0),
        m_pep(pep),
        m_isScaled(false),
        m_isAnnotated(false),
        m_isSortedByMz(true),
        m_isSearchOnly(false),
        m_intensityRanked(NULL),
        m_peakMap(NULL),
        m_signalToNoise(-1.0),
        m_origMaxIntensity(-1.0),
        m_totalIonCurrent(0.0),
        m_precursorIntensity(0.0),
        m_noiseFilterThreshold(0.0),
        m_weight(1.0),
        m_binMagnitude(1.0),
        m_peakMagnitude(0.0),
        m_numAssignedPeaks(0),
        m_bins(NULL),
        m_binIndex(NULL),
        m_scaleMzPower(0.0),
        m_scaleIntensityPower(1.0),
        m_scaleUnassignedFactor(1.0),
        m_fragType("") {

    if (pls.empty()) {
        // should not happen - caller should've checked this
        return;
    }

    m_parentMz = pls[0]->m_parentMz;
    m_parentCharge = pls[0]->m_parentCharge;

    setFragType(pls[0]->m_fragType);

    if (totalNumRep <= 0) {
        totalNumRep = (unsigned int) (pls.size());
    }

    createConsensusSpectrum(pls, totalNumRep, quorum, maxNumPeaks, denoisers, keepRawIntensities);

}

// constructor for theoretical spectrum generation from sequence
SpectraSTPeakList::SpectraSTPeakList(Peptide *pep, string fragType) :
        m_peaks(),
        m_parentMz(pep->monoisotopicMZ()),
        m_parentCharge(pep->charge),
        m_pep(pep),
        m_isScaled(false),
        m_isAnnotated(false),
        m_isSortedByMz(true),
        m_isSearchOnly(false),
        m_intensityRanked(NULL),
        m_peakMap(NULL),
        m_signalToNoise(-1.0),
        m_origMaxIntensity(-1.0),
        m_totalIonCurrent(0.0),
        m_noiseFilterThreshold(0.0),
        m_precursorIntensity(0.0),
        m_weight(1.0),
        m_binMagnitude(1.0),
        m_peakMagnitude(0.0),
        m_numAssignedPeaks(0),
        m_bins(NULL),
        m_binIndex(NULL),
        m_scaleMzPower(0.0),
        m_scaleIntensityPower(1.0),
        m_scaleUnassignedFactor(1.0),
        m_fragType("") {

    setFragType(fragType);
    generateTheoreticalSpectrum();

}

// Copy constructor
SpectraSTPeakList::SpectraSTPeakList(SpectraSTPeakList &other) :
        m_pep(NULL),
        m_intensityRanked(NULL),
        m_peakMap(NULL),
        m_bins(NULL),
        m_binIndex(NULL),
        m_scaleMzPower(0.0),
        m_scaleIntensityPower(1.0),
        m_scaleUnassignedFactor(1.0) {

    (*this) = other;
}

// assignment operator
SpectraSTPeakList &SpectraSTPeakList::operator=(SpectraSTPeakList &other) {

    this->m_parentMz = other.m_parentMz;
    this->m_parentCharge = other.m_parentCharge;
    this->m_pep = other.m_pep;
    this->m_isScaled = other.m_isScaled;
    this->m_isAnnotated = other.m_isAnnotated;
    this->m_isSortedByMz = other.m_isSortedByMz;
    this->m_isSearchOnly = other.m_isSearchOnly,
            this->m_totalIonCurrent = other.m_totalIonCurrent;
    this->m_precursorIntensity = other.m_precursorIntensity;
    this->m_noiseFilterThreshold = other.m_noiseFilterThreshold;
    this->m_weight = other.m_weight;
    this->m_numAssignedPeaks = other.m_numAssignedPeaks;
    this->m_signalToNoise = other.m_signalToNoise;
    this->m_origMaxIntensity = other.m_origMaxIntensity;
    this->m_binMagnitude = other.m_binMagnitude;
    this->m_peakMagnitude = other.m_peakMagnitude;
    this->m_fragType = other.m_fragType;

    if (other.m_isScaled) {
        this->m_scaleMzPower = other.m_scaleMzPower;
        this->m_scaleIntensityPower = other.m_scaleIntensityPower;
        this->m_scaleUnassignedFactor = other.m_scaleUnassignedFactor;
    }

    this->m_peaks.clear();
    for (vector<Peak>::iterator i = other.m_peaks.begin(); i != other.m_peaks.end(); i++) {
        m_peaks.push_back(*i);
    }

    if (other.m_intensityRanked) {
        rankByIntensity(true);
    } else {
        if (this->m_intensityRanked) delete (this->m_intensityRanked);
        this->m_intensityRanked = NULL;
    }

    if (other.m_peakMap) {
        createPeakMap(true);
    } else {
        if (this->m_peakMap) delete (this->m_peakMap);
        this->m_peakMap = NULL;
    }

    if (other.m_bins) {
        if (this->m_bins) delete (this->m_bins);
        this->m_bins = new vector<float>(*(other.m_bins));
    } else {
        if (this->m_bins) delete (this->m_bins);
        this->m_bins = NULL;
    }

    if (other.m_binIndex) {
        if (this->m_binIndex) delete (this->m_binIndex);
        this->m_binIndex = new vector<unsigned int>(*(other.m_binIndex));
    } else {
        if (this->m_binIndex) delete (this->m_binIndex);
        this->m_binIndex = NULL;
    }

    return (*this);
}

// destructor
SpectraSTPeakList::~SpectraSTPeakList() {

    if (m_intensityRanked) {
        delete (m_intensityRanked);
    }
    if (m_peakMap) {
        delete (m_peakMap);
    }
    if (m_bins) {
        delete (m_bins);
    }
    if (m_binIndex) {
        delete (m_binIndex);
    }

}

// getNumPeaks - returns the number of peaks
unsigned int SpectraSTPeakList::getNumPeaks() {
    return ((unsigned int) (m_peaks.size()));
}

// getMaxMz - returns the highest m/z -- this assumes the peaks are sorted by m/z, which they usually are!
double SpectraSTPeakList::getMaxMz() {
    if (m_peaks.empty()) return (0.0);

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    return (m_peaks[m_peaks.size() - 1].mz);
}

// getMinMz - returns the lowest m/z -- this assumes the peaks are sorted by m/z, which they usually are!
double SpectraSTPeakList::getMinMz() {
    if (m_peaks.empty()) return (0.0);

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    return (m_peaks[0].mz);
}


// insert - called to populate the peak list (e.g. as it is read from a file) -- this is the ONLY
// way to set the peaks (to make sure the peak list created is sound). NOTE: there is
// no particular order maintained for the m_peaks vector.
void SpectraSTPeakList::insert(double mz, float intensity, string annotation, string info) {

    if (intensity < (float) m_noiseFilterThreshold) {
        return;
    }

    Peak newPeak;
    newPeak.mz = mz;
    newPeak.intensity = intensity;
    newPeak.annotation = annotation;
    if (!annotation.empty()) {
        m_isAnnotated = true;
        if (annotation[0] != '?') {
            m_numAssignedPeaks++;
        }
    }

    newPeak.info = info;

    if (m_origMaxIntensity < intensity) {
        m_origMaxIntensity = intensity;
    }

    m_totalIonCurrent += intensity;

    if (!m_peaks.empty() && mz < m_peaks[m_peaks.size() - 1].mz) {
        m_isSortedByMz = false;
    }

    m_peaks.push_back(newPeak);

}

// insertForSearch - called to populate the peak list for SEARCH ONLY. Compared to insert,
// this function cuts a few corners to make searching more efficient.
void SpectraSTPeakList::insertForSearch(double mz, float intensity, string annotation) {

    if (intensity < (float) m_noiseFilterThreshold) {
        return;
    }

    Peak newPeak;
    newPeak.mz = mz;
    newPeak.intensity = intensity;
    newPeak.annotation = annotation;
    newPeak.info = "";

    if (m_origMaxIntensity < intensity) {
        m_origMaxIntensity = intensity;
    }

    m_peaks.push_back(newPeak);

}

// setNoiseFilterThreshold - sets the noise filter threshold (used in insert() to decide whether to insert a peak)
void SpectraSTPeakList::setNoiseFilterThreshold(double noiseFilterThreshold) {
    m_noiseFilterThreshold = noiseFilterThreshold;
}

// setPeptidePtr - sets a pointer to the Peptide object, which is NOT the property of SpectraSTPeakList!
void SpectraSTPeakList::setPeptidePtr(Peptide *pep) {
    m_pep = pep;
    // m_parentMz = pep->monoisotopicMZ();
    // m_parentCharge = pep->charge;
}

void SpectraSTPeakList::setFragType(string &fragType) {

    m_fragType = fragType;

}

// writeToFile - just writes to a file in .msp format (The caller needs to write
// the other headers (Name, Comments, etc) before calling this function).
void SpectraSTPeakList::writeToFile(ofstream &libFout) {

    libFout << "NumPeaks: " << m_peaks.size() << endl;

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    vector<Peak>::iterator i;

    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {
        Peak p = *i;

        libFout.precision(4);
        libFout << fixed << p.mz << "\t";
        libFout.precision(1);
        libFout << fixed << p.intensity << "\t";
        libFout << p.annotation << "\t";
        libFout << p.info << endl;
    }

}

// writeToBinaryFile - just writes to a file in binary.
void SpectraSTPeakList::writeToBinaryFile(ofstream &libFout) {

    unsigned int numPeaks = (unsigned int) (m_peaks.size());

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    libFout.write((char *) (&numPeaks), sizeof(unsigned int));

    vector<Peak>::iterator i;

    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {
        double mz = i->mz;
        double intensity = (double) (i->intensity);

        libFout.write((char *) &(mz), sizeof(double));
        libFout.write((char *) &(intensity), sizeof(double));

        libFout << i->annotation << endl;
        libFout << i->info << endl;

    }
}

// writeToDtaFile - similar to writeToFile, except that it won't write
// the Num peak: header, nor the annotations. (The caller needs to write the MW and
// charge before calling this function.)
void SpectraSTPeakList::writeToDtaFile(ofstream &dtaFout) {

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    vector<Peak>::iterator i;

    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {
        Peak p = *i;

        //if (p.annotation.find("p-98", 0) != string::npos ||
        //    p.annotation.find("p-116", 0) != string::npos ||
        //    p.annotation.find("p-196", 0) != string::npos) {
        //  continue;
        //}

        dtaFout.precision(4);
        dtaFout << fixed << p.mz << '\t';
        dtaFout.precision(1);
        dtaFout << fixed << p.intensity << endl;
    }

}

// compare - calculates the dot product of two spectra with the default parameters
double SpectraSTPeakList::compare(SpectraSTPeakList *other) {

    binPeaksWithScaling(0.0, 0.5, 1.0, 1, 0.5, false, true, 0.0);
    other->binPeaksWithScaling(0.0, 0.5, 1.0, 1, 0.5, false, true, 0.0);

    return (calcDot(other));

}

void SpectraSTPeakList::setParentCharge(int parentCharge, bool deleteBins) {

    if (parentCharge == m_parentCharge) return;

    m_parentCharge = parentCharge;

    if (deleteBins && m_bins) {
        delete (m_bins);
    }
}


double SpectraSTPeakList::calcDotNoBinning(SpectraSTPeakList *other, float mzTolerance) {

    if (this->m_peakMagnitude < 0.00001 || other->m_peakMagnitude < 0.00001) {
        return (0.0);
    }

    if (!m_isSortedByMz) {
        sort(this->m_peaks.begin(), this->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        this->m_isSortedByMz = true;
    }
    if (!(other->m_isSortedByMz)) {
        sort(other->m_peaks.begin(), other->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        other->m_isSortedByMz = true;
    }

    vector<Peak>::iterator i = this->m_peaks.begin();
    vector<Peak>::iterator j = other->m_peaks.begin();

    float dot = 0.0;

    while (i != this->m_peaks.end() && j != other->m_peaks.end()) {

        if (fabs(i->mz - j->mz) < mzTolerance) {

            dot += (i->intensity * j->intensity);
            i++;
            j++;

        } else if (i->mz < j->mz - mzTolerance) {

            i++;

        } else {

            j++;

        }

    }

    dot /= (m_peakMagnitude * other->m_peakMagnitude);

    return ((double) dot);
}

double SpectraSTPeakList::calcDotTierwiseOpenModSearch(SpectraSTPeakList *other, float mzTolerance,
                                                       pair<double, string> &openMod, int &numTiersUsed) {

    openMod.first = 0.0;
    openMod.second = "";

    if (this->m_peakMagnitude < 0.00001 || other->m_peakMagnitude < 0.00001) {
        return (0.0);
    }

    int thisParentCharge = this->m_parentCharge;
    if (thisParentCharge == 0) thisParentCharge = other->m_parentCharge;

    double thisParentMass = this->m_parentMz * (double) (thisParentCharge);
    double otherParentMass = other->m_parentMz * (double) (other->m_parentCharge);

    double deltaMass = thisParentMass - otherParentMass;
    double absDeltaMass = fabs(deltaMass);

    if (absDeltaMass >= 300.0) {
        return (0.0);
    }

    // mass defect filter
    double thisMassDefect = thisParentMass / 1.00048 - round(thisParentMass / 1.00048);
    double otherMassDefect = otherParentMass / 1.00048 - round(otherParentMass / 1.00048);
    double deltaMassDefect = thisMassDefect - otherMassDefect;
    if (absDeltaMass >= 3.0 && (deltaMassDefect > 0.15 || deltaMassDefect < -0.15)) {
        return (0.0);
    }

    // begin similarity calculation
    if (!m_isSortedByMz) {
        sort(this->m_peaks.begin(), this->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        this->m_isSortedByMz = true;
    }
    if (!(other->m_isSortedByMz)) {
        sort(other->m_peaks.begin(), other->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
        other->m_isSortedByMz = true;
    }

    int maxTier = other->m_parentCharge;
    int numBins = (int) (this->getMzRange() / (mzTolerance * 2.0));

    vector<float> tierDot(maxTier + 2, 0.0);
    vector<float> tierWeight(maxTier + 2, 0.0);

    vector<int> tierNumPeaks1(maxTier + 2, 0);
    vector<int> tierNumPeaks2(maxTier + 2, 0);

    vector<int> tierNumMatchedPeaks(maxTier + 2, 0);

    int naa = 0;
    if (other->m_pep && other->isAnnotated()) {
        // library spectrum is a peptide and is annotated
        naa = other->m_pep->NAA();
    }

    vector<int> bUnchanged(naa + 1, 0);
    vector<int> bChanged(naa + 1, 0);
    vector<int> yUnchanged(naa + 1, 0);
    vector<int> yChanged(naa + 1, 0);

    vector<Peak> thisMatched;
    vector<Peak> thisUnmatched;
    vector<Peak> otherMatched;
    vector<Peak> otherUnmatched;

    // tier 0, use original peaks
    vector<Peak>::iterator i = this->m_peaks.begin();
    vector<Peak>::iterator j = other->m_peaks.begin();

    int tier = 0;

    while (i != this->m_peaks.end() && j != other->m_peaks.end()) {

        if ((fabs(i->mz - j->mz) <= mzTolerance)) {

            tierDot[tier] += (i->intensity * j->intensity);
            thisMatched.push_back(*i);
            otherMatched.push_back(*j);
            tierNumMatchedPeaks[tier]++;

            i++;
            j++;

        } else if (i->mz - j->mz < -mzTolerance) {
            tierNumPeaks1[tier]++;
            thisUnmatched.push_back(*i);
            i++;

        } else {
            tierNumPeaks2[tier]++;
            otherUnmatched.push_back(*j);
            j++;

        }
    }

    while (i != this->m_peaks.end()) {
        tierNumPeaks1[tier]++;
        thisUnmatched.push_back(*i);
        i++;
    }

    while (j != other->m_peaks.end()) {
        tierNumPeaks2[tier]++;
        otherUnmatched.push_back(*j);
        j++;
    }

    // onto upper tiers NOTE: short-circuit if tier-0 dot is too low?
    for (tier = 1; tier < maxTier; tier++) {

        i = thisUnmatched.begin();
        j = otherUnmatched.begin();

        //calculate the m/z shift
        double mzShift = deltaMass / (double) tier;

        while (i != thisUnmatched.end() && j != otherUnmatched.end()) {

            if ((i->intensity < 0.01)) {
                i++;
                continue;
            }

            if (j->intensity < 0.01) {
                j++;
                continue;
            }

            if ((fabs(i->mz - j->mz - mzShift) <= mzTolerance)) {

                tierDot[tier] += (i->intensity * j->intensity);
                tierNumMatchedPeaks[tier]++;

                // delete these matched peaks -- they will not be matched again in upper tiers
                i->intensity = 0.0;
                j->intensity = 0.0;

                FragmentIon fi(j->mz, j->annotation, 0);
                // cerr << "Matched: " << "T" << tier << " " << j->annotation << "=" << fi.m_ion[0] << fi.m_pos << endl;

                if (fi.m_pos > 0 && !fi.m_bracket && fi.m_isotope == 0 && fi.m_loss == 0) {
                    if (fi.m_ion[0] == 'b') {
                        bChanged[fi.m_pos]++;
                    } else if (fi.m_ion[0] == 'y') {
                        yChanged[fi.m_pos]++;
                    }
                }

                i++;
                j++;

            } else if (i->mz - j->mz - mzShift < -mzTolerance) {
                tierNumPeaks1[tier]++;
                i++;

            } else {
                tierNumPeaks2[tier]++;
                j++;

            }
        }

        while (i != thisUnmatched.end()) {
            if (i->intensity > 0.01) tierNumPeaks1[tier]++;
            i++;
        }

        while (j != otherUnmatched.end()) {
            if (j->intensity > 0.01) tierNumPeaks2[tier]++;
            j++;
        }

    }

    numTiersUsed = 1; // always use Tier 0

    for (tier = 0; tier < maxTier; tier++) {

        int numMatchedPeaksCutoff = SpectraSTPeakList::findMinSignificantNumMatchedPeaks(numBins, tierNumPeaks2[tier],
                                                                                         tierNumPeaks1[tier]);

        // cerr << "T" << tier << ": N=" << numBins << ";M=" << tierNumPeaks2[tier] << ";n=" << tierNumPeaks1[tier] << ";k=" << tierNumMatchedPeaks[tier] << ";k_cutoff=" << numMatchedPeaksCutoff << endl;;

        if (tierNumMatchedPeaks[tier] >= numMatchedPeaksCutoff) {
            tierWeight[tier] = 1.0 + (tierNumMatchedPeaks[tier] - numMatchedPeaksCutoff) / 2.0;
            if (tier > 0) numTiersUsed++;
        } else {
            tierWeight[tier] = 0.0;
        }

    }

    if (tierWeight[0] < 0.0001) {
        tierWeight[0] = 1.0;
    }

    double maxTierWeight = 0.0;
    for (int tier = 0; tier < maxTier; tier++) {
        if (maxTierWeight < tierWeight[tier]) maxTierWeight = tierWeight[tier];
    }

    stringstream tss;

    float dot = 0.0;
    for (int tier = 0; tier < maxTier; tier++) {
        // cerr << "T" << tier << ": " << tierWeight[tier] << ", " << tierDot[tier] << endl;
        float weightedTierDot =
                tierWeight[tier] * tierDot[tier] / (this->m_peakMagnitude * other->m_peakMagnitude); // / maxTierWeight;
        dot += weightedTierDot;

        if (tier > 0) tss << '+';
        tss.precision(2);
        tss << fixed << weightedTierDot;
    }

    openMod.first = deltaMass;
    tss << '|';
    openMod.second = tss.str();

    // do localization if more than one tier used? Or dot passing a certain score?
    if (absDeltaMass > 3.0 && dot > 0.2 && numTiersUsed > 1) {

        // parse annotations of matched peaks at tier 0 (as evidence against mod location)
        for (vector<Peak>::iterator m = otherMatched.begin(); m != otherMatched.end(); m++) {
            FragmentIon fi(m->mz, m->annotation, 0);
            // cerr << "Matched: " << "T0" << " " << m->annotation << "=" << fi.m_ion[0] << fi.m_pos << endl;

            if (fi.m_pos > 0 && !fi.m_bracket && fi.m_isotope == 0 && fi.m_loss == 0) {
                if (fi.m_ion[0] == 'b') {
                    bUnchanged[fi.m_pos]++;
                } else if (fi.m_ion[0] == 'y') {
                    yUnchanged[fi.m_pos]++;
                }
            }
        }

        int maxLocalizationScore = -99999;
        vector<int> bestModPos;

        for (int modPos = 1; modPos <= naa; modPos++) {
            int evidenceFor = 0;
            int evidenceAgainst = 0;
            for (int pos1 = 1; pos1 < naa; pos1++) {
                if (pos1 >= modPos) {
                    evidenceFor += bChanged[pos1]; // cerr << "FOR yC" << naa - pos1 << ":" << yChanged[naa - pos1] << endl;
                    evidenceAgainst += bUnchanged[pos1]; // cerr << "AGT yU" << naa - pos1 << ":" << yUnchanged[naa - pos1] << endl;
                } else {
                    evidenceFor += bUnchanged[pos1];
                    evidenceAgainst += bChanged[pos1];
                }
            }
            for (int pos2 = 1; pos2 < naa; pos2++) {
                if (pos2 >= naa - modPos + 1) {
                    evidenceFor += yChanged[pos2]; // cerr << "FOR bC" << naa - pos2 << ":" << bChanged[naa - pos2] << endl;
                    evidenceAgainst += yUnchanged[pos2]; // cerr << "AGT bU" << pos2 << ":" << bUnchanged[pos2] << endl;
                } else {
                    evidenceFor += yUnchanged[pos2];
                    evidenceAgainst += yChanged[pos2];
                }
            }

            int localizationScore = evidenceFor - evidenceAgainst;
            // cerr << "Loc " << modPos << ": FOR=" << evidenceFor << ", AGT=" << evidenceAgainst << endl;
            if (localizationScore > maxLocalizationScore) {
                maxLocalizationScore = localizationScore;
                bestModPos.clear();
                bestModPos.push_back(modPos);
            } else if (localizationScore == maxLocalizationScore) {
                bestModPos.push_back(modPos);
            }

        }

        if (maxLocalizationScore > 0 && (!(bestModPos.empty()))) {

            stringstream mpss;
            mpss << bestModPos[0];
            for (int mp = 1; mp < bestModPos.size(); mp++) {
                mpss << '/' << bestModPos[mp];
            }
            openMod.second += mpss.str();

        }


    }


    return (dot);

}


// calcDot - the dot product calculation method, without also calculating dot bias. See
// calcDotAndDotBias for more documentation.
double SpectraSTPeakList::calcDot(SpectraSTPeakList *other) {

    if ((!(this->m_bins)) || (!(other->m_bins))) {
        return (0.0);
    }
    if (other->m_binMagnitude < 0.00001 || m_binMagnitude < 0.00001) {
        return (0.0);
    }

    float d = 0.0;
    float dot = 0.0;

    if (!(this->m_binIndex) && (!(other->m_binIndex))) {

        vector<float>::iterator i;
        vector<float>::iterator j;
        for (i = this->m_bins->begin(), j = other->m_bins->begin();
             i != this->m_bins->end() && j != other->m_bins->end();
             i++, j++) {
            d = (*i) * (*j);
            dot += d;
        }

    } else if ((!(this->m_binIndex)) && other->m_binIndex) {

        vector<float>::iterator j;
        vector<unsigned int>::iterator jj;

        for (j = other->m_bins->begin(), jj = other->m_binIndex->begin();
             j != other->m_bins->end() && jj != other->m_binIndex->end();
             j++, jj++) {

            d = (*j) * ((*(this->m_bins))[*jj]);
            dot += d;
        }

    } else if (this->m_binIndex && (!(other->m_binIndex))) {

        vector<float>::iterator i;
        vector<unsigned int>::iterator ii;

        for (i = this->m_bins->begin(), ii = this->m_binIndex->begin();
             i != this->m_bins->end() && ii != this->m_binIndex->end();
             i++, ii++) {

            d = (*i) * ((*(other->m_bins))[*ii]);
            dot += d;
        }

    } else {

        vector<float>::iterator i = this->m_bins->begin();
        vector<unsigned int>::iterator ii = this->m_binIndex->begin();
        vector<float>::iterator j = other->m_bins->begin();
        vector<unsigned int>::iterator jj = other->m_binIndex->begin();

        while (ii != this->m_binIndex->end() && jj != other->m_binIndex->end()) {
            if (*ii == *jj) {
                d = (*i) * (*j);
                dot += d;
                i++;
                ii++;
                j++;
                jj++;
            } else if (*ii < *jj) {
                i++;
                ii++;
            } else {
                j++;
                jj++;
            }
        }
    }


    // normalize to 1 by dividing by the magnitudes
    return ((double) (dot / (m_binMagnitude * other->m_binMagnitude)));
}


// calcDotAndDotBias - calculates the dot product (i.e. cos(theta) where theta is the angle
// between the spectral vectors, and the dot bias. Note that the bins are dotted, not the individual peaks.
double SpectraSTPeakList::calcDotAndDotBias(SpectraSTPeakList *other, double &dotBias) {

    if ((!(this->m_bins)) || (!(other->m_bins))) {
        dotBias = 0.0;
        return (0.0);
    }
    if (other->m_binMagnitude < 0.00001 || m_binMagnitude < 0.00001) {
        dotBias = 0.0;
        return (0.0);
    }

    float d = 0.0;
    float dot = 0.0;
    float sumDotSquares = 0.0;

    // 4 cases: either of the two peak lists can use binIndex or not

    if (!(this->m_binIndex) && (!(other->m_binIndex))) {
        // both don't use a binIndex - easiest, just dot the corresponding bins

        vector<float>::iterator i;
        vector<float>::iterator j;
        for (i = this->m_bins->begin(), j = other->m_bins->begin();
             i != this->m_bins->end() && j != other->m_bins->end();
             i++, j++) {

            d = (*i) * (*j);
            dot += d;
            sumDotSquares += d * d;
        }

    } else if ((!(this->m_binIndex)) && other->m_binIndex) {
        // other uses a binIndex. In this case, go down other->m_binIndex, and use
        // the m/z indices to index into this's bins

        vector<float>::iterator j;
        vector<unsigned int>::iterator jj;

        for (j = other->m_bins->begin(), jj = other->m_binIndex->begin();
             j != other->m_bins->end() && jj != other->m_binIndex->end();
             j++, jj++) {

            d = (*j) * ((*(this->m_bins))[*jj]);
            dot += d;
            sumDotSquares += d * d;
        }

    } else if (this->m_binIndex && (!(other->m_binIndex))) {
        // this uses a binIndex. In this case, go down this->m_binIndex, and use
        // the m/z indices to index into other's bins


        vector<float>::iterator i;
        vector<unsigned int>::iterator ii;

        for (i = this->m_bins->begin(), ii = this->m_binIndex->begin();
             i != this->m_bins->end() && ii != this->m_binIndex->end();
             i++, ii++) {

            d = (*i) * ((*(other->m_bins))[*ii]);
            dot += d;
            sumDotSquares += d * d;
        }

    } else {
        // both uses a binIndex. Then have to do it the see-saw way...

        vector<float>::iterator i = this->m_bins->begin();
        vector<unsigned int>::iterator ii = this->m_binIndex->begin();
        vector<float>::iterator j = other->m_bins->begin();
        vector<unsigned int>::iterator jj = other->m_binIndex->begin();

        while (ii != this->m_binIndex->end() && jj != other->m_binIndex->end()) {
            if (*ii == *jj) {
                d = (*i) * (*j);
                dot += d;
                sumDotSquares += d * d;

                i++;
                ii++;
                j++;
                jj++;
            } else if (*ii < *jj) {
                i++;
                ii++;
            } else {
                j++;
                jj++;
            }
        }
    }

    // normalize to 1 by dividing by the magnitudes
    sumDotSquares /= (double) ((m_binMagnitude * m_binMagnitude * other->m_binMagnitude * other->m_binMagnitude));
    dot /= ((double) (m_binMagnitude * other->m_binMagnitude));
    if (sumDotSquares < 0.0001 || dot < 0.01) {
        dotBias = 0.0;
    } else {
        dotBias = sqrt(sumDotSquares) / dot;
    }
    return (dot);
}

// printPeaks - just cout all the peaks, for debugging only.
void SpectraSTPeakList::printPeaks() {

    vector<Peak>::iterator i;
    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {
        cout << i->mz << '\t' << i->intensity << '\t' << i->annotation << '\t' << i->info << endl;
    }
    cout << "DONE" << endl;

}

// printPeaks - just cout all the bins, for debugging only.
void SpectraSTPeakList::printBins() {

    if (!m_bins) {
        cout << "NOT BINNED!" << endl;
        return;
    }

    if (m_binIndex) {
        vector<float>::iterator i;
        vector<unsigned int>::iterator ii;

        for (ii = m_binIndex->begin(), i = m_bins->begin();
             ii != m_binIndex->end() && i != m_bins->end();
             ii++, i++) {

            cout << (*ii) << "\t" << calcBinMz((unsigned int) *ii) << "\t" << (*i) << endl;
        }

    } else {

        vector<float>::iterator i;

        int binNum = 0;
        for (i = m_bins->begin(); i != m_bins->end(); i++) {
            cout << binNum << "\t" << calcBinMz(binNum) << "\t" << (*i) << endl;
            binNum++;
        }

    }
    cout << "DONE" << endl;

}

// writeMRM - writes peak list as MRM transitions
void SpectraSTPeakList::writeMRM(ofstream &fout, string pre, string post, string format) {

    if (format == "DEFAULT") {
        rankByIntensity();

        for (vector<Peak *>::iterator i = m_intensityRanked->begin(); i != m_intensityRanked->end(); i++) {
            fout << pre << '\t';
            fout.precision(4);
            fout << fixed << (*i)->mz << '\t';
            fout.precision(2);
            fout << fixed << (*i)->intensity << '\t';
            fout << (*i)->annotation << '\t' << post << endl;
        }
    } else if (format == "SHOWINFO") {
        rankByIntensity();
        for (vector<Peak *>::iterator i = m_intensityRanked->begin(); i != m_intensityRanked->end(); i++) {
            fout << pre << '\t';
            fout.precision(4);
            fout << fixed << (*i)->mz << '\t';
            fout.precision(2);
            fout << fixed << (*i)->intensity << '\t';
            fout << (*i)->annotation << '\t' << post << '\t' << (*i)->info << endl;
        }
    }

}


// scalePeaks - scale peaks. WARNING: using this method will change the intensities in the m_peaks array
// Use binPeaks to scale and bin if keeping record of the original intensities is desired.
void SpectraSTPeakList::scalePeaks(double mzPower, double intensityPower, double unassignedFactor,
                                   bool rescale, bool removePrecursors, double removeLightIonsMzCutoff) {

    if (m_isScaled && !rescale) {
        return;
    }

    m_scaleMzPower = mzPower;
    m_scaleIntensityPower = intensityPower;
    m_scaleUnassignedFactor = unassignedFactor;

    vector<Peak>::iterator i;

    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {
        i->intensity = scale(*i, mzPower, intensityPower, unassignedFactor, removePrecursors, removeLightIonsMzCutoff);
    }

    m_isScaled = true;
}

// rankTransform
void SpectraSTPeakList::rankTransform(unsigned int maxRank, bool removePrecursor, double removeLightIonsMzCutoff) {

    rankByIntensity(true, maxRank + 1, removePrecursor, removeLightIonsMzCutoff);

    int numPeaksKept = (maxRank <= m_intensityRanked->size() ? maxRank : m_intensityRanked->size());

    vector<Peak *> retained;

    for (int rank = 1; rank <= m_intensityRanked->size(); rank++) {
        Peak *p = (*m_intensityRanked)[rank - 1];

        if (rank > numPeaksKept) {
            p->intensity = 0.0;
        } else {
            // p->intensity = 10000.0 - 9999.0 * (rank - 1) / ((float)maxRank - 1);
            retained.push_back(p);
            // p->intensity = 10000.0 - 10000.0 * (rank - 1) / ((float)maxRank);
            //      p->intensity = (float)(numPeaksKept - rank + 1) * 200.0; // the 200.0 multipler fits the expected scale in plotspectrast
        }
    }

    double curIntensity = 10000.0;
    double step = 10000.0 / retained.size();
    for (vector<Peak *>::iterator r = retained.begin(); r != retained.end(); r++) {
        (*r)->intensity = curIntensity;
        curIntensity -= step;
    }

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    m_isScaled = true;

}

void SpectraSTPeakList::rankTransformWithQuota(unsigned int maxRank, bool removePrecursor,
                                               double removeLightIonsMzCutoff, int halfRange, int quotaInRange) {

    rankByIntensity(true, 9999999, removePrecursor, removeLightIonsMzCutoff);

    int numPeaksToKeep = (maxRank <= m_intensityRanked->size() ? maxRank : m_intensityRanked->size());

    int numPeaks = 0;
    //int halfRange = 50;
    // int quotaInRange = 5;

    vector<int> quotas(MAX_MZ + 3, quotaInRange);

    vector<Peak *> retained;

    for (int rank = 1; rank <= m_intensityRanked->size(); rank++) {
        Peak *p = (*m_intensityRanked)[rank - 1];
        int location = (int) (p->mz + 0.5);
        if (location > MAX_MZ) location = MAX_MZ;

        if (numPeaks > numPeaksToKeep || quotas[location] <= 0) {
            p->intensity = 0.0;
        } else {

            retained.push_back(p);
            // p->intensity = 10000.0 - 10000.0 * (numPeaks) / ((float)numPeaksToKeep);
            for (int loc = location - halfRange; loc <= location + halfRange; loc++) {
                if (loc > 0 && loc <= MAX_MZ) quotas[loc]--; // within +/-50 Th of this peak, decrement quota by one
            }
            quotas[location] = 0; // block this location (1-mz wide), no other peak allowed
            quotas[location + 1] = 0; // block any +1 isotopic peaks
            quotas[location + 2] = 0; // block any +2 isotopic peaks
            numPeaks++;

        }
    }

    double curIntensity = 10000.0;
    double step = 10000.0 / retained.size();
    for (vector<Peak *>::iterator r = retained.begin(); r != retained.end(); r++) {
        (*r)->intensity = curIntensity;
        curIntensity -= step;
    }

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    m_isScaled = true;

}


// scale - scale peaks.
float
SpectraSTPeakList::scale(Peak &p, double mzPower, double intensityPower, double unassignedFactor, bool removePrecursor,
                         double removeLightIonsMzCutoff) {

    if (m_isScaled) return (p.intensity); // already scaled!

    if (p.intensity < 0.00001) return (0.0);

    if (removePrecursor && isNearPrecursor(p.mz)) {
        return (0.0);
    }

    if (p.mz < removeLightIonsMzCutoff) {
        return (0.0);
    }

    double scaledIntensity = pow((double) (p.mz), mzPower) * pow((double) (p.intensity), intensityPower);
    double factor = 1.0;
    if (!(p.annotation.empty()) && p.annotation[0] == '?') {
        factor = unassignedFactor;
    }
    scaledIntensity *= factor;

    return ((float) scaledIntensity);

}

bool SpectraSTPeakList::isNearPrecursor(double mz) {

    if (m_parentMz < 0.0001) return (false); // m_parentMz not set, can't determine

    if (m_fragType != "ETD" && m_fragType != "ETD-SA" && m_fragType != "ETD-HR") {
        if (mz > m_parentMz - 60.0 && mz < m_parentMz + 20.0) {
            return (true);
        }
    } else {
        int lowCharge = m_parentCharge;
        int highCharge = m_parentCharge;
        //    if (m_parentCharge == 0) {
        // remove all possible charge-reduced precursors for precursor charge up to 6.
        lowCharge = 3;
        highCharge = 6;
        //}

        for (int parentCharge = lowCharge; parentCharge <= highCharge; parentCharge++) {
            double parentMass = m_parentMz * parentCharge;
            for (int charge = parentCharge; charge >= 1; charge--) {
                if (mz >= (parentMass - 20.0) / (double) charge && mz <= (parentMass + 6.0) / (double) charge) {
                    return (true);
                }
            }
        }

    }
    return (false);
}

void SpectraSTPeakList::binPeaksWithScaling(double mzPower, double intensityPower, double unassignedFactor,
                                            int numBinsPerMzUnit, double fractionToNeighbor, bool rebin,
                                            bool removePrecursor, double removeLightIonsMzCutoff) {

    if (m_bins && !rebin) {
        return;
    }
    if (m_bins) {
        delete (m_bins);
    }

    m_numBinsPerMzUnit = numBinsPerMzUnit;

    // calculate the number of bins to divide into
    unsigned int numBins = (MAX_MZ - MIN_MZ + 1) * numBinsPerMzUnit;

    if (m_binIndex) {

        m_bins = new vector<float>;
        m_binIndex->clear();

        int curBinNumber = 0;

        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {

            float inten = scale(*i, mzPower, intensityPower, unassignedFactor, removePrecursor,
                                removeLightIonsMzCutoff);
            if (inten < 0.0001) continue;

            unsigned int binNumber = calcBinNumber(i->mz);
            float neighbor = inten * fractionToNeighbor;
            int shift = binNumber - curBinNumber;

            if (shift < 0) {
                // not sorted. bin it the old-fashioned way
                delete (m_bins);
                m_bins = NULL;
                delete (m_binIndex);
                m_binIndex = NULL;
                binPeaksWithScaling(mzPower, intensityPower, unassignedFactor, numBinsPerMzUnit, fractionToNeighbor,
                                    rebin, removePrecursor, removeLightIonsMzCutoff);
                return;
            }

            int numExistingBins = (int) (m_bins->size());

            if (binNumber > 0) {
                if (shift <= 2 && numExistingBins - 3 + shift >= 0) {
                    (*m_bins)[numExistingBins - 3 + shift] += neighbor;
                } else {
                    m_bins->push_back(neighbor);
                    m_binIndex->push_back(binNumber - 1);
                }
            }

            if (shift <= 1 && numExistingBins - 2 + shift >= 0) {
                (*m_bins)[numExistingBins - 2 + shift] += inten;
            } else {
                m_bins->push_back(inten);
                m_binIndex->push_back(binNumber);
            }

            if (binNumber < numBins - 1) {
                if (shift == 0 && numExistingBins - 1 >= 0) {
                    (*m_bins)[numExistingBins - 1] += neighbor;
                } else {
                    m_bins->push_back(neighbor);
                    m_binIndex->push_back(binNumber + 1);
                }
            }

            curBinNumber = binNumber;
        }

    } else {
        // no binIndex
        m_bins = new vector<float>(numBins, 0.0);

        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {

            float inten = scale(*i, mzPower, intensityPower, unassignedFactor, removePrecursor,
                                removeLightIonsMzCutoff);
            if (inten < 0.0001) continue;

            int binNumber = (int) (calcBinNumber(i->mz));

            // put the full intensity into the bin closest to the m/z of the peak
            (*m_bins)[binNumber] += inten;

            // put a fraction of the intensity into the neigboring bins. This "overflowing"
            // of intensities into the neighbor helps to widen the m/z tolerance when comparing peaks
            float neighbor = inten * fractionToNeighbor;

            if (binNumber > 0) {
                (*m_bins)[binNumber - 1] += neighbor;
            }
            if (binNumber < (int) (numBins - 1)) {
                (*m_bins)[binNumber + 1] += neighbor;
            }

        }

    }

    // calculate the magnitude of the bin vector
    float binsSumOfSquares = 0.0;
    for (vector<float>::iterator j = m_bins->begin(); j != m_bins->end(); j++) {
        binsSumOfSquares += (*j) * (*j);
    }
    float binMagnitude = 0.0;
    if (binsSumOfSquares < 0.0001) {
        m_binMagnitude = 0.01;
    } else {
        m_binMagnitude = sqrt(binsSumOfSquares);
    }

}

// binPeaks - put peaks in bins. (ASSUMING ALREADY SCALED)
void SpectraSTPeakList::binPeaks(int numBinsPerMzUnit, double fractionToNeighbor, bool rebin) {

    if (m_bins && !rebin) {
        return;
    }
    if (m_bins) {
        delete (m_bins);
    }

    m_numBinsPerMzUnit = numBinsPerMzUnit;

    // calculate the number of bins to divide into
    m_numBins = (MAX_MZ - MIN_MZ + 1) * numBinsPerMzUnit;

    if (m_binIndex) {

        m_bins = new vector<float>;
        m_binIndex->clear();

        int curBinNumber = 0;

        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {

            float inten = i->intensity;
            if (inten < 0.0001) continue;

            unsigned int binNumber = calcBinNumber(i->mz);
            float neighbor = inten * fractionToNeighbor;
            int shift = binNumber - curBinNumber;

            if (shift < 0) {
                // not sorted. bin it the old-fashioned way
                delete (m_bins);
                m_bins = NULL;
                delete (m_binIndex);
                m_binIndex = NULL;

                binPeaks(numBinsPerMzUnit, fractionToNeighbor, rebin);
                return;
            }

            int numExistingBins = (int) (m_bins->size());

            if (binNumber > 0) {
                if (shift <= 2 && numExistingBins - 3 + shift >= 0) {
                    (*m_bins)[numExistingBins - 3 + shift] += neighbor;
                } else {
                    m_bins->push_back(neighbor);
                    m_binIndex->push_back(binNumber - 1);
                }
            }

            if (shift <= 1 && numExistingBins - 2 + shift >= 0) {
                (*m_bins)[numExistingBins - 2 + shift] += inten;
            } else {
                m_bins->push_back(inten);
                m_binIndex->push_back(binNumber);
            }

            if (binNumber < m_numBins - 1) {
                if (shift == 0 && numExistingBins - 1 >= 0) {
                    (*m_bins)[numExistingBins - 1] += neighbor;
                } else {
                    m_bins->push_back(neighbor);
                    m_binIndex->push_back(binNumber + 1);
                }
            }

            curBinNumber = binNumber;
        }

    } else {
        // no binIndex
        m_bins = new vector<float>(m_numBins, 0.0);

        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {

            float inten = i->intensity;
            if (inten < 0.0001) continue;

            int binNumber = (int) (calcBinNumber(i->mz));

            // put the full intensity into the bin closest to the m/z of the peak
            (*m_bins)[binNumber] += inten;

            // put a fraction of the intensity into the neigboring bins. This "overflowing"
            // of intensities into the neighbor helps to widen the m/z tolerance when comparing peaks
            float neighbor = inten * fractionToNeighbor;

            if (binNumber > 0) {
                (*m_bins)[binNumber - 1] += neighbor;
            }
            if (binNumber < (int) (m_numBins - 1)) {
                (*m_bins)[binNumber + 1] += neighbor;
            }

        }

    }

    // calculate the magnitude of the bin vector, count the number of occupied (non-zero) bins
    float binsSumOfSquares = 0.0;
    m_numOccupiedBins = 0;
    for (vector<float>::iterator j = m_bins->begin(); j != m_bins->end(); j++) {
        binsSumOfSquares += (*j) * (*j);
        if (*j > 0.01) m_numOccupiedBins++;
    }
    float binMagnitude = 0.0;
    if (binsSumOfSquares < 0.0001) {
        m_binMagnitude = 0.01;
    } else {
        m_binMagnitude = sqrt(binsSumOfSquares);
    }

}


// passFilter - run a filter similar to dtafilter to remove the garbage
// spectra. the options available are similar to that of dtafilter (specified in SpectraSTSearchParams)
bool SpectraSTPeakList::passFilter(SpectraSTSearchParams &searchParams) {

    if (searchParams.filterMaxIntensityBelow > 0.01 && m_origMaxIntensity < searchParams.filterMaxIntensityBelow) {
        return (false);
    }

    vector<Peak>::iterator i;
    unsigned int peakCount = 0;
    float totalIntensity = 0.0;
    float intensityBelowMz = 0.0;
    // vector<float*> peaksAt515;
    // float intensityAt515 = 0.0;
    double minMz = 999999.0;
    double maxMz = 0.0;


    for (i = m_peaks.begin(); i != m_peaks.end(); i++) {

        float intensity = i->intensity;
        if (intensity >= searchParams.filterCountPeakIntensityThreshold &&
            i->mz >= searchParams.filterLightIonsMzThreshold) {
            // peak big enough to be counted
            peakCount++;
            totalIntensity += intensity;
            // if (i->mz < 515.8 && i->mz > 514.8) {
            // remember the 515.3 peak
            //  peaksAt515.push_back(&(i->intensity));
            //  intensityAt515 += intensity;
            // }

            if (i->mz < searchParams.filterAllPeaksBelowMz) {
                intensityBelowMz += intensity;
            }

            if (i->mz > maxMz) maxMz = i->mz;
            if (i->mz < minMz) minMz = i->mz;

        }
    }

    if (maxMz - minMz < searchParams.filterMinMzRange) {
        return (false);
    }

/*
  if (searchParams.filterRemoveHuge515Threshold > 0.001) {
    if (!peaksAt515.empty() && (intensityAt515 > searchParams.filterRemoveHuge515Threshold * totalIntensity)) {
      for (vector<float*>::iterator p515 = peaksAt515.begin(); p515 != peaksAt515.end(); p515++) {
        *(*p515) = 0.0;
      // cout << "Kill a 515.3 peak!" << endl;
      // just zero the peak now, should probably just screen out the spectrum too... who knows
      }
    }
  }
  */

    if (peakCount >= searchParams.filterMinPeakCount && intensityBelowMz < 0.95 * totalIntensity) {
        // that is, with at least searchParams.filterMinPeakCount significant peaks
        // AND not almost all peaks are below m/z = searchParams.filterAllPeaksBelowMz
        return (true);
    } else {
        return (false);
    }

}

// calcBinNumber - returns the bin number for a given m/z
unsigned int SpectraSTPeakList::calcBinNumber(double mz) {

//  if (!m_bins) {
//    return (0);
//  }

    unsigned int mzInt = (unsigned int) (mz);
    if (mzInt < MIN_MZ) {
        return (0);
    }
    if (mzInt > MAX_MZ) {
        return ((MAX_MZ - MIN_MZ + 1) * m_numBinsPerMzUnit - 1);
    }

    return ((unsigned int) (mz * m_numBinsPerMzUnit) - MIN_MZ * m_numBinsPerMzUnit);

}

// calcBinMz - returns the m/z value the correspond to a given bin
double SpectraSTPeakList::calcBinMz(unsigned int binNum) {

//  if (!m_bins) {
//    return (0.0);
//  } else {
    return ((binNum + MIN_MZ * m_numBinsPerMzUnit) / (double) (m_numBinsPerMzUnit));
//  }
}

// calcFractionUnassigned - calculates the fraction of total intensity of the top <maxRank> peaks that is unassigned.
// numAssigned is set to the number of peaks unassigned among the top <maxRank> peaks
double
SpectraSTPeakList::calcFractionUnassigned(unsigned int maxRank, unsigned int &numUnassigned, unsigned int &numAssigned,
                                          bool ignoreNearParentRegion, bool ignoreRareAnnotations) {

    rankByIntensity();

    unsigned int NAA = 0;
    if (m_pep) NAA = m_pep->NAA();

    float totalIntensity = 0.0;
    float unassignedIntensity = 0.0;
    numUnassigned = 0;
    numAssigned = 0;

    for (unsigned int rank = 0; rank < maxRank && rank < m_intensityRanked->size(); rank++) {
        double mz = (*m_intensityRanked)[rank]->mz;
        if (ignoreNearParentRegion && isNearPrecursor(mz)) {
            maxRank++; // skip over this peak, but will allow another to be considered, so as to keep the number of peaks considered constant
            continue;
        }

        float inten = (*m_intensityRanked)[rank]->intensity;
        if (!isAssigned((*m_intensityRanked)[rank]->annotation, NAA, ignoreRareAnnotations)) {
            numUnassigned++;
            unassignedIntensity += inten;
        } else {
            numAssigned++;
        }
        totalIntensity += inten;
    }

    return ((double) (unassignedIntensity / totalIntensity));
}

// isAssigned - given an annotation, determine if that peak is properly assigned to a common ion fragment
// The NAA is given so that fragments of length NAA - 1 are given special treatment (their neutral losses tend to be more common)
bool SpectraSTPeakList::isAssigned(string annotation, unsigned int NAA, bool ignoreRareAnnotations) {

    if (!ignoreRareAnnotations) {
        // consider assigned if annotation does not start with a ?
        return (!annotation.empty() && annotation[0] != '?');
    } else {
        if (annotation.empty() || annotation[0] == '?') return (false);

        string::size_type dummy = 0;
        // take off bracket; chop off piece after slash
        annotation = nextToken(annotation, 0, dummy, " /]\t\r\n", "[");

        // consider all 'a' ions rare
        if (annotation[0] == 'a') return (false);

        if (annotation[0] == 'y' || annotation[0] == 'b') {
            int numAA = atoi(nextToken(annotation.substr(1), 0, dummy, " -^/]\t\r\n").c_str());
            if (numAA == NAA - 1) {
                // no matter what the neutral loss is, a loss from fragment of precursor minus one AA is considered prominent
                return (true);
            }
            string::size_type minusPos = annotation.find('-');
            if (minusPos != string::npos) {
                // is a neutral loss
                string nl = nextToken(annotation, minusPos + 1, dummy, " ^/]\t\r\n");
                if (nl != "17" && nl != "18" && nl != "80" && nl != "98") {
                    return (false);
                } else {
                    return (true);
                }
            }
        }
        return (true);
    }

}

// isBorY - returns if a peak annotation is a common ion
bool SpectraSTPeakList::isBorY(string annotation) {

    if (annotation.empty() || annotation[0] == '?') return (false);

    string::size_type dummy = 0;
    // take off bracket; chop off piece after slash
    annotation = nextToken(annotation, 0, dummy, " /]\t\r\n", "[");

    // consider all 'a' ions rare
    if (annotation[0] != 'y' && annotation[0] != 'b') return (false);

    string::size_type minusPos = annotation.find_first_of("+-");
    if (minusPos != string::npos) {
        return (false);
    }

//  string::size_type iPos = annotation.find('i');
//  if (iPos != string::npos) {
//    return (false);
//  }

    return (true);

}

// normalizeTo - set the top (base) peak to basePeakValue, and scale all others accordingly
void SpectraSTPeakList::normalizeTo(float basePeakValue, float maxDynamicRange) {

    if (m_peaks.size() == 0) return; // empty spectrum

    float maxIntensity = 0.0;
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (maxIntensity < (*i).intensity) {
            maxIntensity = (*i).intensity;
        }
    }

    if (maxIntensity < 0.000001) return; // likely no peak, something's wrong

    float minIntensity = 0.0;
    if (maxDynamicRange > 1.0) {
        minIntensity = basePeakValue / maxDynamicRange;
    }

    vector<Peak> newPeaks;
    m_totalIonCurrent = 0.0;

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        float origIntensity = i->intensity;

        i->intensity *= (basePeakValue / maxIntensity);

        // in order not to have tiny tiny peaks after scaling
        if (i->intensity >= minIntensity) {
            m_totalIonCurrent += origIntensity;
            newPeaks.push_back(*i);
        }
    }

    m_peaks.clear();
    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
    }

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    if (m_intensityRanked) delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

}

// annotate - the spectra given the peptide ion. if redo is true, will annotate again
// even though it has already been annotated. if fixMz is true, it will set the peak's m/z to the first annotation's theoretical
// value
void SpectraSTPeakList::annotate(bool redo, bool fixMz) {

    if (m_isAnnotated && !redo) {
        return;
    }

    if (!m_pep) {
        return;
    }

    // clear all annotation
    for (vector<Peak>::iterator p = m_peaks.begin(); p != m_peaks.end(); p++) {
        (*p).annotation = "";
    }

    vector<FragmentIon *> ions;
    m_pep->generateFragmentIons(ions, m_fragType);

    for (vector<FragmentIon *>::iterator fi = ions.begin(); fi != ions.end(); fi++) {
        Peak *assigned = annotateIon(*fi, fixMz);
        if (assigned) {
            annotateIsotopicIons(assigned, *fi, fixMz);
        }
    }

    if (m_fragType == "CID-QTOF" || m_fragType == "HCD") {
        for (vector<FragmentIon *>::iterator fi = ions.begin(); fi != ions.end(); fi++) {
            if ((*fi)->m_assigned) {
                annotateInternalFragments(*fi, fixMz);
            }
        }
    }


    for (vector<FragmentIon *>::iterator fi = ions.begin(); fi != ions.end(); fi++) {
        delete (*fi);
    }

    // mark as '?' for all unassigned peaks
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if ((*i).annotation.empty()) {
            (*i).annotation = "?";
        } else {
            m_numAssignedPeaks++;
        }
    }

    m_isAnnotated = true;

}


// annotateIon - annotate a particular ion
Peak *SpectraSTPeakList::annotateIon(FragmentIon *fi, bool fixMz) {

    double mz = fi->m_mz;
    unsigned int charge = fi->m_charge;

    double annotateMzAccuracy = 1.0;

    double extraTolerance = 0.0;

    if (m_fragType == "CID-QTOF" || m_fragType == "HCD" || m_fragType == "ETD-HR") {
        annotateMzAccuracy = 0.02;
    } else {
        if (mz * (double) charge > 1400.0) {
            extraTolerance = 1.025 / (double) charge; // approx. mass of neutron + buffer
        }
    }

    createPeakMap();

    // find all peaks within the tolerable range
    map<double, pair<Peak *, double> >::iterator low = m_peakMap->upper_bound(mz - annotateMzAccuracy);
    map<double, pair<Peak *, double> >::iterator high = m_peakMap->lower_bound(
            mz + annotateMzAccuracy + extraTolerance);

    map<double, pair<Peak *, double> >::iterator best = findBestPeakToAssign(mz, low, high,
                                                                             annotateMzAccuracy + extraTolerance);

    if (best == high) {
        return (NULL);
    }

    fi->m_assigned = true;

    Peak *bestPeak = best->second.first;

    if (fixMz && (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[')) {

        bestPeak->mz = mz;
        bestPeak->annotation = fi->m_ion + "/0.000";
        best->second.second = (double) (fi->m_prominence);

    } else {


        double mzDiff = bestPeak->mz - mz;
        double adjProminence = (double) (fi->m_prominence) - fabs(mzDiff);
        stringstream ss0;
        ss0 << fi->m_ion << '/';
        ss0.precision(3);
        ss0 << fixed << mzDiff;

        if (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[') {
            bestPeak->annotation = ss0.str();
            best->second.second = adjProminence;
        } else {
            if (best->second.second < adjProminence) {
                // this annotation is better than the best existing one, put this one in front
                string oldAnnotation = bestPeak->annotation;
                bestPeak->annotation = ss0.str() + "," + oldAnnotation;
                best->second.second = adjProminence;
            } else {
                bestPeak->annotation += "," + ss0.str();
            }
        }

    }

    for (map<double, pair<Peak *, double> >::iterator i = low; i != high; i++) {
        Peak *p = i->second.first;
        if (p->intensity > 1000.0) {
            if (i != best && p->annotation.empty()) {
                double mzDiffBigPeak = p->mz - mz;
                stringstream ssb;
                ssb << '[' << fi->m_ion << '/';
                ssb.precision(3);
                ssb << fixed << mzDiffBigPeak << ']';
                p->annotation = ssb.str();
                i->second.second = 0.0;
            }

        }
    }

    return (bestPeak);
}

// findBestPeakToAssign - finds the best peak within (low, high) to assign to the ion with theoretical m/z = mzTh.
map<double, pair<Peak *, double> >::iterator
SpectraSTPeakList::findBestPeakToAssign(double mzTh, map<double, pair<Peak *, double> >::iterator low,
                                        map<double, pair<Peak *, double> >::iterator high, double mzAccuracy) {

    map<double, pair<Peak *, double> >::iterator i;
    map<double, pair<Peak *, double> >::iterator max = low;
    map<double, pair<Peak *, double> >::iterator best = high;

    if (low == high) {
        // nothing within (low, high)
        return (high);
    }

    // finds the biggest peak within (low, high)
    for (i = low; i != high; i++) {
        if (i->second.first->intensity > max->second.first->intensity) {
            max = i;
        }
    }

    float bestAlignScore = 0.0;

    for (i = low; i != high; i++) {
        double delta = fabs(i->second.first->mz - mzTh) / mzAccuracy;
        // alignScore = sqrt(intensity of this peak / intensity of the biggest peak within (low, high)) - delta^2
        // i.e. it increases with the relative size of the peak and decreases with the distance of it from mzTh
        float alignScore = sqrt(i->second.first->intensity / max->second.first->intensity) - delta * delta;
        if (bestAlignScore < alignScore && i->second.first->intensity >= 1.0) {
            best = i;
            bestAlignScore = alignScore;
        }

    }

    return (best);

}


// annotateIsotopticIons - for each peak assigned, try to find any isotopic ions of that ion and assign them too
// only consider +1 Da (and if +1 Da is found, +2 Da, too). Higher isotopic peaks are not considered.
void SpectraSTPeakList::annotateIsotopicIons(Peak *assignedMonoisotopic, FragmentIon *fi, bool fixMz) {

    double baseMz = assignedMonoisotopic->mz;
    float baseIntensity = assignedMonoisotopic->intensity;
    unsigned int charge = fi->m_charge;

    double isotopicPeakMzTolerance = 0.3;

    if (m_fragType == "CID-QTOF" || m_fragType == "HCD" || m_fragType == "ETD-HR") {
        isotopicPeakMzTolerance = 0.02;
    }

    createPeakMap();
    // first isotopic peak
    map<double, pair<Peak *, double> >::iterator begin = m_peakMap->upper_bound(
            baseMz + (1.00 / (double) charge) - isotopicPeakMzTolerance);
    map<double, pair<Peak *, double> >::iterator end = m_peakMap->lower_bound(
            baseMz + (1.00 / (double) charge) + isotopicPeakMzTolerance);

    if (begin == end) {
        return;
    }

    map<double, pair<Peak *, double> >::iterator max = begin;
    bool found = false;

    //bool heavy = ((baseMz * (double)charge) > 1200.0);
    // true if ion so heavy that the monoisotopic peak is no longer the biggest in the envelope

    bool heavy = true; // don't care if isotopic peak is bigger than the monoisotopic one regardless of m/z

    map<double, pair<Peak *, double> >::iterator i;

    // only assign isotopic peaks to unassigned peaks
    for (i = begin; i != end; i++) {
        Peak *p = i->second.first;
        if ((p->annotation.empty() || p->annotation[0] == '[')
            && (heavy || p->intensity < baseIntensity)
            && p->intensity >= 1.0) {
            if (found) {
                if (p->intensity > max->second.first->intensity) {
                    max = i;
                }
            } else {
                max = i;
            }

            found = true;
        }
    }

    if (!found) {
        return;
    }

    Peak *bestPeak = max->second.first;
    float firstIsotopicPeakIntensity = bestPeak->intensity;

    if (fixMz && (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[')) {
        bestPeak->mz = fi->m_mz + 1.00 / (double) charge;
        // NOTE: due to many different elemental compositions, the distance between the monoisotopic and first higher isotopic peak cannot be more accurate than that
        bestPeak->annotation = fi->m_ion + "i/0.000";
        max->second.second = (double) (fi->m_prominence);

    } else {

        double mzDiff = bestPeak->mz - (fi->m_mz + 1.00 / (double) charge);
        double adjProminence = (double) (fi->m_prominence) - fabs(mzDiff);
        stringstream ss1;
        ss1 << fi->m_ion << 'i' << '/';
        ss1.precision(3);
        ss1 << fixed << mzDiff;

        if (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[') {
            bestPeak->annotation = ss1.str();
            max->second.second = adjProminence;
        } else {
            if (max->second.second < adjProminence) {
                // this annotation is better than the best existing one, put this one in front
                string oldAnnotation = bestPeak->annotation;
                bestPeak->annotation = ss1.str() + "," + oldAnnotation;
                max->second.second = adjProminence;
            } else {
                bestPeak->annotation += "," + ss1.str();
            }
        }
    }

    // second isotopic peak
    begin = m_peakMap->upper_bound(baseMz + (2.00 / (double) charge) - isotopicPeakMzTolerance);
    end = m_peakMap->lower_bound(baseMz + (2.00 / (double) charge) + isotopicPeakMzTolerance);

    if (begin == end) {
        return;
    }

    max = begin;
    found = false;

    for (i = begin; i != end; i++) {
        Peak *p = i->second.first;
        if ((p->annotation.empty() || p->annotation[0] == '[')
            && (heavy || p->intensity < firstIsotopicPeakIntensity)
            && p->intensity >= 1.0) {
            if (found) {
                if (p->intensity > max->second.first->intensity) {
                    max = i;
                }
            } else {
                max = i;
            }

            found = true;
        }
    }

    if (!found) {
        return;
    }

    bestPeak = max->second.first;

    if (fixMz && (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[')) {

        bestPeak->mz = fi->m_mz + 2.00 / (double) charge;
        bestPeak->annotation = fi->m_ion + "i/0.000";
        max->second.second = (double) (fi->m_prominence);

    } else {

        double mzDiff = bestPeak->mz - (fi->m_mz + 2.00 / (double) charge);
        double adjProminence = (double) (fi->m_prominence) - fabs(mzDiff);
        stringstream ss2;
        ss2 << fi->m_ion << 'i' << '/';
        ss2.precision(3);
        ss2 << fixed << mzDiff;

        if (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[') {
            bestPeak->annotation = ss2.str();
            max->second.second = adjProminence;
        } else {
            if (max->second.second < adjProminence) {
                // this annotation is better than the best existing one, put this one in front
                string oldAnnotation = bestPeak->annotation;
                bestPeak->annotation = ss2.str() + "," + oldAnnotation;
                max->second.second = adjProminence;
            } else {
                bestPeak->annotation += "," + ss2.str();
            }
        }
    }


}


void SpectraSTPeakList::annotateInternalFragments(FragmentIon *fi, bool fixMz) {

    char ionType = fi->m_ion[0];

    if ((fi->m_ion[0] != 'y' && fi->m_ion[0] != 'b') || fi->m_loss != 0) {
        return;
    }

    int pepNumAA = m_pep->NAA();
    int fragNumAA = fi->m_pos;

    double baseMz = fi->m_mz;
    unsigned int charge = fi->m_charge;

    double mzTolerance = 0.02;

    createPeakMap();

    double ifMz = baseMz;

    if (ionType == 'b') {
        if (m_pep->isModsSet && !(m_pep->nTermMod.empty())) {
            ifMz -= Peptide::getModMonoisotopicMass(m_pep->nTermMod) / (double) charge;
        }
    }

    if (ionType == 'y') {
        ifMz -= (*Peptide::AAMonoisotopicMassTable)['!'];
        if (m_pep->isModsSet && !(m_pep->cTermMod.empty())) {
            ifMz -= Peptide::getModMonoisotopicMass(m_pep->cTermMod) / (double) charge;
        }
    }

    for (int numLostAA = 1; numLostAA < fragNumAA - 1; numLostAA++) {

        stringstream ionss;
        ionss << 'm';

        if (ionType == 'b') {
            ifMz -= m_pep->monoisotopicMZResidue(numLostAA, charge);
            ionss << numLostAA + 1 << ':' << fragNumAA;
        } else if (ionType == 'y') {
            ifMz -= m_pep->monoisotopicMZResidue(pepNumAA - numLostAA + 1, charge);
            ionss << pepNumAA - fragNumAA + 1 << ':' << pepNumAA - numLostAA;
        }

        map<double, pair<Peak *, double> >::iterator low = m_peakMap->upper_bound(ifMz - mzTolerance);
        map<double, pair<Peak *, double> >::iterator high = m_peakMap->lower_bound(ifMz + mzTolerance);

        map<double, pair<Peak *, double> >::iterator best = findBestPeakToAssign(ifMz, low, high, mzTolerance);

        if (best == high) {
            continue;
        }

        Peak *bestPeak = best->second.first;
        FragmentIon *internalFi = new FragmentIon(ionss.str(), 0, 0, ifMz, charge, 3);

        if (fixMz && (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[')) {

            bestPeak->mz = ifMz;
            bestPeak->annotation = internalFi->m_ion + "/0.000";
            best->second.second = 3.0;

            annotateIsotopicIons(bestPeak, internalFi, fixMz);

        } else {

            double mzDiff = bestPeak->mz - ifMz;
            double adjProminence = 3.0 - fabs(mzDiff);
            stringstream ss0;
            ss0 << internalFi->m_ion << '/';
            ss0.precision(3);
            ss0 << fixed << mzDiff;

            if (bestPeak->annotation.empty() || bestPeak->annotation[0] == '[') {

                bestPeak->annotation = ss0.str();
                best->second.second = adjProminence;

                annotateIsotopicIons(bestPeak, internalFi, fixMz);

            } else if (bestPeak->annotation.find('m') != string::npos) {
                // already have an internal fragment annotation, don't add
            } else {
                if (best->second.second < adjProminence) {
                    // this annotation is better than the best existing one, put this one in front
                    string oldAnnotation = bestPeak->annotation;
                    bestPeak->annotation = ss0.str() + "," + oldAnnotation;
                    best->second.second = adjProminence;
                } else {
                    bestPeak->annotation += "," + ss0.str();
                }
            }
        }

        delete (internalFi);

    }

}


// plot - plot the peak list, for debugging mostly
void SpectraSTPeakList::plot(string name, string label) {

    string specFileName(name + ".spec");
    ofstream specFout;
    if (!myFileOpen(specFout, specFileName)) {
        g_log->error("CREATE", "Cannot open file \"" + specFileName +
                               "\" for writing spectra file for gnuplot. Plotting skipped.");
        return;
    }

    string plotFileName(name + ".gp");
    ofstream plotFout;
    if (!myFileOpen(plotFout, plotFileName)) {
        g_log->error("CREATE",
                     "Cannot open file \"" + specFileName + "\" for writing gnuplot plot file. Plotting skipped.");
        return;
    }

    string pngFileName(name + ".png");

    float maxIntensity = 0.0;
    double minMz = getMinMz();
    double maxMz = getMaxMz();
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (maxIntensity < (*i).intensity) {
            maxIntensity = (*i).intensity;
        }
    }

    minMz = 0.0;

    //maxMz = 1200.0;
    maxMz = 100.0 * ((int) (maxMz / 100) + 1); // round up to the next 100

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (!((*i).annotation.empty()) && i->annotation[0] != '?') {
            specFout << (*i).mz << ' ' << (*i).intensity / maxIntensity << ' ' << 0.0 << ' ' << 0.0 << endl;
        } else {
            // don't label
            specFout << (*i).mz << ' ' << 0.0 << ' ' << (*i).intensity / maxIntensity << ' ' << 0.0 << endl;
        }

    }

    specFout.close();

    plotFout << "set terminal png" << endl;
    plotFout << "set output \"" << pngFileName << "\"" << endl;
    plotFout << "set size " << 1.0 << "," << 1.0 << endl;
    plotFout << "set nokey" << endl;
    plotFout << "set border 1" << endl;
    plotFout << "set xtics border nomirror" << endl;
    plotFout << "set xzeroaxis linetype -1 linewidth 1.0" << endl;
    plotFout << "set yzeroaxis linetype -1 linewidth 1.0" << endl;
    plotFout << "set mxtics" << endl;
    plotFout << "set noytics" << endl;

    plotFout << "set origin 0.0,0.0" << endl;
    plotFout << "set yrange [0.0:1.1]" << endl;

    plotFout.precision(3);
    plotFout << "set label \"" << maxIntensity << "\" at " << minMz << ",1.04 front" << endl;
    plotFout << "set label \"" << label << "\" at " << maxMz << ",1.04 right front" << endl;

    plotFout.precision(3);

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (!((*i).annotation.empty()) && i->annotation[0] != '?' && i->annotation[0] != '[') {
            string::size_type commaPos = 0;
            string labelText("");
            while (commaPos < (*i).annotation.length()) {
                string ion = nextToken((*i).annotation, commaPos, commaPos, ",\t\r\n");
                string::size_type slashPos = ion.find('/', 0);
                ion = ion.substr(0, slashPos);
                if ((ion[0] == 'p' || ion[0] == 'I' || ion.find('-', 0) == string::npos) &&
                    (ion.find('i', 0) == string::npos)) {
                    if (!labelText.empty()) labelText += ",";
                    labelText += ion;
                }
                commaPos++;
            }
            if (!labelText.empty()) {
                plotFout << "set label \"" << labelText << "\" at ";
                plotFout.precision(3);
                plotFout << (*i).mz << ",";
                plotFout.precision(3);
                plotFout << (*i).intensity / maxIntensity + 0.02 << " left rotate front";
                plotFout << " font \"Helvetica,2\" tc lt 1" << endl;
            }
        }

    }

    plotFout.precision(3);
    plotFout << "plot [" << minMz << ":" << maxMz << "] \"" << specFileName << "\" using 1:3 with impulse lc " << 3;
    plotFout << ", \"" << specFileName << "\" using 1:2 with impulse lc " << 1;
    plotFout << ", \"" << specFileName << "\" using 1:4 with impulse lc " << -1 << endl;


    // PLOT!
    myGnuplot(plotFileName);

    // delete the plot file and spec file
//  removeFile(plotFileName);
//  removeFile(specFileName);


}

// createConsensusSpectrum - creates the peak list from multiple peak lists. work in progress
void
SpectraSTPeakList::createConsensusSpectrum(vector<SpectraSTPeakList *> &pls, unsigned int totalNumRep, double quorum,
                                           unsigned int maxNumPeaks, vector<SpectraSTDenoiser *> *denoisers,
                                           bool keepRawIntensities) {

    unsigned int minNumRepWithPeak = (unsigned int) (totalNumRep * quorum - 0.00001) + 1;
    if (minNumRepWithPeak < 1) minNumRepWithPeak = 1;

    double alignMzAccuracy = 1.0;
    if (m_fragType == "CID-QTOF" || m_fragType == "HCD" || m_fragType == "ETD-HR") {
        alignMzAccuracy = 0.1;
    }

    // first, rank all peak lists by intensity
    for (vector<SpectraSTPeakList *>::iterator pl = pls.begin(); pl != pls.end(); pl++) {
        (*pl)->rankByIntensity();
    }

    unsigned int plCount = 0;
    for (vector<SpectraSTPeakList *>::iterator cur = pls.begin(); cur != pls.end(); cur++, plCount++) {
        // for each peak list...

        if (totalNumRep == (unsigned int) (pls.size()) && (unsigned int) (pls.size()) - plCount < minNumRepWithPeak) {
            // not enough spectra to vote the peak in even if all remaining spectra contain that peak
            break;
        }

        vector<Peak *> *pl1 = (*cur)->m_intensityRanked;

        unsigned int pl1max = (unsigned int) (pl1->size()) > maxNumPeaks ? maxNumPeaks : (unsigned int) (pl1->size());

        for (unsigned int r1 = 0; r1 < pl1max; r1++) {
            // r1 is the rank of the peak
            // for each peak in the current peak list...

            if ((*pl1)[r1]->info == "A" || (*pl1)[r1]->info == "S") {
                // if (!((*pl1)[r1])) {
                // already aligned
                continue;
            }

            vector<Peak *> aligned;

            // the percentile is 0 if the peak is the largest, and 1 if the peak is the smallest
            double percentile1 = (double) (r1) / 200.0;
            if (percentile1 > 0.9) percentile1 = 0.9;

            double mz1 = (*pl1)[r1]->mz;
            float intensity1 = (*pl1)[r1]->intensity;

            // parse out the number of replicates used for this peak list
            // this is in case the peak list is a previous consensus. in this case,
            // it has to count more in the voting
            string info1 = (*pl1)[r1]->info;
            unsigned int numRep1 = 1;
            if (!info1.empty()) {
                string::size_type slashPos1 = info1.find('/', 0);
                if (slashPos1 != string::npos) {
                    numRep1 = atoi((info1.substr(0, slashPos1)).c_str());
                }
            }

            double weight1 = (*cur)->m_weight;

            // initializes the m/z and intensity of the new peak
            double wMz1 = weight1 * mz1;
            float wIntensity1 = (float) weight1 * intensity1;
            double sumMz = wMz1;
            double sumMzSq = mz1 * wMz1;
            float sumIntensity = wIntensity1;
            float sumIntensitySq = intensity1 * wIntensity1;
            double sumW = weight1;
            double mzAve = mz1;

            // counting the number of replicates containing this peak
            unsigned int numRepWithPeak = numRep1;

            for (vector<SpectraSTPeakList *>::iterator other = cur + 1; other != pls.end(); other++) {
                // for all subsequent peak lists...

                vector<Peak *> *pl2 = (*other)->m_intensityRanked;

                unsigned int pl2max =
                        (unsigned int) (pl2->size()) > maxNumPeaks ? maxNumPeaks : (unsigned int) (pl2->size());

                for (unsigned int r2 = 0; r2 < pl2max; r2++) {

                    if ((*pl2)[r2]->info == "A" || (*pl2)[r2]->info == "S") {
                        // if (!((*pl2)[r2])) {
                        // already aligned
                        continue;
                    }

                    double percentile2 = (double) (r2) / 200.0;
                    if (percentile2 > 0.9) percentile2 = 0.9;

                    // percentileMax is the "worse" of the two peaks we are trying to align
                    // say if peak 1's percentile is 0/150 (top out of 150 peaks)
                    // and peak 2's percentile is 5/100 (6th peak out of of 100 peaks),
                    // percentileMax gets 5/100.
                    double percentileMax = (percentile1 > percentile2 ? percentile1 : percentile2);

                    double mz2 = (*pl2)[r2]->mz;

                    // the difference in m/z, normalized
                    // gets a value of 1 if the mzDiff is equal to ALIGN_MZ_TOLERANCE
                    double mzDiff = fabs(mzAve - mz2) / alignMzAccuracy;

                    // the align score -- the larger percentileMax is, the smaller the tolerance of alignment is!
                    // e.g. if both peak 1 and peak 2 are the top peaks in their respective spectra, percentileMax = 0,
                    // so they will be aligned if they are within ALIGN_MZ_TOLERANCE.
                    // this trick is to make sure small peaks (possibly noise) won't get too easily aligned
                    double alignScore = (1.0 - percentileMax) - mzDiff;

                    if (alignScore >= 0) {
                        // found a peak in peak list 2 that can be aligned

                        float intensity2 = (*pl2)[r2]->intensity;
                        unsigned int numRep2 = 1;
                        string info2 = (*pl2)[r2]->info;
                        if (!info2.empty()) {
                            string::size_type slashPos2 = info2.find('/', 0);
                            if (slashPos2 != string::npos) {
                                numRep2 = atoi((info2.substr(0, slashPos2)).c_str());
                            }
                        }

                        // weighted average the m/z and intensity together
                        double weight2 = (*other)->m_weight;
                        double wMz2 = weight2 * mz2;
                        float wIntensity2 = (float) weight2 * intensity2;
                        sumMz += wMz2;
                        sumIntensity += wIntensity2;
                        sumMzSq += mz2 * wMz2;
                        sumIntensitySq += intensity2 * wIntensity2;
                        sumW += weight2;
                        mzAve = sumMz / sumW;
                        // mark this peak as aligned already
                        (*pl2)[r2]->info = "A";
                        aligned.push_back((*pl2)[r2]);
                        // (*pl2)[r2] = NULL;

                        // count the votes
                        numRepWithPeak += numRep2;
                        break;
                    }
                }
            }

            // done going through all peak lists searching for this peak

            (*pl1)[r1]->info = "A";
            aligned.push_back((*pl1)[r1]);
            //       (*pl1)[r1] = NULL;

            if (numRepWithPeak >= minNumRepWithPeak) {
                // make the peak quorum, voted in. will include this peak in the final consensus
                Peak newPeak;
                newPeak.mz = sumMz / sumW;
                newPeak.intensity = sumIntensity / sumW;
                newPeak.annotation = "";

                double mzVar = sumMzSq / sumW - newPeak.mz * newPeak.mz;
                double mzStDev = 0.0;
                if (mzVar > 0.000001) {
                    mzStDev = sqrt(mzVar);
                }

                float intensityVar = sumIntensitySq / sumW - newPeak.intensity * newPeak.intensity;
                float intensityStDev = 0.0;
                if (intensityVar > 0.000001) {
                    intensityStDev = sqrt(intensityVar);
                }

                // create the info string
                stringstream ss;
                ss << numRepWithPeak << '/' << minNumRepWithPeak;
                ss.precision(4);
                ss << " " << fixed << mzStDev;
                ss.precision(2);
                ss << "|" << fixed << intensityStDev / newPeak.intensity;

                // scale intensity? let's not do it for now
                // newPeak.intensity *= ((float)numRepWithPeak / (float)totalNumRep);

                newPeak.info = ss.str();
                m_peaks.push_back(newPeak);

                for (vector<Peak *>::iterator ap = aligned.begin(); ap != aligned.end(); ap++) {
                    (*ap)->info = "S"; // enough votes, these are signals in the original replicates
                }
            }
        }

    }

    // sort the peaks by m/z again
    sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
    m_isSortedByMz = true;

    // re-normalize
    if (!keepRawIntensities) {
        normalizeTo(10000.0);
    }

    if (denoisers) {
        // submit the replicates to denoiser for training

        for (vector<SpectraSTPeakList *>::iterator pl = pls.begin(); pl != pls.end(); pl++) {
            unsigned int charge = getParentCharge();
            if (charge > MAX_CHARGE) charge = MAX_CHARGE;

            // add to the all-charge denoiser
            (*denoisers)[0]->addTrainingSpectrum(*pl);

            // add to charge-specific denoiser
            if (charge > 0) (*denoisers)[charge]->addTrainingSpectrum(*pl);

        }
    }

    annotate();

}

// rankByIntensity - sort peaks by intensity and place in m_intensityRanked object
void SpectraSTPeakList::rankByIntensity(bool redo, unsigned int maxRank, bool removePrecursor,
                                        double removeLightIonsMzCutoff) {

    if (m_intensityRanked && !redo) {
        return;
    }

    if (m_intensityRanked) {
        delete (m_intensityRanked);
    }

    m_intensityRanked = new vector<Peak *>;

    for (vector<Peak>::iterator p = m_peaks.begin(); p != m_peaks.end(); p++) {
        if ((removePrecursor && isNearPrecursor(p->mz)) || (p->mz < removeLightIonsMzCutoff)) {
            p->intensity = 0.0;
            continue;
        }
        m_intensityRanked->push_back(&(*p));
    }

    if (m_intensityRanked->empty()) return;

    if (maxRank < m_intensityRanked->size() - 1) {
        partial_sort(m_intensityRanked->begin(), m_intensityRanked->begin() + maxRank, m_intensityRanked->end(),
                     SpectraSTPeakList::sortPeakPtrsByIntensityDesc);
    } else {
        sort(m_intensityRanked->begin(), m_intensityRanked->end(), SpectraSTPeakList::sortPeakPtrsByIntensityDesc);
    }
}

void SpectraSTPeakList::createPeakMap(bool redo) {

    if (m_peakMap && !redo) {
        return;
    }

    if (m_peakMap) delete (m_peakMap);

    m_peakMap = new map<double, pair<Peak *, double> >;

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        pair<Peak *, double> p;
        p.first = &(*i);
        p.second = 0.0;
        (*m_peakMap)[(*i).mz] = p;
    }

}

// calcSignalToNoise - calculates a "signal-to-noise" for this peak list
double SpectraSTPeakList::calcSignalToNoise() {

    if (m_signalToNoise > 0.00001) {
        return (m_signalToNoise);
    }

    if (m_peaks.size() <= 10) return (1.0);  // too few peaks - garbage spectrum

    rankByIntensity();

    unsigned int numPeaks = (unsigned int) (m_intensityRanked->size());
    unsigned int numSignalPeaks = 4;

    if (m_intensityRanked->size() < 40) {
        numSignalPeaks = numPeaks / 10 + 1;
    }

    // consider the average of (3rd, 4th, 5th, 6th) as signal
    float sumSignal = 0.0;
    unsigned int numUsedSignal = 0;
    unsigned int numSignal = 0;
    for (unsigned int r = 0; r < numPeaks; r++) {
        Peak *p = (*m_intensityRanked)[r];
        if (!isNearPrecursor(p->mz)) {
            numSignal++;
            if (numSignal > 2) {
                sumSignal += p->intensity;
                numUsedSignal++;
            }
        }
        if (numUsedSignal >= numSignalPeaks) break;
    }

    if (numUsedSignal == 0) return (1.0); // can't find a signal, not good

    float signal = sumSignal / (float) (numUsedSignal);

    unsigned int noiseIndex = numPeaks - 1;

    if (numPeaks > 40) {
        noiseIndex = 40 + (numPeaks - 40) / 2;   // median after first 40 peaks
    }

    // problem with using median as noise level is that spectra with many small (noise) peaks will
    // get an exaggerated S/N value, as its median is much deeper into the noise regime
    // so, for spectra with >200 peaks, cap the noise level at the 200th peak
    if (noiseIndex > 200) noiseIndex = 200;

    float noise = (*m_intensityRanked)[noiseIndex]->intensity;

    m_signalToNoise = (double) (signal / noise);

    if (m_signalToNoise >= 400.0) {
        m_signalToNoise = 400.0;
    }

    return m_signalToNoise;
}

// calcXrea - calculates the "Xrea" of a peak list (REF), a measure of spectral quality (NOT used in SpectraST)
double SpectraSTPeakList::calcXreaOld() {

    if (m_peaks.size() < 6) return (0.0);

    rankByIntensity();

    float slope = m_totalIonCurrent / (float) m_intensityRanked->size();

    float cumInten = 0.0;
    float diagonal = 0.0;
    float xrea = 0.0;
    float triangle = 0.0;

    for (int rank = (int) m_intensityRanked->size() - 1; rank >= 0; rank--) {
        diagonal += slope;
        cumInten += (*m_intensityRanked)[rank]->intensity;
        xrea += diagonal - cumInten;
        triangle += diagonal;
    }

    xrea = xrea / triangle;

    return ((double) (xrea));
}

// calcXrea - calculates the "Xrea" of a peak list (REF), a measure of spectral quality
double SpectraSTPeakList::calcXrea(bool ignorePrecursorRegion) {

    if (m_peaks.size() < 10) return (0.0);

    rankByIntensity();

    // float slope = m_totalIonCurrent / (float) m_intensityRanked->size();

    float cumInten = 0.0;
    float cumCumInten = 0.0;
    int numPeaksUsed = 0;
    float intensity = 0.0;
//  int maxRank = 300; // cap at 300 peaks
//  if (maxRank > (int)(m_intensityRanked->size())) maxRank = (int)(m_intensityRanked->size());

    for (int rank = (int) (m_intensityRanked->size()) - 1; rank >= 0; rank--) {
        if (ignorePrecursorRegion && isNearPrecursor((*m_intensityRanked)[rank]->mz)) {

        } else {
            numPeaksUsed++;
            intensity = (*m_intensityRanked)[rank]->intensity;
            cumInten += intensity;
            cumCumInten += cumInten;
        }
    }

    if (numPeaksUsed < 10) return (0.0);

    float triangle = (cumInten * (float) (numPeaksUsed + 1)) * 0.5;
    float xrea = (triangle - cumCumInten) / (triangle /* + intensity */); // alpha = intensity of highest peak

    return ((double) (xrea));
}

// calcXCorr - calculates the SEQUEST-like cross correlation (not quite the same as SEQUEST)
// the preprocessing in SEQUEST is slightly different
double SpectraSTPeakList::calcXCorr() {

    if (!m_pep) return (0.0);

    m_parentMz = m_pep->monoisotopicMZ();

    // do not use binIndex
    if (m_binIndex) delete (m_binIndex);

    binPeaksWithScaling(0.0, 1.0, 1.0, 1, 0.5, true, true, 0.0); // force rebinning

    unsigned int numBins = (unsigned int) (m_bins->size());
    unsigned int regionSize = numBins / 10 + 1;

    float maxIntInRegion[10];
    bool emptyRegion[10];
    for (unsigned int r = 0; r < 10; r++) {
        maxIntInRegion[r] = 0.1;
        emptyRegion[r] = true;
    }

    for (unsigned int b = 0; b < numBins; b++) {
        unsigned int region = b / regionSize;
        if (maxIntInRegion[region] < (*m_bins)[b]) {
            maxIntInRegion[region] = (*m_bins)[b];
            emptyRegion[region] = false;
        }
    }

    for (unsigned int b = 0; b < numBins; b++) {
        unsigned int region = b / regionSize;
        if (!emptyRegion[region]) {
            (*m_bins)[b] = (*m_bins)[b] * 50.0 / maxIntInRegion[region];
        } else {
            (*m_bins)[b] = 0.0;
        }
    }

    map<int, float> theoSpec;
    m_pep->SEQUESTTheoreticalSpectrum(theoSpec);
    vector<float> theoBins;
    theoBins.assign(numBins, 0.0);

    float theoNormFac = 0.0;
    for (map<int, float>::iterator t = theoSpec.begin(); t != theoSpec.end(); t++) {
        theoBins[calcBinNumber((double) (t->first))] += t->second;
    }

    float sumDots = 0.0;
    float xcorr = 0.0;
    for (int tau = -75; tau <= 75; tau++) {

        float dot = 0.0;
        for (unsigned int b = 0; b < numBins; b++) {
            if (b + tau < (int) numBins) {
                dot += (*m_bins)[b] * theoBins[b + tau] / 10000.0;
            }
        }

        sumDots += dot;

        if (tau == 0) {
            xcorr = dot;
        }
    }
    xcorr = xcorr - (sumDots / 151.0);

    // for unknown reasons, xcorr calculated this way is smaller than SEQUEST xcorr.
    // empirically, SEQUEST xcorr is about 1.45 times bigger
    xcorr *= 1.45;

    // don't reuse the bins -- they're not binned in the regular way
    delete (m_bins);
    m_bins = NULL;

    return ((double) (xcorr));
}


// setWeight
void SpectraSTPeakList::setWeight(double weight) {
    m_weight = weight;
}

// removeInquoratePeaks - go through the peak list and remove all peaks without enough replicates
void SpectraSTPeakList::removeInquoratePeaks(unsigned int minNumRepWithPeak) {


    bool changed = false;
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (!(i->info.empty())) {

            string::size_type spacePos = 0;
            string repStr = nextToken(i->info, spacePos, spacePos, " \t\r\n");

            string::size_type slashPos = repStr.find('/');

            if (slashPos == string::npos || slashPos >= repStr.length() - 1) {
                continue;
            }

            unsigned int origMinNumRepWithPeak = atoi(repStr.substr(slashPos + 1).c_str());
            if (minNumRepWithPeak <= origMinNumRepWithPeak) {
                // peak quorum the same as before, or even looser, so no need to remove
                // we're assuming the rest of the peaks will have same quorum - so just quit altogether
                return;
            }
            unsigned int N = atoi(repStr.substr(0, slashPos).c_str());

            if (N < minNumRepWithPeak) {
                // inquorate -- remove this peak. for efficiency,
                // instead of deleting the peak, set the intensity to zero
                // will rewrite the entire peak list without the zero peaks later
                i->intensity = 0.0;
                changed = true;
            }
        }
    }

    if (changed) {

        // rewriting
        vector<Peak> newPeaks;
        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
            if (i->intensity > 0.001) {
                if (!i->info.empty()) {
                    string::size_type spacePos = 0;
                    string repStr = nextToken(i->info, spacePos, spacePos, " \t\r\n");
                    string rest = i->info.substr(spacePos);

                    string::size_type slashPos = repStr.find('/');
                    stringstream newInfoss;
                    newInfoss << repStr.substr(0, slashPos) << '/' << minNumRepWithPeak << rest;
                    i->info = newInfoss.str();
                }
                newPeaks.push_back(*i);
            }
        }

        m_peaks.clear();
        for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
            m_peaks.push_back(*j);
        }

        m_isScaled = false;

        if (m_bins) delete (m_bins);
        m_bins = NULL;

        if (m_intensityRanked) delete (m_intensityRanked);
        m_intensityRanked = NULL;

        if (m_peakMap) delete (m_peakMap);
        m_peakMap = NULL;

    }


}

// removeITRAQPeaks - go through the peak list and remove the ITRAQ quantitation peaks at 113 -121 Th
void SpectraSTPeakList::removeITRAQReporterPeaks(bool setToZeroOnly) {

    bool changed = false;
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (i->mz > 112.0 && i->mz < 122.0) {
            i->intensity = 0.0;
            changed = true;
        }
    }

    if (setToZeroOnly) {
        return;
    }

    if (changed) {

        // rewriting
        vector<Peak> newPeaks;
        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
            if (i->intensity > 0.001) {
                newPeaks.push_back(*i);
            }
        }

        m_peaks.clear();
        for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
            m_peaks.push_back(*j);
        }

        m_isScaled = false;

        if (m_bins) delete (m_bins);
        m_bins = NULL;

        if (m_intensityRanked) delete (m_intensityRanked);
        m_intensityRanked = NULL;

        if (m_peakMap) delete (m_peakMap);
        m_peakMap = NULL;

    }
}

// removeTMTReporterPeaks - go through the peak list and remove the TMT quantitation peaks at 126 -132 Th
void SpectraSTPeakList::removeTMTReporterPeaks(bool setToZeroOnly) {

    bool changed = false;
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (i->mz > 126.0 && i->mz < 132.0) {
            i->intensity = 0.0;
            changed = true;
        }
    }

    if (setToZeroOnly) {
        return;
    }

    if (changed) {

        // rewriting
        vector<Peak> newPeaks;
        for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
            if (i->intensity > 0.001) {
                newPeaks.push_back(*i);
            }
        }

        m_peaks.clear();
        for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
            m_peaks.push_back(*j);
        }

        m_isScaled = false;

        if (m_bins) delete (m_bins);
        m_bins = NULL;

        if (m_intensityRanked) delete (m_intensityRanked);
        m_intensityRanked = NULL;

        if (m_peakMap) delete (m_peakMap);
        m_peakMap = NULL;

    }
}

void SpectraSTPeakList::quickSimplify(unsigned int maxNumPeaks, double maxDynamicRange,
                                      bool removePrecursor, double removeLightIonsMzCutoff) {

    if ((unsigned int) (m_peaks.size()) <= maxNumPeaks && m_origMaxIntensity <= maxDynamicRange && !removePrecursor &&
        removeLightIonsMzCutoff < 0.001)
        return;

    rankByIntensity(true, maxNumPeaks + 1, removePrecursor, removeLightIonsMzCutoff);

    if (m_intensityRanked->empty()) return;

    double basePeakIntensity = (*m_intensityRanked)[0]->intensity;
    double minIntensity = basePeakIntensity / maxDynamicRange;

    for (unsigned int r = 1; r <= m_intensityRanked->size(); r++) {
        Peak *p = (*m_intensityRanked)[r - 1];
        if (r > maxNumPeaks || p->intensity < minIntensity) {
            p->intensity = 0.0; // instead of deleting the peak, it just sets it to zero. this is faster
        }
    }

    if (m_bins) delete (m_bins);
    m_bins = NULL;

}

// simplify - trims the peak list such that there are at most maxNumPeaks and the dynamic range (highest intensity / lowest intensity)
// is at most maxDynamicRange. If deisotope is true, peaks that are likely higher isotopic peaks are removed. If excludeParent is true,
// peaks near the precursorMz are removed. Returns the fraction of total intensity that is retained after simplifying.
double SpectraSTPeakList::simplify(unsigned int maxNumPeaks, double maxDynamicRange,
                                   bool deisotope, bool removePrecursor, double removeLightIonsMzCutoff) {

    if ((unsigned int) (m_peaks.size()) <= maxNumPeaks
        && m_origMaxIntensity <= maxDynamicRange
        && !deisotope && !removePrecursor && removeLightIonsMzCutoff < 0.001)
        return (1.0);

    vector<Peak> newPeaks;

    unsigned int numPeaks;

    if (!deisotope) {

        rankByIntensity(true, maxNumPeaks + 1, removePrecursor, removeLightIonsMzCutoff);

        if (m_intensityRanked->empty()) return (1.0);

        double basePeakIntensity = (*m_intensityRanked)[0]->intensity;
        double minIntensity = basePeakIntensity / maxDynamicRange;
        numPeaks = (unsigned int) (m_intensityRanked->size());

        for (unsigned int r = 0; r < maxNumPeaks && r < numPeaks; r++) {
            Peak p = *((*m_intensityRanked)[r]);

            if (p.intensity < minIntensity) break;
            newPeaks.push_back(p);
        }

    } else {

        rankByIntensity(true, maxNumPeaks * 3, removePrecursor,
                        removeLightIonsMzCutoff); // sort more peaks in case they are isotopes

        if (m_intensityRanked->empty()) return (1.0);

        double basePeakIntensity = (*m_intensityRanked)[0]->intensity;
        double minIntensity = basePeakIntensity / maxDynamicRange;
        numPeaks = (unsigned int) (m_intensityRanked->size());

        for (unsigned int r = 0; r < maxNumPeaks && r < numPeaks; r++) {
            Peak p = *((*m_intensityRanked)[r]);

            if (p.intensity < minIntensity) break;

            if (p.annotation.empty()) {
                bool isIsotope = false;
                for (vector<Peak>::iterator k = newPeaks.begin(); k != newPeaks.end(); k++) {
                    double diff = p.mz - k->mz;
                    if (diff > 0.0 && diff < 2.3) {
                        isIsotope = true;
                        break;
                    }
                }
                if (isIsotope) {
                    maxNumPeaks++;
                } else {
                    newPeaks.push_back(p);
                }
            } else {
                string::size_type pos = 0;
                if (nextToken(p.annotation, 0, pos, ",/\t\r\n").find('i') != string::npos) {
                    maxNumPeaks++;
                } else {
                    newPeaks.push_back(p);
                }
            }
        }
    }
    // sort the new peaks by m/z
    sort(newPeaks.begin(), newPeaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);

    double totalScaledIntensity = 0.0;
    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if ((removePrecursor && isNearPrecursor(i->mz)) || (i->mz < removeLightIonsMzCutoff)) {
            continue;
        }
        string::size_type pos = 0;
        if (deisotope && nextToken(i->annotation, 0, pos, ",/\t\r\n").find('i') != string::npos) {
            continue;
        }
        totalScaledIntensity += sqrt(i->intensity);
    }

    double retainedScaledIntensity = 0.0;
    m_peaks.clear();
    m_totalIonCurrent = 0.0;
    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
        m_totalIonCurrent += j->intensity;
        retainedScaledIntensity += sqrt(j->intensity);
    }

    // peaks are inserted in m/z order
    m_isSortedByMz = true;

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

    return (retainedScaledIntensity / totalScaledIntensity);
}

// clean - performs cleaning of the spectrum. It does three things -- each can be turned on or off:
// (i) deisotope, (ii) remove near-precursor ions, and (iii) remove light ions. Similar to simplify() but
// does not remove peaks based on intensities.
double SpectraSTPeakList::clean(bool deisotope, bool removePrecursor, double removeLightIonsMzCutoff) {

    if (!deisotope && !removePrecursor && removeLightIonsMzCutoff < 0.001) return (1.0);

    vector<Peak> newPeaks;

    rankByIntensity(true, 999999, removePrecursor, removeLightIonsMzCutoff);
    unsigned int numPeaks = (unsigned int) (m_intensityRanked->size());

    if (m_intensityRanked->empty()) return (1.0);

    for (unsigned int r = 0; r < numPeaks; r++) {

        Peak p = *((*m_intensityRanked)[r]);

        if (p.annotation.empty()) {

            // try to find if there is any larger peak than this peak that is within -2.3 to -0.0 Th of this peak
            // quite inefficient!
            for (vector<Peak>::iterator k = newPeaks.begin(); k != newPeaks.end(); k++) {
                double diff = p.mz - k->mz;
                if (diff > 0.0 && diff < 2.3) {
                    newPeaks.push_back(p);
                    break;
                }
            }

        } else {
            string::size_type pos = 0;
            if (nextToken(p.annotation, 0, pos, ",/\t\r\n").find('i') == string::npos) {
                newPeaks.push_back(p);
            }
        }
    }
    // sort the new peaks by m/z
    sort(newPeaks.begin(), newPeaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);

    double totalIntensity = m_totalIonCurrent;
    double retainedIntensity = 0.0;

    m_peaks.clear();
    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
        retainedIntensity += j->intensity;
    }

    // peaks are inserted in m/z order
    m_isSortedByMz = true;

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

    m_totalIonCurrent = retainedIntensity;
    return (retainedIntensity / totalIntensity);
}


// getFracUnassignedStr - returns a string containing fraction unassigned information
string SpectraSTPeakList::getFracUnassignedStr() {

    float fracUnassigned = 0.0;
    float fracUnassignedTop20 = 0.0;
    float fracUnassignedTop5 = 0.0;
    unsigned int numUnassignedPeaksTop20 = 0;
    unsigned int numUnassignedPeaksTop5 = 0;
    unsigned int numUnassignedPeaks = 0;

    rankByIntensity();
    annotate();

    if (!m_pep) return ("");

    unsigned int numAA = m_pep->NAA();

    unsigned int top5 = 5;
    unsigned int top20 = 20;
    bool isNearParent = false;
    bool assigned = false;
    bool top20done = false;

    float totalIntensity = 0.0;
    float totalIntensityTop20 = 0.0;
    float totalIntensityTop5 = 0.0;

    for (unsigned int r = 0; r < (unsigned int) (m_intensityRanked->size()); r++) {
        Peak *p = (*m_intensityRanked)[r];

        totalIntensity += p->intensity;
        if (p->annotation.empty() || p->annotation == "?") {
            numUnassignedPeaks++;
            fracUnassigned += p->intensity;
        }

        if (!top20done) {
            isNearParent = isNearPrecursor(p->mz);
            assigned = isAssigned(p->annotation, numAA, true);
        }

        if (r < top5) {
            if (isNearParent) {
                // don't count peak near parent in the top 5.
                top5++;
            } else {
                totalIntensityTop5 += p->intensity;
                if (!assigned) {
                    numUnassignedPeaksTop5++;
                    fracUnassignedTop5 += p->intensity;
                }
            }
        }

        if (r < top20) {
            if (isNearParent) {
                // don't count peak near parent in the top 20
                top20++;
            } else {
                totalIntensityTop20 += p->intensity;
                if (!assigned) {
                    numUnassignedPeaksTop20++;
                    fracUnassignedTop20 += p->intensity;
                }
            }
        } else {
            top20done = true;
        }
    }

    if (totalIntensity > 0.00001) fracUnassigned /= totalIntensity;
    if (totalIntensityTop5 > 0.00001) fracUnassignedTop5 /= totalIntensityTop5;
    if (totalIntensityTop20 > 0.00001) fracUnassignedTop20 /= totalIntensityTop20;

    stringstream fracUnassignedss;
    fracUnassignedss.precision(2);
    fracUnassignedss << fixed << showpoint << fracUnassignedTop5 << ',' << numUnassignedPeaksTop5 << "/5;";
    fracUnassignedss << fixed << showpoint << fracUnassignedTop20 << ',' << numUnassignedPeaksTop20 << "/20;";
    fracUnassignedss << fixed << showpoint << fracUnassigned << ',' << numUnassignedPeaks << '/' << m_peaks.size();

    return (fracUnassignedss.str());

}

// getNthLargestPeak - gets the Nth largest peak of the peak list
void SpectraSTPeakList::getNthLargestPeak(unsigned int n, Peak &p) {

    rankByIntensity(false, 99999, false, 0.0);
    // note: this is the most basic sort, no removing peaks.
    // But since "redo" is set to false, if the caller wants other behavior, he/she can call rankByIntensity with his/her desired arguments before this.

    if (n < 1 || n > (unsigned int) (m_intensityRanked->size())) {
        p.mz = 0.0;
        p.intensity = 0.0;
        p.annotation = "";
        p.info = "";
        return;
    }

    p = *((*m_intensityRanked)[n - 1]);

}

// repositionPeaks - reposition the peaks based on the peptide identification (m_pep)
void SpectraSTPeakList::repositionPeaks(bool keepEffectivePeakCountConstant) {

    //  rankByIntensity();

    int isotope = 0;
    // unsigned int numMovedMinorPeaks = 0;

    // for (unsigned int r = 0; r < (unsigned int)(m_intensityRanked->size()); r++) {

    //   Peak* p = (*m_intensityRanked)[r];
    for (vector<Peak>::iterator p = m_peaks.begin(); p != m_peaks.end(); p++) {
        string annotation = p->annotation;

        if (p->mz < 400.0) continue;

        if (keepEffectivePeakCountConstant && isNearPrecursor(p->mz)) {
            continue;
        }

        /*
    if (annotation.empty() || annotation[0] == '?') {
      double shift = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
      p->mz += shift;
      continue;
    }
    */

        if (annotation.empty() || annotation[0] == '?') {
            // don't move
            continue;
        }

        if (annotation[0] == '[') {
            annotation = annotation.substr(1);
        }

        if (annotation[0] == 'I') {
            // TODO: don't know how to deal with this
            continue;
        }

        char ionType = annotation[0];

        string::size_type pos = 0;
        unsigned int numAA = atoi(nextToken(annotation, 1, pos, "^i-+:/]*, \t\r\n", " \t\r\n").c_str());

        unsigned int startAAPos = 0;
        unsigned int endAAPos = 0;
        if (ionType == 'm' && annotation[pos] == ':') {
            startAAPos = numAA;
            endAAPos = atoi(nextToken(annotation, pos + 1, pos, "^i-+/]*, \t\r\n", " \t\r\n").c_str());
        }

        int neutralLoss = 0;
        int charge = 1;
        double mzDiff = 0.0;

        if (annotation.find('i') == string::npos) {
            isotope = 0;
        } else {
            isotope++;
        }

        if (pos < annotation.length() && (annotation[pos] == '-' || annotation[pos] == '+')) {
            // neutral loss
            neutralLoss = atoi(nextToken(annotation, pos + 1, pos, "^i*/], \t\r\n", " \t\r\n").c_str());

            if (annotation[pos] == '+') neutralLoss *= -1;
        }

        bool isNISTisotope = false;
        if (pos < annotation.length() &&
            (annotation[pos] == 'i' || annotation[pos] == '*')) { // NIST format, i before charge
            pos++;
        }

        if (pos < annotation.length() && (annotation[pos] == '^')) {
            charge = atoi(nextToken(annotation, pos + 1, pos, "i*/], \t\r\n", " \t\r\n").c_str());
        }

        if (pos < annotation.length() && (annotation[pos] == 'i')) { // SpectraST format, i after charge
            pos++;
        }

        if (pos < annotation.length() && annotation[pos] == '/') {
            mzDiff = atof(nextToken(annotation, pos + 1, pos, " ,\t\r\n", " \t\r\n").c_str());
        }

        double newMz = 0.0;
        if (ionType == 'p') {
            newMz = m_pep->monoisotopicMH() / (double) charge;
        } else if (ionType == 'm') {
            newMz = m_pep->monoisotopicMZInternalFragment(startAAPos, endAAPos, charge);
        } else {
            newMz = m_pep->monoisotopicMZFragment(ionType, numAA, charge);
        }

        newMz -= (neutralLoss / (double) charge);

        newMz += (isotope / (double) charge);
        newMz += mzDiff;

        if (newMz < 400.0) {
            continue;
        }

        if (keepEffectivePeakCountConstant && newMz < 100.0) {
            // don't move this peak if it will be below 100 Th
            continue;
        }

        // in case the peptide changes in mass
        m_parentMz = m_pep->monoisotopicMZ();
        m_parentCharge = m_pep->charge;

        if (keepEffectivePeakCountConstant && isNearPrecursor(newMz)) {
            continue;
        }

        //if (numMovedMinorPeaks >= 10 && (neutralLoss > 0 || ionType == 'a' || isotope > 1)) {
        //  continue;
        // }

        string firstAnnotation = nextToken(p->annotation, 0, pos, " ,\t\r\n", " \t\r\n");

        p->mz = newMz;
        p->annotation = firstAnnotation; // + "  MOVED";

        //if (neutralLoss > 0 || ionType == 'a' || isotope > 1) {
        //  numMovedMinorPeaks++;
        //}


    }

    sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
    m_isSortedByMz = true;
}

// sortPeakByIntensity - comparison function used by sort() to sort peaks
bool SpectraSTPeakList::sortPeaksByIntensityDesc(Peak a, Peak b) {

    return (a.intensity > b.intensity);

}

// sortPeakByIntensity - comparison function used by sort() to sort peaks
bool SpectraSTPeakList::sortPeakPtrsByIntensityDesc(Peak *a, Peak *b) {

    return (a->intensity > b->intensity);

}

// sortPeakByMzAsc - comparison function used by sort() to sort peaks
bool SpectraSTPeakList::sortPeaksByMzAsc(Peak a, Peak b) {

    return (a.mz < b.mz);
}


// sortByMScoreDesc - comparison function for sorting "MScore" (suitability as MRM transition)
bool SpectraSTPeakList::sortByMScoreDesc(pair<int, unsigned int> a, pair<int, unsigned int> b) {

    if (a.first > b.first) {
        return (true);
    } else if (a.first == b.first) {
        return (a.second < b.second);
    } else {
        return (false);
    }

}

// generateTheoreticalSpectrum - experimenting
void SpectraSTPeakList::generateTheoreticalSpectrum() {

    if (!m_pep) return;

    vector<FragmentIon *> ions;
    m_pep->generateFragmentIons(ions);

    for (vector<FragmentIon *>::iterator fi = ions.begin(); fi != ions.end(); fi++) {

        if ((*fi)->m_prominence >= 7) {
            float intensity = 10000.0;
            if ((*fi)->m_prominence == 8) intensity = 5000.0;
            if ((*fi)->m_prominence == 7) intensity = 500.0;

            insert((*fi)->m_mz, intensity, (*fi)->getAnnotation(), "");
        }

        delete (*fi);
    }

    sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
    m_isSortedByMz = true;
}

// reduce - reduces the spectra to at most maxNumPeaks, using MRM transition selection rules. work in progress...
double SpectraSTPeakList::reduce(unsigned int maxNumPeaks, double minMz, double maxMz, unsigned int numRepsUsed) {

    if (m_peaks.size() <= maxNumPeaks) {
        return (1.0);
    }

    annotate();

    rankByIntensity();

    vector<pair<int, unsigned int> > mScores(m_intensityRanked->size());

    unsigned int rankWithoutBadPeaks = 0;

    double maxIntensity = 10000.0;
    if (m_intensityRanked->size() > 0) maxIntensity = (*m_intensityRanked)[0]->intensity;

    for (unsigned int r = 0; r < m_intensityRanked->size(); r++, rankWithoutBadPeaks++) {
        mScores[r].first = 0;
        mScores[r].second = r;
        Peak *p = (*m_intensityRanked)[r];

        string::size_type pos = 0;
        string ion = nextToken(p->annotation, 0, pos, "/, \t\r\n");

        // outside mass range -- completely useless
        if (p->mz < minMz || p->mz > maxMz) {
            rankWithoutBadPeaks--;
            continue;
        }

        // requires annotation that's not an isotopic (or slight mass shifted ion)
        if (ion[0] != '?' && ion[0] != '[' && ion.find('i') == string::npos) {
            mScores[r].first += 100000000;
        } else {
            // don't count this one in ranking
            rankWithoutBadPeaks--;
        }

        // annotation does not contain precursor and non-backbone losses,
        // and fragment m/z not within 5 Th of the precursor m/z (wider?)
        if (p->annotation.find_first_of("pI") == string::npos && fabs(p->mz - m_parentMz) > 5.0) {
            mScores[r].first += 10000000;
        }

        // not due to a cleavage less than 3 AA's from either terminus
        int aa = atoi(nextToken(ion, 1, pos, "^i-/, \t\r\n").c_str());
        if (aa >= 3) { // && (m_pep && m_pep->NAA() - aa >= 3)) {
            mScores[r].first += 1000000;
        }

        // top quartile in intensity or at least 20% of base intensity
        if (rankWithoutBadPeaks < m_intensityRanked->size() / 4 || (p->intensity / maxIntensity > 0.2)) {
            mScores[r].first += 100000;
        }

        // not a neutral loss
        if (ion.find_first_of("-") == string::npos) {
            mScores[r].first += 10000;
        }

        // y ion
        if (p->annotation[0] == 'y') {
            mScores[r].first += 1000;
        }

        // ion heavier than precursor m/z (by at least 20 Th)
        if (p->mz > m_parentMz + 20.0) {
            mScores[r].first += 100;
        }

        // reproducility considerations
        if (!p->info.empty()) {
            unsigned int nreps = atoi(nextToken(p->info, 0, pos, "/ \t\r\n").c_str());
            // int quorum = atoi(nextToken(p->info, pos + 1, pos, " \t\r\n").c_str());
            unsigned int quorum = (int) ((double) numRepsUsed * 0.9 +
                                         0.5); // need 1/1, 2/2, 3/3, 4/4, 5/5, 5/6, 6/7, 7/8 etc.

            // more replicates than bare minimum -- mainly useful for small Nreps spectra
            if (nreps > quorum) {
                mScores[r].first += 10;
            }

            string mzVarStr = nextToken(p->info, pos + 1, pos, "/ \t\r\n");
            string intVarStr = nextToken(p->info, pos + 1, pos, " \t\r\n");

            if (!intVarStr.empty()) {
                // intensity does not vary among replicates more than 20% of value
                double intVar = atof(intVarStr.c_str());
                if (intVar < 0.2) {
                    mScores[r].first += 1;
                }
            } else {
                // NIST format
                double intVar = atof(mzVarStr.c_str());
                if (intVar < 5.0) {
                    mScores[r].first += 1;
                }
            }

        }

        // tags on the MScore at the end of each peak's info
        stringstream msss;
        msss << p->info << ' ' << mScores[r].first;
        p->info = msss.str();
    }


    sort(mScores.begin(), mScores.end(), SpectraSTPeakList::sortByMScoreDesc);

    vector<Peak> newPeaks;
    for (unsigned int m = 0; m < maxNumPeaks; m++) {
        Peak newp = *((*m_intensityRanked)[mScores[m].second]);

        string::size_type pos = 0;
        string firstAnnotation = nextToken(newp.annotation, 0, pos, ",\t\r\n");
        string::size_type slashPos = firstAnnotation.find('/');
        if (slashPos != string::npos) {
            double mzShift = atof(firstAnnotation.substr(slashPos + 1).c_str());
            newp.mz -= mzShift;
        }

        newPeaks.push_back(newp);
    }

    double retainedIntensity = 0.0;
    m_peaks.clear();
    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
        retainedIntensity += j->intensity;
    }

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

    return (retainedIntensity / m_totalIonCurrent);
}

// centroid - centroiding method...
void SpectraSTPeakList::centroid(string instrument) {

    // presumed resolution - should be conservative?
    double res400 = 10000.0; // for TOF
    if (instrument == "FT") {
        res400 = 100000.0;
    } else if (instrument == "Orbitrap") {
        res400 = 50000.0;
    }

    if (!m_isSortedByMz) {
        sort(m_peaks.begin(), m_peaks.end(), sortPeaksByMzAsc);
        m_isSortedByMz = true;
    }

    // determine smallest m/z-interval between peaks
    // this should tell us the frequency with which the
    // mass spectrometer takes readings
    //
    // this step is helpful because in many profile spectra,
    // the peak is omitted completely if the intensity is below
    // a certain threshold. So the peak list has jumps in m/z values,
    // and neighboring peaks are not necessarily close in m/z.

    double minInterval = 1000000.0;
    double minIntensity = 1000000.0;

    for (int p = 0; p < (int) (m_peaks.size()) - 1; p++) {
        double interval = m_peaks[p + 1].mz - m_peaks[p].mz;
        if (interval < 0.0) {
            cerr << "Peak list not sorted by m/z. No centroiding done." << endl;
            return;
        }
        if (minInterval > interval) minInterval = interval;
        if (minIntensity > m_peaks[p].intensity) minIntensity = m_peaks[p].intensity;
    }

    vector<Peak> *allPeaks = new vector<Peak>;
    for (int i = 0; i < (int) m_peaks.size() - 1; i++) {
        allPeaks->push_back(m_peaks[i]);
        double gap = m_peaks[i + 1].mz - m_peaks[i].mz;
        double curMz = m_peaks[i].mz;
        int numZeros = 0;
        while (gap > 1.9 * minInterval) {
            if (numZeros < 3 || curMz > m_peaks[i + 1].mz - 3.1 * minInterval) {
                curMz += minInterval;
            } else {
                curMz = m_peaks[i + 1].mz - 3.0 * minInterval;
                gap = 4.0 * minInterval;
            }
            Peak pk;
            pk.mz = curMz;
            pk.intensity = 0.0;
            allPeaks->push_back(pk);
            gap -= minInterval;
            numZeros++;
        }

    }

    //  m_peaks.clear();
    // for (vector<Peak>::iterator al = allPeaks.begin(); al != allPeaks.end(); al++) {
    //   m_peaks.push_back(*al);
    // }
    // plot(m_pep->getSafeName() + "_filled", "");

    vector<Peak> *smoothedPeaks = new vector<Peak>;
    for (int i = 0; i < (int) allPeaks->size(); i++) {

        int weight = 6;
        Peak smoothPeak;
        smoothPeak.mz = (*allPeaks)[i].mz;
        smoothPeak.intensity = 6 * (*allPeaks)[i].intensity;

        if (i >= 2) {
            weight += 1;
            smoothPeak.intensity += (*allPeaks)[i - 2].intensity;
        }
        if (i >= 1) {
            weight += 4;
            smoothPeak.intensity += 4 * (*allPeaks)[i - 1].intensity;
        }
        if (i < (int) m_peaks.size() - 1) {
            weight += 4;
            smoothPeak.intensity += 4 * (*allPeaks)[i + 1].intensity;
        }
        if (i < (int) m_peaks.size() - 2) {
            weight += 1;
            smoothPeak.intensity += (*allPeaks)[i + 2].intensity;
        }

        smoothPeak.intensity /= (float) weight;

        smoothedPeaks->push_back(smoothPeak);

    }

    delete (allPeaks);

    //  m_peaks.clear();
    // for (vector<Peak>::iterator sm = smoothedPeaks.begin(); sm != smoothedPeaks.end(); sm++) {
    //   m_peaks.push_back(*sm);
    // }
    // plot(m_pep->getSafeName() + "_smoothed", "");

    int j;
    float maxIntensity;
    int bestPeak;
    bool bLastPos;

    int nextBest;
    double FWHM;

    bLastPos = false;

    vector<Peak> newPeaks;
    //step along each point in spectrum
    for (int i = 0; i < (int) smoothedPeaks->size() - 1; i++) {

        //check for change in direction
        if ((*smoothedPeaks)[i].intensity < (*smoothedPeaks)[i + 1].intensity) {

            bLastPos = true;
            continue;

        } else {

            if (bLastPos) {
                bLastPos = false;

                //find max
                //This is an artifact of using a window of length n (user variable) for identifying
                //a better peak apex on low-res data. Feel free to remove this section if desired.
                //Replace with:
                //  bestPeak=j;
                //  maxIntensity=s[j].intensity;

                maxIntensity = 0;
                for (j = i; j < i + 1; j++) {
                    if ((*smoothedPeaks)[j].intensity > maxIntensity) {
                        maxIntensity = (*smoothedPeaks)[j].intensity;
                        bestPeak = j;
                    }
                }

                //Best estimate of Gaussian centroid
                //Get 2nd highest point of peak
                if (bestPeak == smoothedPeaks->size() - 1) {
                    nextBest = bestPeak - 1;
                } else if ((*smoothedPeaks)[bestPeak - 1].intensity > (*smoothedPeaks)[bestPeak + 1].intensity) {
                    nextBest = bestPeak - 1;
                } else {
                    nextBest = bestPeak + 1;
                }

                //Get FWHM of Orbitrap
                //This is the function you must change for each type of instrument.
                //For example, the FT would be:
                //  FWHM = s[bestPeak].mz * s[bestPeak].mz / (400*res400);
                // Orbitrap
                // FWHM = (*smoothedPeaks)[bestPeak].mz * sqrt((*smoothedPeaks)[bestPeak].mz) / (20 * res400);

                if (instrument == "FT") {
                    FWHM = (*smoothedPeaks)[bestPeak].mz * (*smoothedPeaks)[bestPeak].mz / (400 * res400);
                } else if (instrument == "Orbitrap") {
                    FWHM = (*smoothedPeaks)[bestPeak].mz * sqrt((*smoothedPeaks)[bestPeak].mz) / (20 * res400);
                } else {
                    // for TOF
                    FWHM = (*smoothedPeaks)[bestPeak].mz / res400;
                }

                Peak centroid;
                //Calc centroid MZ (in three lines for easy reading)
                // centroid.mz = pow(FWHM , 2) * log((*smoothedPeaks)[bestPeak].intensity / (*smoothedPeaks)[nextBest].intensity);
                // centroid.mz /= 8 * log(2.0) * ((*smoothedPeaks)[bestPeak].mz - (*smoothedPeaks)[nextBest].mz);
                // centroid.mz += ((*smoothedPeaks)[bestPeak].mz + (*smoothedPeaks)[nextBest].mz) / 2;
                centroid.mz = (*smoothedPeaks)[bestPeak].mz;

                //Calc centroid intensity
                // double exponent = pow(((*smoothedPeaks)[bestPeak].mz - centroid.mz) / FWHM, 2) * (4 * log(2.0));
                //
                // if (exponent > 1.0) {
                //  centroid.intensity = (*smoothedPeaks)[bestPeak].intensity;
                // } else {
                //  centroid.intensity = (*smoothedPeaks)[bestPeak].intensity / exp(-exponent);
                // }

                centroid.intensity = (*smoothedPeaks)[bestPeak].intensity;

                //Hack until I put in mass ranges
                //Another fail-safe can be made for inappropriate intensities
                if (centroid.mz < 0 || centroid.mz > 2000 || centroid.intensity < 0.99 * minIntensity) {
                    //do nothing if invalid mz
                } else {
                    newPeaks.push_back(centroid);
                }

            }

        }
    }

    delete (smoothedPeaks);

    m_peaks.clear();
    m_origMaxIntensity = 0.0;
    m_totalIonCurrent = 0.0;
    m_isAnnotated = false;

    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
        if (m_origMaxIntensity < j->intensity) m_origMaxIntensity = j->intensity;
        m_totalIonCurrent += j->intensity;
    }

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;
    //  plot(m_pep->getSafeName() + "_centroided", "");
}


// hasConsecutiveIonSeries - determines if the peak list contains consecutive ion series (work in progress)
bool SpectraSTPeakList::hasConsecutiveIonSeries() {

    if (!m_pep) return (false);

    annotate();

    rankByIntensity();

    vector<int> bSeries(m_pep->NAA(), 999);
    vector<int> ySeries(m_pep->NAA(), 999);

    for (int r = 0; r < 150 && r < (int) m_intensityRanked->size(); r++) { // top 150 peaks only

        string annotation = ((*m_intensityRanked)[r])->annotation;
        string::size_type commaPos = 0;
        while (commaPos < annotation.length()) {
            string ion = nextToken(annotation, commaPos, commaPos, ",\t\r\n");
            if (ion[0] != 'b' && ion[0] != 'y') {
                commaPos++;
                continue;
            }
            string::size_type iPos = ion.find('i', 0);
            if (iPos != string::npos) {
                commaPos++;
                continue;
            }
            int charge = 1;
            int pos = 0;
            string::size_type hatPos = ion.find('^', 0);
            string::size_type slashPos = ion.find('/', 0);

            string::size_type minusPos = ion.find('-', 0);
            if (minusPos != string::npos && minusPos < slashPos) {
                commaPos++;
                continue;
            }

            if (hatPos == string::npos) {
                charge = 1;
                pos = atoi(ion.substr(1, slashPos - 1).c_str());
            } else {
                charge = atoi((ion.substr(hatPos + 1, slashPos - hatPos - 1)).c_str());
                pos = atoi(ion.substr(1, hatPos - 1).c_str());
            }

            if (pos < 1 || pos > (int) m_pep->NAA() - 1) {
                commaPos++;
                continue;
            }

            if (ion[0] == 'b') {
                if (charge < bSeries[pos]) bSeries[pos] = charge;
            } else if (ion[0] == 'y') {
                if (charge < ySeries[pos]) ySeries[pos] = charge;
            }
            commaPos++;
        }

    }

    unsigned int bScore = 0;
    unsigned int yScore = 0;

    cout << ' ';
    for (int i = 1; i < (int) m_pep->NAA(); i++) {
        cout << (bSeries[i] < 10 ? bSeries[i] : 0);
        if (bSeries[i] < 10 && bSeries[i - 1] < 10 && bSeries[i - 1] <= bSeries[i]) {
            bScore++;

        }
    }
    cout << ' ';
    for (int i = 1; i < (int) m_pep->NAA(); i++) {
        cout << (ySeries[i] < 10 ? ySeries[i] : 0);

        if (ySeries[i] < 10 && ySeries[i - 1] < 10 && ySeries[i - 1] <= ySeries[i]) {
            yScore++;

        }
    }

    cout << ' ';
    double minScore = (double) (m_pep->NAA()) / 3.0 + 1.0;
    if (minScore > 7.0) minScore = 7.0;

    cout << "    b=" << bScore << " y=" << yScore << " min=" << minScore << " combMin=" << 1.4 * minScore;
    if ((double) bScore >= minScore || (double) yScore >= minScore || (double) (bScore + yScore) >= 1.6 * minScore) {
        cout << " GOOD" << endl;
        return (true);
    } else {
        cout << " BAD" << endl;
        return (false);
    }

}

void SpectraSTPeakList::shiftAllPeaks(double mzShift, double randomizeRange) {

    vector<Peak> newPeaks;

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {

        double newMzShift = mzShift;

        if (randomizeRange > 0.001) {
            srand(time(NULL));
            newMzShift += ((double) rand() / (double) RAND_MAX * 2.0 * randomizeRange - randomizeRange);
        }

        Peak p;
        p.mz = i->mz + newMzShift;
        p.intensity = i->intensity;
        p.annotation = i->annotation;
        p.info = i->info;
        newPeaks.push_back(p);
    }

    m_peaks.clear();

    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
    }

    sort(m_peaks.begin(), m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
    m_isSortedByMz = true;

    m_isScaled = false;

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    if (m_intensityRanked) delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

}

void SpectraSTPeakList::flattenAllPeaks() {

    for (vector<Peak>::iterator j = m_peaks.begin(); j != m_peaks.end(); j++) {
        j->intensity = 10000.0;
    }

}

void SpectraSTPeakList::removeNoncanonicalPeaks() {

    for (vector<Peak>::iterator j = m_peaks.begin(); j != m_peaks.end(); j++) {
        if (!(isBorY(j->annotation))) {
            j->intensity = 0.0;
        }
    }
}

Peak *SpectraSTPeakList::findPeakPtr(double mz, double tolerance, map<double, double> *foundMZs) {

    createPeakMap();

    map<double, pair<Peak *, double> >::iterator low = m_peakMap->upper_bound(mz - tolerance);
    map<double, pair<Peak *, double> >::iterator high = m_peakMap->lower_bound(mz + tolerance);

    if (low == m_peakMap->end() || low == high) {
        // can't find any peak within the tolerance of m/z
        return NULL;
    }

    map<double, pair<Peak *, double> >::iterator i;
    map<double, pair<Peak *, double> >::iterator max = low;

    // finds the biggest peak within (low, high)
    for (i = low; i != high; i++) {
        if (foundMZs->find(i->second.first->mz) == foundMZs->end() &&
            i->second.first->intensity > max->second.first->intensity) {
            max = i;
        }
    }

    return (max->second.first);

}

double SpectraSTPeakList::findPeak(double mz, double tolerance) {

    createPeakMap();

    map<double, pair<Peak *, double> >::iterator low = m_peakMap->upper_bound(mz - tolerance);
    map<double, pair<Peak *, double> >::iterator high = m_peakMap->lower_bound(mz + tolerance);

    if (low == m_peakMap->end() || low == high) {
        // can't find any peak within the tolerance of m/z
        return (0.0);
    }

    map<double, pair<Peak *, double> >::iterator i;
    map<double, pair<Peak *, double> >::iterator max = low;

    // finds the biggest peak within (low, high)
    for (i = low; i != high; i++) {
        if (i->second.first->intensity > max->second.first->intensity) {
            max = i;
        }
    }

    return (max->second.first->intensity);

}

bool SpectraSTPeakList::isSinglyCharged() {

    if (m_parentCharge > 1) return (false);

    float totalIntensity = 0.0;
    float totalIntensityAboveParent = 0.0;
    unsigned int numPeaksAboveParent = 0;

    for (vector<Peak>::iterator i = m_peaks.begin(); i != m_peaks.end(); i++) {
        if (!isNearPrecursor(i->mz)) totalIntensity += i->intensity;
        if (i->mz > m_parentMz + 20.0) {
            totalIntensityAboveParent += i->intensity;
            numPeaksAboveParent++;
        }
    }

    if (totalIntensityAboveParent < 0.2 * totalIntensity || numPeaksAboveParent <= 3) {
        return (true);
    }

    return (false);
}

bool SpectraSTPeakList::passFilterUnidentified(SpectraSTCreateParams &params) {

    if (m_parentMz < 350.0) {
        return (false);
    }

    if (getMzRange() < 350.0) {
        return (false);
    }

    if (getNumPeaks() < (unsigned) params.unidentifiedMinimumNumPeaksToInclude) {
        return (false);
    }

    if (params.unidentifiedRemoveSinglyCharged && isSinglyCharged()) {
        return (false);
    }

    return (true);

}

void SpectraSTPeakList::prepareForSearch(SpectraSTSearchParams &params, bool isLibrarySpectrum) {

    if (m_isSearchOnly) return;

    if (params.filterITRAQReporterPeaks) {
        removeITRAQReporterPeaks(true);
    }

    if (params.filterTMTReporterPeaks) {
        removeTMTReporterPeaks(true);
    }

    if (!(params.useSp4Scoring)) {

        scalePeaks(params.peakScalingMzPower, params.peakScalingIntensityPower, params.peakScalingUnassignedPeaks,
                   false, true, params.filterLightIonsMzThreshold);

        if (params.useRankTransformWithQuota) {

            rankTransformWithQuota(params.filterLibMaxPeaksUsed, true, params.filterLightIonsMzThreshold,
                                   params.useRankTransformWithQuotaWindowSize,
                                   params.useRankTransformWithQuotaNumberOfPeaks);

        } else {

            // use the same number of peaks for library and query spectra
            rankTransform(params.filterLibMaxPeaksUsed, true, params.filterLightIonsMzThreshold);
        }

    } else {

        scalePeaks(params.peakScalingMzPower, params.peakScalingIntensityPower, params.peakScalingUnassignedPeaks,
                   false, true, params.filterLightIonsMzThreshold);

        if (isLibrarySpectrum) {
            quickSimplify(params.filterLibMaxPeaksUsed, 999999, true, params.filterLightIonsMzThreshold);
        } else {
            quickSimplify(params.filterMaxPeaksUsed, params.filterMaxDynamicRange, true,
                          params.filterLightIonsMzThreshold);
        }
    }

    if (params.peakNoBinning) {
        prepareForNoBinningDot();

    } else {
        binPeaks(params.peakBinningNumBinsPerMzUnit, params.peakBinningFractionToNeighbor, false);
        m_peaks.clear(); // this saves memory -- all dot product calculations only need the bins
    }


    if (m_intensityRanked) delete (m_intensityRanked);
    m_intensityRanked = NULL;

    if (m_peakMap) delete (m_peakMap);
    m_peakMap = NULL;

    m_isSearchOnly = true;

}

void SpectraSTPeakList::prepareForNoBinningDot() {

    vector<Peak> newPeaks;

    float sumSqInt = 0.0;

    for (vector<Peak>::iterator p = m_peaks.begin(); p != m_peaks.end(); p++) {
        if (p->intensity < 0.01) continue;

        newPeaks.push_back(*p);
        sumSqInt += (p->intensity * p->intensity);
    }

    m_peaks.clear();

    for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
        m_peaks.push_back(*j);
    }

    m_peakMagnitude = 0.0;
    if (sumSqInt > 0.01) {
        m_peakMagnitude = sqrt(sumSqInt);
    }

    if (m_bins) delete (m_bins);
    m_bins = NULL;

    // NOTE: This function does not change the order of peaks in m_peaks (it merely deletes some of them)
    // So m_isSortedByMz is not touched.

}


void SpectraSTPeakList::useBinIndex() {

    if (m_binIndex) return;

    m_binIndex = new vector<unsigned int>;
    if (m_bins) delete (m_bins);

}


void SpectraSTPeakList::generatePoissonCutoffTable(int maxNumTrials, double alpha) {

    // Given k, solve equation ln(alpha * k!) = kln(L) -L for L using Newton's method

    if (poissonCutoffTable) delete (poissonCutoffTable);

    poissonCutoffTable = new vector<double>(maxNumTrials + 1, 1e20);

    double lhs = log(alpha);
    bool alwaysSignificant = false;

    for (int k = 1; k <= maxNumTrials; k++) {

        double lnk = log((double) k);
        lhs += lnk;

        // cerr << "c=" << lhs << ";k=" << k << ";maxc=" << k * log(k) - k;

        if (alwaysSignificant || k * lnk - k < lhs) {

            // cerr << "| NO WAY" << endl;
            // no chance p(k) will reach alpha for any k greater than or equal to this
            (*poissonCutoffTable)[k] = 1e20; // set to some huge number that will always be greater than the real lambda
            alwaysSignificant = true;
            continue;
        }

        double L = exp(lhs / (double) k); // first guess
        double change = 1.0; // convergence criteria, change <= 0.001

        while (change > 0.001) {

            double newL = L - (k * log(L) - L - lhs) / (k / L - 1.0);
            // cerr << "|" << newL;

            change = fabs(newL - L);
            L = newL;
        }

        // cerr << '|' << L << endl;
        (*poissonCutoffTable)[k] = L;
    }

}

int SpectraSTPeakList::findMinSignificantNumMatchedPeaks(int numBins, int numLibraryPeaks, int numQueryPeaks) {

    if (!SpectraSTPeakList::poissonCutoffTable) return (0);

    double L = (double) numLibraryPeaks * (double) numQueryPeaks / (double) numBins;

    for (int k = 1; k <= SpectraSTPeakList::poissonCutoffTable->size(); k++) {
        if ((*poissonCutoffTable)[k] > L) return (k);
    }
    return (SpectraSTPeakList::poissonCutoffTable->size() + 1);

}

