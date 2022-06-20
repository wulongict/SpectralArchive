#include "SpectraSTSearch.hpp"
#include "SpectraSTLibEntry.hpp"
#include "FileUtils.hpp"

#include <algorithm>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>


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

/* Class: SpectraSTSearch
 * 
 * Implements the search for one query spectrum.
 *  
 */


using namespace std;


extern bool g_verbose;

// constructor
SpectraSTSearch::SpectraSTSearch(SpectraSTQuery *query, SpectraSTSearchParams &params, SpectraSTSearchOutput *output) :
        m_query(query),
        m_params(params),
        m_output(output),
        m_candidates() {

    m_candidates.clear();

}

// destructor
SpectraSTSearch::~SpectraSTSearch() {

    if (m_query) delete m_query;

    for (vector<SpectraSTCandidate *>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
        delete (*i);
    }
}

// search - main function to perform one search
void SpectraSTSearch::search(SpectraSTLib *lib) {

    double precursorMz = m_query->getPrecursorMz();

    m_query->getPeakList()->prepareForSearch(m_params, false);

    if (g_verbose) {
        cout << endl;
        cout << "Now searching query: " << m_query->getName() << " (PrecursorMZ = " << precursorMz;

        if (m_query->getDefaultCharge() != 0) {
            cout << "; PrecursorCharge = " << m_query->getDefaultCharge();
            for (int ch = 1; ch <= 7; ch++) {
                if (m_query->isPossibleCharge(ch) && ch != m_query->getDefaultCharge()) cout << "," << ch;
            }
        }

        cout << ")" << endl;
    }

    // retrieves all entries from the library within the tolerable m/z range
    vector<SpectraSTLibEntry *> entries;

    // get more than needed to catch isotope error or (for low resolution data) observed average mass
    double lowMz = precursorMz - m_params.precursorMzTolerance - 2.0;
    double highMz = precursorMz + m_params.precursorMzTolerance;
    bool shortAnnotation = true;

    if (m_params.useTierwiseOpenModSearch) {
        // need full annotation to locate modification
        shortAnnotation = false;
    }

    lib->retrieve(entries, lowMz, highMz, shortAnnotation);

    // for all retrieved entries, do the necessary filtering, add the good ones to m_candidates
    for (vector<SpectraSTLibEntry *>::iterator i = entries.begin(); i != entries.end(); i++) {

        //   if (!(m_params.ignoreChargeOneLibSpectra && (*i)->getCharge() == 1) &&
        //       !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "ICAT_cl" && (*i)->isCleavableICAT()) &&
        //       !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "ICAT_uc" && (*i)->isUncleavableICAT()) &&
        //       !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "CAM" && (*i)->isCAMCysteine()) &&
        //       !(m_params.ignoreAbnormalSpectra && (*i)->getStatus() != "Normal") &&
        //       !(m_params.ignoreSpectraWithUnmodCysteine && (*i)->hasUnmodifiedCysteine()) &&

        if (isWithinPrecursorTolerance(*i)) {

            (*i)->getPeakList()->prepareForSearch(m_params, true);

            SpectraSTCandidate *newCandidate = new SpectraSTCandidate(*i, m_params);
            m_candidates.push_back(newCandidate);


        }
    }

    if (g_verbose) {
        cout << "\tFound " << m_candidates.size() << " candidate(s)... " << " Comparing... ";
        cout.flush();
    }

    unsigned int numHits = 0;

    // compare query to each candidate by calculating the dot product
    for (vector<SpectraSTCandidate *>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {

        SpectraSTLibEntry *entry = (*i)->getEntry();
        int charge = entry->getCharge();

        double dot = 0.0;
        double dotBias = 0.0;
        int numTiersUsed = 0;

        if (m_params.useSp4Scoring) {
            dot = m_query->getPeakList(charge)->calcDotAndDotBias(entry->getPeakList(), dotBias);
        } else {
            if (m_params.useTierwiseOpenModSearch) {
                dot = m_query->getPeakList(charge)->calcDotTierwiseOpenModSearch(entry->getPeakList(), 0.5 /
                                                                                                       (double) (m_params.peakBinningNumBinsPerMzUnit),
                                                                                 (*i)->getSimScoresRef().openMod,
                                                                                 numTiersUsed);
                dotBias = (double) numTiersUsed;
            } else if (m_params.peakNoBinning) {
                dot = m_query->getPeakList(charge)->calcDotNoBinning(entry->getPeakList(), 0.5 /
                                                                                           (double) (m_params.peakBinningNumBinsPerMzUnit));
                dotBias = (double) (m_query->getPeakList()->getNumPeaks()) / 100.0;
            } else {
                dot = m_query->getPeakList(charge)->calcDot(entry->getPeakList());
            }

        }

        (*i)->getSimScoresRef().dot = dot;
        (*i)->getSimScoresRef().dotBias = dotBias;
        (*i)->setSortKey(dot);

        double openModMz = 0.0;
        //  double openModMz = (*i)->getSimScoresRef().openMod.first / (double)charge;

        if (m_params.precursorMzUseAverage) {
            ((*i)->getSimScoresRef()).precursorMzDiff = precursorMz - entry->getAveragePrecursorMz() - openModMz;
        } else {
            ((*i)->getSimScoresRef()).precursorMzDiff = precursorMz - entry->getPrecursorMz() - openModMz;
        }

        if (dot > 0.01) numHits++;
    }

    // sort the hits by the sort key
    // (in this case, the value of "dot" returned by the SpectraSTPeakList::compare function)
    sort(m_candidates.begin(), m_candidates.end(), SpectraSTCandidate::sortPtrsDesc);

    if (m_params.detectHomologs > 1) {
        detectHomologs();
    } else {
        if (m_candidates.size() > 1) {
            (m_candidates[0]->getSimScoresRef()).firstNonHomolog = 2;
        }
    }

    if (m_params.usePValue) {
        finalizeScoreUsePValue();
    } else {
        finalizeScoreLinearCombination(numHits);
    }

    if (m_params.useSp4Scoring) {

        // sort again by F value (due to dot bias, sorting by dot and by F value could be different
        sort(m_candidates.begin(), m_candidates.end(), SpectraSTCandidate::sortPtrsDesc);
    }

}

bool SpectraSTSearch::isWithinPrecursorTolerance(SpectraSTLibEntry *entry) {

    double mzDiff = m_query->getPrecursorMz() - entry->getPrecursorMz();
    double mzTol = m_params.precursorMzTolerance;

    if (mzTol < 0.5) { // high mass accuracy data

        int charge = entry->getCharge();
        if (charge == 0) charge = 1;

        double mzDiffFirstIsotope = mzDiff - 1.0 / charge;
        //  double mzDiffSecondIsotope = mzDiff - 2.0 / charge;

        if ((fabs(mzDiff) <= mzTol || fabs(mzDiffFirstIsotope) <= mzTol) && // || fabs(mzDiffSecondIsotope) <= mzTol) &&
            (m_params.searchAllCharges || m_query->isPossibleCharge(entry->getCharge()))) {
            return (true);
        } else {
            return (false);
        }

    } else { // low mass accuracy data

        if (fabs(mzDiff) <= mzTol &&
            (m_params.searchAllCharges || m_query->isPossibleCharge(entry->getCharge()))) {
            return (true);
        } else {
            return (false);
        }
    }

}

// finalizeScoreSP4 - calculates the deltas, fval and hit stats for SP5
void SpectraSTSearch::finalizeScoreLinearCombination(unsigned int numHits) {

    if (m_candidates.empty()) {
        return;
    }

    // unsigned int numHits = (unsigned int)(m_candidates.size());
    double totalDot = 0;
    double totalSqDot = 0;

    for (unsigned int rank = 0; rank < (unsigned int) (m_candidates.size()); rank++) {
        SpectraSTSimScores &s = m_candidates[rank]->getSimScoresRef();
        totalDot += s.dot;
        totalSqDot += s.dot * s.dot;

        s.hitsNum = numHits;

        if (s.firstNonHomolog > 0) {
            // delta is dot - dot(first non-homologous hit lower than this one)
            s.delta = s.dot - (m_candidates[s.firstNonHomolog - 1]->getSimScoresRef()).dot;
        } else {
            s.delta = 0.0;
        }

        double fval = s.calcFval(m_params.fvalFractionDelta, m_params.fvalUseDotBias);

        m_candidates[rank]->setSortKey(fval);
    }

    double mean = totalDot / (double) numHits;
    double stdev = totalSqDot / (double) numHits - mean * mean;
    if (stdev > 0.000001) stdev = sqrt(stdev);

    for (vector<SpectraSTCandidate *>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
        SpectraSTSimScores &s = (*i)->getSimScoresRef();
        s.hitsMean = mean;
        s.hitsStDev = stdev;
    }

}

// finalizeScoreSP5 - calculates the deltas, fval and hit stats for SP5
void SpectraSTSearch::finalizeScoreUsePValue() {

    if (m_candidates.empty()) {
        return;
    }

    SpectraSTSimScores &top = m_candidates[0]->getSimScoresRef();
    if (top.dot < 0.2) {
        return;
    }

    double gaussianMean = 0.40;
    double gaussianStDev = 0.08;

    int numHits = (int) (m_candidates.size());

    // int startRank = top.firstNonHomolog - 1;
    // if (startRank < 0) startRank = 0;
    int startRank = 5;
    int endRank = numHits - 5;

    int sampleSize = 0;

    double totalDot = 0.0;
    double totalSqDot = 0.0;
    double mean = 0.0;
    double variance = 0.0;
    double stdev = 0.0;
    double KSScore = 0.0;
    double ADScore = 0.0;

    for (int rank = startRank; rank < endRank; rank++) {

        SpectraSTSimScores &s = m_candidates[rank]->getSimScoresRef();

        // if (s.dot < 0.001) s.dot = 0.001;
        if (s.dot < 0.04) break;

        double sqrtDot = sqrt(s.dot);  // sqrt of dot for Gaussian fit

        totalDot += sqrtDot;
        totalSqDot += (sqrtDot * sqrtDot);
        sampleSize++;

    }

    if (sampleSize > 0) {

        mean = totalDot / (double) sampleSize;
        double meanSq = totalSqDot / (double) sampleSize;
        variance = meanSq - mean * mean;
        if (variance > 0.000001) stdev = sqrt(variance);

    } else {

        mean = 0.40; // default

    }

    bool useDefault = false;

    if (sampleSize < 50 || stdev < 0.05) {

        useDefault = true;

        // not enough to fit, just use default above
        // cerr << "DEFAULT n= " << sampleSize << " Mean= " << mean << " StDev= " << stdev << endl;

    }

    if (useDefault) {

        gaussianMean = mean; // 0.40;
        gaussianStDev = 0.08;

    } else {

        gaussianMean = mean;
        gaussianStDev = stdev;

        // cerr << "FITTED n= " << sampleSize << " Mean= " << mean << " StDev= " << stdev << endl;

    }

    // done estimating distribution, now calculate Fval and record all mean/stdev in lower hits

    // also estimate KS and AD scores for the fitted data points
    KSScore = 0.0;
    ADScore = 0.0;
    vector<double> allPValues;

    for (int rank = 0; rank < numHits; rank++) {

        SpectraSTSimScores &s = m_candidates[rank]->getSimScoresRef();

        if (s.firstNonHomolog > 0) {
            // delta is dot - dot(first non-homologous hit lower than this one)
            s.delta = s.dot - (m_candidates[s.firstNonHomolog - 1]->getSimScoresRef()).dot;
        } else {
            s.delta = 0.0;
        }

        // Gumbel
        //if (s.pValue < 0) s.pValue = 1 - exp(-exp(-(s.dot - gumbelMu) / gumbelBeta));

        // Gamma
        // s.pValue = gsl_sf_gamma_inc_Q(gammaAlpha, s.dot / gammaBeta);

        // Weibull
        // if (s.pValue < 0) s.pValue = exp(-pow(s.dot / weibullLambda, weibullK));

        double sqrtDot = sqrt(s.dot);
        long double z = (long double) (sqrtDot - gaussianMean) / (long double) gaussianStDev;

        s.pValue = erfc(z * 0.70710678) * 0.5000000000;

        if (rank >= startRank && sqrtDot > 0.1 && sampleSize > 0) {
            // ignore zero dots when calculating K-S score or A-D score
            double empiricalPValue = (double) (rank - startRank + 1) / (double) sampleSize;
            double diffPValue = fabs(s.pValue - empiricalPValue);

            if (diffPValue > KSScore) KSScore = diffPValue;
            // cerr << "rank = " << rank << ", startRank = " << startRank << ", p = " << s.pValue << ", " << "ip = " << empiricalPValue << ", KS = " << KSScore << endl;

            // for Anderson-Darling
            // allPValues.push_back(s.pValue);

        }
    }

    // at a significance of 0.01, the KS-test succeeds if KSScore < 1.63 / sqrt(n) (for n large enough)
    if (sampleSize >= 20) {
        KSScore *= (sqrt((double) sampleSize) / 1.63);
    } else {
        KSScore = 0.0;
    }
    top.KSScore = KSScore;

/* (don't check fit for now

  // check fit
  if (!useDefault && KSScore > 1.5) {
    // bad fit, re-calculate using default
    
    // cerr << "REVERT TO DEFAULT KS=" << KSScore << " n= " << sampleSize << " Mean= " << mean << " StDev= " << stdev << endl;
    
    gaussianMean = mean; // 0.40;
    gaussianStDev = 0.08;
    
    for (int rank = 0; rank < numHits; rank++) {
      
      SpectraSTSimScores& s = m_candidates[rank]->getSimScoresRef();
      
      double sqrtDot = sqrt(s.dot);
      long double z = (long double)(sqrtDot - gaussianMean) / (long double)gaussianStDev;
      s.pValue = erfc(z * 0.70710678) * 0.5000000000;
    }
    
  }
*/

    for (int rank = 0; rank < numHits; rank++) {

        SpectraSTSimScores &s = m_candidates[rank]->getSimScoresRef();

        double fval = s.calcSP5Fval(numHits);

        s.hitsNum = sampleSize;
        s.hitsMean = mean;  // mean of the hits used to fit model, not all candidates
        s.hitsStDev = stdev;  // stdev of the hits used to fit model, not all candidates

        m_candidates[rank]->setSortKey(fval);


    }

    // ADScore = calcAndersonDarlingScore(allPValues);
    // cerr << sampleSize << '\t' << KSScore << '\t' << ADScore << endl;

    // Wenguang: remove the queries with bad fit; it also intends to remove all the query with extremely limited search space (below 50)-- because
    // these cases are fitted with the default paras, so that they can never be well fitted...
    // if (top.KSScore > 2.0) {
    //  top.fval = 0.0;
    // }


}

double SpectraSTSearch::calcAndersonDarlingScore(vector<double> &allPValues) {

    if (allPValues.empty()) return (0.0);

    int n = allPValues.size();
    double sum = 0.0;

    for (int rank = n - 1; rank >= 0; rank--) {

        int k = n - rank;

        double p = allPValues[k - 1];
        double p_mirror = allPValues[n - k];


        cerr << "p = " << p << "; p_mirror = " << p_mirror << "; k = " << k << endl;
        sum -= (double) (2 * k + 1) * (log(1 - p) + log(p_mirror)) / (double) n;

    }

    return (sum);

}


/*
// calcHitsStats - calculates such things as the mean and stdev of the dots of all candidates
void SpectraSTSearch::calcHitsStats() {
  
  unsigned int numHits;
  double mean;
  double stdev;
  
  if (m_candidates.empty()) {
    numHits = 0;
    mean = 0.0;
    stdev = 0.0;
  } else {
    
    double totalDot = 0;
    double totalSqDot = 0;
    for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
      SpectraSTSimScores& s = (*i)->getSimScoresRef();
      totalDot += s.dot;
      totalSqDot += s.dot * s.dot;
    }
    numHits = (unsigned int)(m_candidates.size());
    mean = totalDot / (double)numHits;
    stdev = totalSqDot / (double)numHits - mean * mean;
    if (stdev > 0.000001) stdev = sqrt(stdev);
  }
  for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
    SpectraSTSimScores& s = (*i)->getSimScoresRef();
    s.hitsNum = numHits;
    s.hitsMean = mean;
    s.hitsStDev = stdev;
    
  }
  	
}
*/

// detectHomologs - try to see if the lower hits are homologous (or identical) to the first one
void SpectraSTSearch::detectHomologs() {

    if (m_candidates.empty() || m_candidates.size() == 1) return;

    Peptide *topHit = m_candidates[0]->getEntry()->getPeptidePtr();

    if (!topHit) {
        // not a peptide. no notion of homology
        (m_candidates[0]->getSimScoresRef()).firstNonHomolog = 2;
        return;
    }

    bool homologFound = false;
    unsigned int curRank = 0;

    // go down the hit lists until hitting a nonhomologous hit
    do {
        homologFound = false;
        curRank++;

        Peptide *thisHit = m_candidates[curRank]->getEntry()->getPeptidePtr();

        if (!thisHit) {
            // not a peptide. definitely nonhomologous
            break;
        }

        int identity = 0;
        if (*topHit == *thisHit) {
            // identical!
            homologFound = true;
            //} else if (fabs(m_candidates[0]->getEntry()->getPrecursorMz() - m_candidates[curRank]->getEntry()->getPrecursorMz()) < 0.01) {
            //  // exactly same mass, consider as homolog (possibly target/decoy pair)
            //  homologFound = true;
        } else if ((topHit->stripped.length() > thisHit->stripped.length() && topHit->isSubsequence(*thisHit, true)) ||
                   (thisHit->stripped.length() < topHit->stripped.length() && thisHit->isSubsequence(*topHit, true))) {
            // one is subsequence of the other!
            homologFound = true;
        } else if (topHit->isHomolog(*thisHit, 0.7, identity)) {
            homologFound = true;

        }

    } while (homologFound && curRank < m_params.detectHomologs - 1 &&
             curRank < (unsigned int) (m_candidates.size()) - 1);

    // setting the field firstNonHomolog for all the homologs found
    for (unsigned int rank = 0; rank < curRank; rank++) {
        (m_candidates[rank]->getSimScoresRef()).firstNonHomolog = curRank + 1;
    }

}


// print - prints out the search result
void SpectraSTSearch::print() {


    if (m_candidates.empty() || (!(m_candidates[0]->passTopHitFvalThreshold()))) {
        // no hit (or no hit above threshold
        if (g_verbose) {
            cout << " DONE! Top hit: NO_MATCH";
            cout.flush();
        }
        if (m_params.hitListExcludeNoMatch) {
            // told to exclude all NO_MATCH's from the final output, so just return
            return;
        }
    }

    // set the assumedCharge to that of the top hit, if any
    int assumedCharge = 0;
    if (!m_candidates.empty() && m_candidates[0]->passTopHitFvalThreshold()) {
        assumedCharge = m_candidates[0]->getEntry()->getCharge();
    }


    string name = m_query->getName();
    double precursorMz = m_query->getPrecursorMz();
    double rt = m_query->getRetentionTime();
    // prints all the query information
    m_output->printStartQuery(name, precursorMz, assumedCharge, rt);


    if (m_candidates.empty() || !(m_candidates[0]->passTopHitFvalThreshold())) {
        // no hit or no hit above threshold
        SpectraSTSimScores emptyScore;

        m_output->printHit(name, 1, NULL, emptyScore);

    } else {
        // print the top hit!
        if (g_verbose) {
            cout << " DONE! Top hit: " << m_candidates[0]->getEntry()->getFullName() << " (Fval = "
                 << m_candidates[0]->getSortKey() << ")";
            cout.flush();
        }

        m_output->printHit(name, 1, m_candidates[0]->getEntry(), m_candidates[0]->getSimScoresRef());

        if (m_params.hitListShowMaxRank > 1) {
            // told to print the lower hits too

            unsigned int firstNonHomolog = (m_candidates[0]->getSimScoresRef()).firstNonHomolog;

            for (unsigned int rank = 2; rank <= (unsigned int) (m_candidates.size()); rank++) {
                if ((rank <= m_params.hitListShowMaxRank && m_candidates[rank - 1]->passLowerHitsFvalThreshold()) ||
                    (m_params.hitListShowHomologs && rank < firstNonHomolog)) {
                    m_output->printHit(name, rank, m_candidates[rank - 1]->getEntry(),
                                       m_candidates[rank - 1]->getSimScoresRef());
                } else {
                    // below threshold now, stop printing
                    break;
                }
            }
        }

    }

    // print the closing tags to finish up
    m_output->printEndQuery(name);


}

// isLikelyGood. Simply returns if the F value is above 0.5, indicating a likely good hit (not always!)
bool SpectraSTSearch::isLikelyGood() {
    if (m_candidates.size() > 0 && m_candidates[0]->getSortKey() >= 0.5) {
        return (true);
    } else {
        return (false);
    }
    return (false);
}

// isLikelyBad. Simply returns if the F value is below 0.2.
bool SpectraSTSearch::isLikelyBad() {
    if (m_candidates.size() > 0 && m_candidates[0]->getSortKey() < 0.2) {
        return (true);
    } else {
        return (false);
    }
    return (false);
}

// isDecoy. returns if the hit is a decoy entry
bool SpectraSTSearch::isDecoy(unsigned int rank) {
    if ((unsigned int) (m_candidates.size()) >= rank) {
        string spec("");
        string rem("");
        return ((m_candidates[rank - 1]->getEntry()->getOneComment("Spec", spec) && spec == "Decoy") ||
                (m_candidates[rank - 1]->getEntry()->getOneComment("Remark", rem) && rem.substr(0, 5) == "DECOY"));
    } else {
        return (false);
    }
    return (false);
}

bool SpectraSTSearch::isSingleton(unsigned int rank) {
    if ((unsigned int) (m_candidates.size()) >= rank) {
        string nreps("");
        return (m_candidates[rank - 1]->getEntry()->getNrepsUsed() == 1);
    }
    return (false);
}

double SpectraSTSearch::upperIncompleteGamma(double x, double a) {


    double sum = 0.0;
    double term = 1.0 / a;
    int n = 1;

    while (term > 0) {
        sum += term;
        term *= (x / (a + (double) n));
        n++;
    }

    return (pow(x, a) * exp(-x) * sum);


}

