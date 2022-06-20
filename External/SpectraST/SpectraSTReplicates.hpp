#ifndef SPECTRASTREPLICATES_HPP_
#define SPECTRASTREPLICATES_HPP_

#include "SpectraSTLibEntry.hpp"
#include "SpectraSTCreateParams.hpp"
#include "SpectraSTDenoiser.hpp"
#include "SpectraSTConstants.hpp"
#include "Peptide.hpp"

#include <map>
#include <string>
#include <vector>

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

/* Class: SpectraSTReplicates
 * 
 * Manages a set of "replicate spectra" -- spectra that are spectrally similar and are identified to the same peptide ion.
 * Performs function like weighting, consensus creation, etc. 
 */

using namespace std;

// a type of one replicate
typedef struct _replicate {
    SpectraSTLibEntry *entry; // points to the entry itself
    unsigned int numUsed;
    unsigned int numTotal;
    int status; // 1 = used, 0 = not used, neg. numbers = cluster number
    double prob;
    double sn;
    double xcorr;
} Replicate;

class SpectraSTReplicates {

public:
    SpectraSTReplicates(vector<SpectraSTLibEntry *> &entries, SpectraSTCreateParams &params,
                        vector<SpectraSTDenoiser *> *denoisers = NULL);

    virtual ~SpectraSTReplicates();

    SpectraSTLibEntry *makeConsensusSpectrum();

    SpectraSTLibEntry *findBestReplicate();

    void setPlotPath(string plotPath) { m_plotPath = plotPath; }

    void setRecordRawSpectra(bool value) { m_recordRawSpectra = value; }

    static double getMedian(vector<double> &v);

private:

    SpectraSTCreateParams &m_params;
    Peptide *m_pep; // The Peptide object representing the identification of these replicates. NOT a property of this class
    string m_name;
    vector<Replicate> m_reps;
    unsigned int m_numUsed;
    unsigned int m_numTotal;

    bool m_missingXCorr; // whether some replicates don't have xcorr values

    string m_plotPath;

    bool m_recordRawSpectra;

    vector<SpectraSTDenoiser *> *m_denoisers;

    // hashes containing sequence search, sample source, and instrument information (parsed from Comment fields)
    map<char, map<string, pair<double, double> > *> m_seqs;
    map<string, pair<unsigned int, unsigned int> > m_samples;
    map<string, pair<unsigned int, unsigned int> > m_instruments;


    void addEntry(SpectraSTLibEntry *entry);

    bool removeDissimilarReplicates();

    void aggregateStats(SpectraSTLibEntry *final, SpectraSTPeakList *pl, vector<Replicate *> &usedReps,
                        bool isConsensus = true);

    void plotSingle(SpectraSTLibEntry *single);

    void
    plotConsensus(SpectraSTLibEntry *consensus, SpectraSTPeakList *consensusPeakList, vector<Replicate *> &usedReps);

    void
    processConsensus(SpectraSTLibEntry *consensus, SpectraSTPeakList *consensusPeakList, vector<Replicate *> &usedReps);

    void processBestReplicate(SpectraSTLibEntry *best, vector<Replicate *> &usedReps, bool isRaw);

    void processSingle(SpectraSTLibEntry *single, bool isRaw);

    static bool sortReplicatesByWeight(Replicate a, Replicate b);

};

#endif /*SPECTRASTREPLICATES_HPP_*/
