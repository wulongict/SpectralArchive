#include "SpectraSTReplicates.hpp"
#include "SpectraSTPeakList.hpp"
#include "SpectraSTSpLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <math.h>
#include <sstream>
#include <algorithm> // for std::sort

// #define PLOT_ALL_CONSENSUS 
// #define PLOT_BAD_CONSENSUS 1

//static double g_fracUnassignedTop20;
//static double g_fracUnassignedTop5;
//static int    g_numUnassignedTop20;

extern SpectraSTLog *g_log;

//extern double g_retained;
//extern double g_retainedSq;
//extern int g_retainedCount;

SpectraSTReplicates::SpectraSTReplicates(vector<SpectraSTLibEntry *> &entries, SpectraSTCreateParams &params,
                                         vector<SpectraSTDenoiser *> *denoisers) :
        m_params(params),
        m_reps(),
        m_numUsed(0),
        m_numTotal(0),
        m_seqs(),
        m_samples(),
        m_instruments(),
        m_missingXCorr(false),
        m_recordRawSpectra(false),
        m_denoisers(denoisers) {

    m_recordRawSpectra = m_params.recordRawSpectra;

    if (!entries.empty()) {
        m_pep = entries[0]->getPeptidePtr();
        m_name = entries[0]->getName();

    }


    for (vector<SpectraSTLibEntry *>::iterator i = entries.begin(); i != entries.end(); i++) {
        addEntry(*i);
    }

    //  m_calcXCorr = true;
    if (m_missingXCorr && m_params.replicateWeight == "XCORR") {

        // at least one replicate does not have a SEQUEST-calculated xcorr,
        // need to calculate xcorr anew for all replicates (since our xcorr function cannot
        // reproduce SEQUEST's value exactly, making it problematic to compare)

        for (vector<Replicate>::iterator r = m_reps.begin(); r != m_reps.end(); r++) {
            SpectraSTPeakList *pl = r->entry->getPeakList();
            double calcXCorr = pl->calcXCorr();
            r->xcorr = calcXCorr;
            double xcorrCubed = r->xcorr * r->xcorr * r->xcorr;
            if (xcorrCubed > 100.0) xcorrCubed = 100.0;
            pl->setWeight(xcorrCubed);
        }
    }

    // sort(m_reps.begin(), m_reps.end(), SpectraSTReplicates::sortReplicatesByWeight);


}

// destructor
SpectraSTReplicates::~SpectraSTReplicates() {

    for (map<char, map<string, pair<double, double> > *>::iterator j = m_seqs.begin(); j != m_seqs.end(); j++) {
        if (j->second) {
            delete (j->second);
        }
    }

}

// addEntry - adds a replicate to the list, performs weighting 
void SpectraSTReplicates::addEntry(SpectraSTLibEntry *entry) {

    Replicate r;
    r.entry = entry;
    r.sn = 0.0;
    r.xcorr = 0.0;
    r.prob = 0.0;
    r.numUsed = 1;
    r.numTotal = 1;

    // hijack -- removes anything whose massdiff is not 0
    //string mdStr("");
    //if (entry->getOneComment("MassDiff", mdStr)) {
    //  double massdiff = atof(mdStr.c_str());
    //  int massdiffInt = (int)(massdiff + 0.5);
    //  if (massdiffInt != 0) return;
    //}
    // END hijack

    // parse out num reps (The replicate can be a previous consensus spectrum, too)
    string numRepStr("");
    if (entry->getOneComment("Nreps", numRepStr)) {
        string::size_type slashPos = numRepStr.find('/', 0);
        r.numUsed = atoi(numRepStr.substr(0, slashPos).c_str());
        r.numTotal = atoi(numRepStr.substr(slashPos + 1).c_str());
    }

    string specType("Raw");
    entry->getOneComment("Spec", specType);

    // add to the overall numUsed and numTotal counts
    m_numUsed += r.numUsed;
    m_numTotal += r.numTotal;

    // get probability
    r.prob = entry->getProb();

    // parse out xcorr - only available from SEQUEST of course
    string xcorrStr("");
    if (entry->getOneComment("XCorr", xcorrStr)) {
        r.xcorr = atof(xcorrStr.c_str());
    } else {
        r.xcorr = 1.0;
        m_missingXCorr = true;
    }

    // parse out S/N, if not available, calculates it
    // note that previous consensus spectra will already have an SN field, so that
    // we don't have to calcSignalToNoise for them, as it won't work (They already had noise peaks removed).
    string snStr("");
    if (entry->getOneComment("SN", snStr)) {
        r.sn = atof(snStr.c_str());
    } else {
        r.sn = entry->getPeakList()->calcSignalToNoise();
    }

    // parse out precursor intensity.
    string preIntStr("");
    double preInt = 10000.0; // meaningless default value. avoid this option if precursor intensity is not consistently present
    if (entry->getOneComment("PrecursorIntensity", preIntStr)) {
        preInt = atof(preIntStr.c_str());
        if (preInt < 0.1) {
            if (entry->getOneComment("TotalIonCurrent", preIntStr)) {
                preInt = atof(preIntStr.c_str());
            }
        }
    }
    entry->getPeakList()->setPrecursorIntensity(preInt);

    // weighting. there are 4 options:
    // XCORR: weight = (XCorr)^3, capped at 100
    // PROB: weight = (-log10(1.00001 - Prob) * 20, capped at 100. (This emphasizes the region near 1.0)
    // SN: weight = sqrt(SN) * 5.0, capped at 100 (since SN is capped at 400).
    // INTP: weight = sqrt(precursor intensity)
    // all weights are multiplied by numUsed; therefore a previous consensus spectrum of multiple spectra will
    // count proportionally more
    double weight = 1.0;
    if (m_params.replicateWeight == "XCORR") {
        double xcorrCubed = r.xcorr * r.xcorr * r.xcorr;
        if (xcorrCubed > 100.0) xcorrCubed = 100.0;
        weight = r.numUsed * xcorrCubed;
    } else if (m_params.replicateWeight == "PROB") {
        weight = r.numUsed * (-log10(1.00001 - r.prob)) * 20.0;
    } else if (m_params.replicateWeight == "SN") {
        weight = r.numUsed * sqrt(r.sn) * 5.0;
    } else if (m_params.replicateWeight == "INTP") {
        weight = r.numUsed * sqrt(preInt);
    } else if (m_params.replicateWeight == "NONE") {
        weight = r.numUsed;
    } else {
        weight = r.numUsed;
    }

    entry->getPeakList()->setWeight(weight);

    // parse out the Se and Sample fields of the Comment and add them to the hashes
    entry->getInstrumentInfo(m_instruments, "TOTAL_ONLY");
    entry->getSampleInfo(m_samples, "TOTAL_ONLY");

    r.status = 1;

    m_reps.push_back(r);


}

// findBestReplicate - finds the best replicate based on weight
SpectraSTLibEntry *SpectraSTReplicates::findBestReplicate() {

    if (m_reps.empty()) {
        stringstream logss;
        logss << "Removing spectrum: " << m_name << " (NO GOOD REPLICATE)";
        g_log->log("CREATE", logss.str());
        return (NULL);
    }

    // if minimumNumReplicates is specified, it will only pick a best replicates from
    // enough replicates. if not enough replicates are present for a peptide ion,
    // no "best-replicate" spectrum will be included for that peptide ion in the final library
    if (m_numUsed < m_params.minimumNumReplicates) {
        stringstream logss;
        logss << "Removing spectrum: " << m_name << " (NOT ENOUGH REPLICATES)";
        g_log->log("CREATE", logss.str());
        return (NULL);
    }

    // sort replicates by weight now
    sort(m_reps.begin(), m_reps.end(), SpectraSTReplicates::sortReplicatesByWeight);

    // m_reps are already sorted by weight, so m_reps[0] is the best replicate
    SpectraSTLibEntry *best = m_reps[0].entry;
    bool isRaw = (m_reps[0].numUsed == 1);

    // also check for dissimilar replicates. If the "best" replicate is considerably different from everyone else
    // then it's probably not so good after all. Will actually picked the best replicate that is in the largest similarity cluster
    if (m_params.removeDissimilarReplicates) {
        if (!removeDissimilarReplicates()) {
            stringstream logss;
            logss << "Removing spectrum: " << m_name << " (NOT ENOUGH SIMILAR REPLICATES)";
            g_log->log("CREATE", logss.str());
            return (NULL);
        } else {
            for (vector<Replicate>::iterator r = m_reps.begin(); r != m_reps.end(); r++) {
                if (r->status == 1) {
                    best = r->entry;
                    isRaw = (r->numUsed == 1);
                    break;
                }
            }
        }
    }

    m_numUsed = 0;
    vector<Replicate *> usedReps;
    for (vector<Replicate>::iterator rep = m_reps.begin(); rep != m_reps.end(); rep++) {

        if (rep->status == 1) {
            m_numUsed += rep->numUsed;

            // update the "USED" statistics of the Se, Sample and Instrument fields
            rep->entry->getSeqInfo(m_seqs);
            rep->entry->getInstrumentInfo(m_instruments, "USED_ONLY");
            rep->entry->getSampleInfo(m_samples, "USED_ONLY");

            usedReps.push_back(&(*rep));

        }
    }

    if (m_numUsed == 1) {
        stringstream logss;
        logss << "Keeping single spectrum: " << m_name;
        g_log->log("CREATE", logss.str());

    } else {
        stringstream logss;
        logss << "Selecting best replicate among " << m_numUsed << " (of " << m_numTotal << ") replicates: " << m_name;
        g_log->log("CREATE", logss.str());

    }

    processBestReplicate(best, usedReps, isRaw);

    return (best);
}

// makeConsensusSpectrum - creates a consensus spectrum of the replicates
SpectraSTLibEntry *SpectraSTReplicates::makeConsensusSpectrum() {

    if (m_reps.empty()) {
        stringstream logss;
        logss << "Removing spectrum: " << m_name << " (NO GOOD REPLICATE)";
        g_log->log("CREATE", logss.str());
        return (NULL);
    }

    if (m_numUsed < m_params.minimumNumReplicates) {
        stringstream logss;
        logss << "Removing spectrum: " << m_name << " (NOT ENOUGH REPLICATES)";
        g_log->log("CREATE", logss.str());
        return (NULL);
    }

    // only one replicate (which could be a consensus of replicates previously!) -
    // this will be our "consensus"
    if (m_reps.size() == 1) {
        SpectraSTLibEntry *single = m_reps[0].entry;

        if (m_numUsed == 1) {
            // a real single
            processSingle(single, true);

            stringstream logss;
            logss << "Keeping single spectrum: " << m_name;
            g_log->log("CREATE", logss.str());

        } else {
            // actually a spectrum from a consensus of many previously
            processSingle(single, false);

            stringstream logss;
            logss << "Keeping previous consensus spectrum: " << m_name;
            g_log->log("CREATE", logss.str());

        }
        return (single);
    }

    // sort replicates by weight now
    sort(m_reps.begin(), m_reps.end(), SpectraSTReplicates::sortReplicatesByWeight);

    // remove dissimilar replicates
    if (m_params.removeDissimilarReplicates) {
        if (!removeDissimilarReplicates()) {
            stringstream logss;
            logss << "Removing spectrum: " << m_name << " (NOT ENOUGH SIMILAR REPLICATES)";
            g_log->log("CREATE", logss.str());
            return (NULL);
        }
    }

    // reset num used counter; go down the list of replicates to select replicates that we will use
    m_numUsed = 0;
    vector<Replicate *> usedReps;
    vector<SpectraSTPeakList *> pls;
    bool enough = false;
    for (vector<Replicate>::iterator rep = m_reps.begin(); rep != m_reps.end(); rep++) {

        if (enough) {
            rep->status = 0;
        }

        if (rep->status == 1) {
            m_numUsed += rep->numUsed;

            // update the "USED" statistics of the Se, Sample and Instrument fields
            rep->entry->getSeqInfo(m_seqs);
            rep->entry->getInstrumentInfo(m_instruments, "USED_ONLY");
            rep->entry->getSampleInfo(m_samples, "USED_ONLY");

            usedReps.push_back(&(*rep));
            pls.push_back(rep->entry->getPeakList());

            if (m_numUsed >= m_params.maximumNumReplicates) {
                // if we are capping the number of replicates used, and we already hit the maximum
                // we are done, not need to add more replicates
                enough = true;
            }
        }
    }

    if (usedReps.size() == 1) {
        // single (we got here because removeDissimilarReplicates() got rid of all but one replicate spectrum)
        // numTotal might have been increased, need to set Nreps
        stringstream numRepss;
        numRepss << m_numUsed << '/' << m_numTotal;
        usedReps[0]->entry->setOneComment("Nreps", numRepss.str());

        if (m_numUsed == 1) {
            processSingle(usedReps[0]->entry, true);

            stringstream logss;
            logss << "Keeping single spectrum: " << m_name;
            g_log->log("CREATE", logss.str());
        } else {
            processSingle(usedReps[0]->entry, false);

            stringstream logss;
            logss << "Keeping previous consensus spectrum: " << m_name;
            g_log->log("CREATE", logss.str());
        }
        return (usedReps[0]->entry);

    } else {

        // now we actually have more than one spectra that we need to "consensus" together
        stringstream logss;
        logss << "Creating consensus of " << m_numUsed << " (of " << m_numTotal << ") replicates: " << m_name;
        g_log->log("CREATE", logss.str());

        // calling the consensus-forming constructor
        SpectraSTPeakList *cpl = new SpectraSTPeakList(pls, usedReps[0]->entry->getPeptidePtr(), m_numUsed,
                                                       m_params.peakQuorum, m_params.maximumNumPeaksUsed, m_denoisers,
                                                       m_params.keepRawIntensities);

        // use top entry as holder for consensus -- note that consensus->getPeakList() still points to the best replicate;
        // it has not been set to the actual consensus spectrum (cpl) yet
        SpectraSTLibEntry *consensus = usedReps[0]->entry;

        processConsensus(consensus, cpl, usedReps);

        // up tp this point, consensus->getPeakList() stills holds the best replicate (which is also pointed to by the usedReps vector)
        // now that we are done with evaluating the replicates, we finally can overwrite the best replicate with the consensus spectrum
        consensus->setPeakList(cpl);

        return (consensus);
    }
}

// processConsensus - finish up a consensus creation
void SpectraSTReplicates::processConsensus(SpectraSTLibEntry *consensus, SpectraSTPeakList *consensusPeakList,
                                           vector<Replicate *> &usedReps) {


    string specType("");

    // set comments
    specType = "Consensus";
    consensus->setOneComment("Spec", "Consensus");

    consensus->setStatus("Normal"); // reset

    consensus->setInstrumentInfo(m_instruments);
    consensus->setSeqInfo(m_seqs);
    consensus->setSampleInfo(m_samples);

    string raw("");
    if (consensus->getOneComment("RawSpectrum", raw)) {
        consensus->setOneComment("BestRawSpectrum", raw);
    }

    if (m_recordRawSpectra) {
        // NOTE: no duplicate removal here. if two identical raw spectrum names show up, they will both be included
        string rawSpectra("");
        for (vector<Replicate *>::iterator r = usedReps.begin(); r != usedReps.end(); r++) {
            string spectra = (*r)->entry->getRawSpectra();
            if (spectra.empty()) continue;
            if (!rawSpectra.empty()) rawSpectra += ",";
            rawSpectra += spectra;
        }
        consensus->setOneComment("RawSpectra", rawSpectra);
    }

    consensus->deleteOneComment("RawSpectrum");

    stringstream numRepss;
    numRepss << m_numUsed << '/' << m_numTotal;
    consensus->setOneComment("Nreps", numRepss.str());

    if (m_params.maximumNumPeaksKept > 0) {
        double retained = 1.0;
        retained = consensusPeakList->simplify(m_params.maximumNumPeaksKept, 99999.0, false, false, 0.0);
    }

    aggregateStats(consensus, consensusPeakList, usedReps, true);

    if (m_params.plotSpectra == "ALL" || m_params.plotSpectra == consensus->getStatus() ||
        m_params.plotSpectra == specType) {
        plotConsensus(consensus, consensusPeakList, usedReps);
    }

}

// processBestReplicate - finish up a best replicate selection
void SpectraSTReplicates::processBestReplicate(SpectraSTLibEntry *best, vector<Replicate *> &usedReps, bool isRaw) {


    best->setOneComment("Spec", "BestReplicate");

    best->setInstrumentInfo(m_instruments);
    best->setSeqInfo(m_seqs);
    best->setSampleInfo(m_samples);

    stringstream numRepss;
    numRepss << m_numUsed << '/' << m_numTotal;
    best->setOneComment("Nreps", numRepss.str());

    processSingle(best, isRaw);

    aggregateStats(best, best->getPeakList(), usedReps, false);
}

// processSingle - finish up including a single spectrum
void SpectraSTReplicates::processSingle(SpectraSTLibEntry *single, bool isRaw) {

    if (m_numUsed == 1) {
        single->setOneComment("Spec", "Single");
    }
    single->setStatus("Normal"); // reset

    string dummy("");
    if (!isRaw && !(single->getOneComment("Probcorr", dummy))) {
        // enforce the peak quorum  -- for this only entry (note that this is necessary
        // in the case where the entry is previously created with a lower quorum)
        // hack - don't do this for NIST spectra, which have this Probcorr field, since the
        // quorum specification has different format
        int minNumRepWithPeak = (int) ((double) m_numUsed * m_params.peakQuorum - 0.00001) + 1;
        if (minNumRepWithPeak < 1) minNumRepWithPeak = 1;
        single->getPeakList()->removeInquoratePeaks(minNumRepWithPeak);
    }

    if (isRaw && m_denoisers && (*m_denoisers)[0]->isFilterReady()) {
        // select the denoiser of the right charge
        unsigned int charge = single->getCharge();
        if (charge > MAX_CHARGE) charge = MAX_CHARGE;
        if (!((*m_denoisers)[charge]->isFilterReady())) {
            // not ready, maybe lack of training data for that charge. Use the denoiser trained with all data
            charge = 0;
        }
        (*m_denoisers)[charge]->filter(single->getPeakList(), m_params.maximumNumPeaksKept,
                                       m_params.denoiserMinimumSignalProb);
    }

    if (m_params.maximumNumPeaksKept > 0) {
        double retained = 1.0;
        retained = single->getPeakList()->simplify(m_params.maximumNumPeaksKept, 99999.0, false, false, 0.0);
    }

    single->annotatePeaks();

    //   single->getPeakList()->flattenAllPeaks();
    //  single->getPeakList()->removeNoncanonicalPeaks();

    string specType("");
    if (m_params.plotSpectra == "ALL" || m_params.plotSpectra == single->getStatus() ||
        (single->getOneComment("Spec", specType) && !m_params.plotSpectra.empty() &&
         m_params.plotSpectra == specType)) {
        plotSingle(single);
    }

}

// removeDissimilarReplicates - uses a lazy clustering method to remove replicates that are not similar to everyone else.
// Idea: Start from the top-weighted (highest-quality) replicate, put it in cluster zero. For every replicate following it, compare it
// to all the existing clusters starting from cluster zero; if similar enough, add to the cluster, else create its own cluster.
// At the end, the cluster having the most replicates are retained, everything else is thrown out. 
bool SpectraSTReplicates::removeDissimilarReplicates() {

    vector<SpectraSTLibEntry *> clusters;
    vector<unsigned int> clusterSize;

    for (vector<Replicate>::iterator rep = m_reps.begin(); rep != m_reps.end(); rep++) {

        SpectraSTPeakList *pl = rep->entry->getPeakList();
        int mdInt = rep->entry->getMassDiffInt();

        for (vector<SpectraSTLibEntry *>::size_type cl = 0; cl < clusters.size(); cl++) {

            // hijack -- do not cluster replicates or different massdiff integer values
            // if (mdInt != clusters[cl]->getMassDiffInt()) {
            // do not attempt to cluster
            //	continue;
            //      }
            // end hijack

            double dot = clusters[cl]->getPeakList()->compare(pl);
            if (dot >= 0.9999) {
                // fishy. probably identical spectra being included twice!
                string query1("");
                string query2("");
                if (rep->entry->getOneComment("RawSpectrum", query1) &&
                    clusters[cl]->getOneComment("RawSpectrum", query2)) {
                    g_log->log("CREATE", "Identical replicates: " + query1 + " and " + query2);
                }
                // will not use this spectrum
                rep->status = 0;
                break;
            }
            if (dot >= 0.6) {
                // similar enough to cluster cl. add to it and we're done with this replicate.
                rep->status = -((int) cl + 1);
                clusterSize[cl] += rep->numUsed;
                break;
            }
        }

        if (rep->status > 0) {
            // still doesn't belong to any cluster. start a new one
            rep->status = -((int) (clusters.size()) + 1);
            clusters.push_back(rep->entry);
            clusterSize.push_back(rep->numUsed);
        }
    }

    unsigned int largestClusterSize = 0;
    vector<SpectraSTLibEntry *>::size_type bestCluster = 0;
    for (vector<SpectraSTLibEntry *>::size_type cl = 0; cl < clusters.size(); cl++) {
        if (clusterSize[cl] > largestClusterSize) {
            largestClusterSize = clusterSize[cl];
            bestCluster = cl;
        }
    }

    // check to see if this largest cluster becomes too small to create a consensus/best-replicate.
    // if so, we will quit and return false.
    if (largestClusterSize < m_params.minimumNumReplicates) {
        return (false);
    }

    m_numUsed = largestClusterSize;

    // set the status of all replicates in the "best" cluster to 1 (used)
    // everything else will have status 0 (not used)
    for (vector<Replicate>::iterator rep = m_reps.begin(); rep != m_reps.end(); rep++) {
        if (rep->status == -((int) bestCluster + 1)) {
            rep->status = 1;
        } else {
            rep->status = 0;
        }
    }
    return (true);

}

// aggregateStats - calculates some statistics of the consensus spectrum created
void SpectraSTReplicates::aggregateStats(SpectraSTLibEntry *final, SpectraSTPeakList *pl, vector<Replicate *> &usedReps,
                                         bool isConsensus) {

    double sumDot = 0.0;
    double sumSqDot = 0.0;
    unsigned int problematicDot = 0;

    double maxSN = 0.0;
    double maxProb = 0.0;
    double minProb = 1.0;
    double aveProb = 0.0;
    double adjProb = 1.0;
    double maxXCorr = 0.0;

    unsigned int sumNumPeaks = 0;
    unsigned int sumSqNumPeaks = 0;
    double sumFracAssigned = 0.0;
    unsigned int numAnnotated = 0;

    double maxRt = -1.0;
    double minRt = 999999.9;
    double maxIRT = -1.0;
    double minIRT = 999999.9;
    double maxImb = -1.0;
    double minImb = 999999.9;

    // double sumRt = 0.0;
    // unsigned int sumNumUsedWithRt = 0;
    vector<double> RTs;
    vector<double> iRTs;
    vector<double> imbs;

//  SpectraSTPeakList* pl = final->getPeakList();

    double sumTic = 0.0;
    double sumPrecInt = 0.0;
    double sumOrigMaxInt = 0.0;

    unsigned int ticCount = 0;
    unsigned int precIntCount = 0;
    unsigned int origMaxIntCount = 0;

    // loop over all replicates for various evaluation data
    for (vector<Replicate *>::iterator r = usedReps.begin(); r != usedReps.end(); r++) {

        // check dot product between consensus and each replicate
        double dot = pl->compare((*r)->entry->getPeakList());
//    cerr << dot << endl;
        sumDot += dot;
        sumSqDot += dot * dot;
        if (dot < 0.7) {
            problematicDot++;
        }

        // check the number of peaks (and how many are assigned) in the replicates
        unsigned int numPeaks = (*r)->entry->getPeakList()->getNumPeaks();
        sumNumPeaks += numPeaks;
        sumSqNumPeaks += numPeaks * numPeaks;
        unsigned int numAssignedPeaks = (*r)->entry->getPeakList()->getNumAssignedPeaks();
        sumFracAssigned += (double) numAssignedPeaks / (double) numPeaks;
        numAnnotated++;

        // aggregate the quality measures of all replicates
        if ((*r)->sn > maxSN) maxSN = (*r)->sn;
        if ((*r)->xcorr > maxXCorr) maxXCorr = (*r)->xcorr;
        if ((*r)->prob > maxProb) maxProb = (*r)->prob;
        if ((*r)->prob < minProb) minProb = (*r)->prob;
        aveProb += (*r)->prob;
        adjProb *= (1 - (*r)->prob);

        // aggregate the retention times and iRTs
        double rt = 0.0;
        string rtstr("");
        if ((*r)->entry->getOneComment("RetentionTime", rtstr)) {
            string::size_type commaPos = 0;
            double rtmedian = atof(nextToken(rtstr, commaPos, commaPos, ", \t\r\n").c_str());
            double rtmax = atof(nextToken(rtstr, commaPos + 1, commaPos, ", \t\r\n").c_str());
            double rtmin = atof(nextToken(rtstr, commaPos + 1, commaPos, ", \t\r\n").c_str());

            RTs.insert(RTs.end(), (*r)->numUsed, rtmedian);

            //      sumRt += rtmedian * (*r)->numUsed;
            //      sumNumUsedWithRt += (*r)->numUsed;
            if (rtmax > maxRt) maxRt = rtmax;
            if (rtmin < minRt) minRt = rtmin;
        }

        string iRTStr("");
        if ((*r)->entry->getOneComment("iRT", iRTStr)) {
            string::size_type commaPos = 0;
            double iRTMedian = atof(nextToken(iRTStr, commaPos, commaPos, ", \t\r\n").c_str());
            double iRTMax = atof(nextToken(iRTStr, commaPos + 1, commaPos, ", \t\r\n").c_str());
            double iRTMin = atof(nextToken(iRTStr, commaPos + 1, commaPos, ", \t\r\n").c_str());

            iRTs.insert(iRTs.end(), (*r)->numUsed, iRTMedian);

            if (iRTMax > maxIRT) maxIRT = iRTMax;
            if (iRTMin < minIRT) minIRT = iRTMin;
        }

        string imbstr("");
        if ((*r)->entry->getOneComment("IonMobility", imbstr)) {
            string::size_type commaPos = 0;
            double imbMedian = atof(nextToken(imbstr, commaPos, commaPos, ", \t\r\n").c_str());
            double imbMax = atof(nextToken(imbstr, commaPos + 1, commaPos, ", \t\r\n").c_str());
            double imbMin = atof(nextToken(imbstr, commaPos + 1, commaPos, ", \t\r\n").c_str());

            imbs.insert(imbs.end(), (*r)->numUsed, imbMedian);

            //      sumRt += rtmedian * (*r)->numUsed;
            //      sumNumUsedWithRt += (*r)->numUsed;
            if (imbMax > maxImb) maxImb = imbMax;
            if (imbMin < minImb) minImb = imbMin;
        }


        string ticStr("");
        if ((*r)->entry->getOneComment("TotalIonCurrent", ticStr)) {
            sumTic += (atof(ticStr.c_str()) * (*r)->numUsed);
            ticCount += (*r)->numUsed;
        }

        string precIntStr("");
        if ((*r)->entry->getOneComment("PrecursorIntensity", precIntStr)) {
            sumPrecInt += (atof(precIntStr.c_str()) * (*r)->numUsed);
            precIntCount += (*r)->numUsed;
        }

        double origMaxInt = 0.0;
        string origMaxIntStr("");
        if ((*r)->entry->getOneComment("OrigMaxIntensity", origMaxIntStr)) {
            sumOrigMaxInt += (atof(origMaxIntStr.c_str()) * (*r)->numUsed);
            origMaxIntCount += (*r)->numUsed;
        }

    }

    // record replicates' number of peaks
    double meanNumPeaks = (double) sumNumPeaks / (double) (usedReps.size());
    double stDevNumPeaks = 0.0;
    double varNumPeaks = (double) sumSqNumPeaks / (double) (usedReps.size()) - meanNumPeaks * meanNumPeaks;
    if (varNumPeaks > 0.0001) stDevNumPeaks = sqrt(varNumPeaks);
    stringstream repPeakss;
    repPeakss.precision(1);
    repPeakss << fixed << showpoint << meanNumPeaks << '/' << stDevNumPeaks;
    final->setOneComment("RepNumPeaks", repPeakss.str());

    // record replicates' dot products to final (consensus)
    double meanDot = sumDot / (double) (usedReps.size());
    double stDevDot = 0.0;
    double varDot = sumSqDot / (double) (usedReps.size()) - meanDot * meanDot;
    if (varDot > 0.0000001) stDevDot = sqrt(varDot);
    stringstream dotFinalss;
    dotFinalss.precision(2);
    dotFinalss << fixed << showpoint << meanDot << ',' << stDevDot << ';' << problematicDot << '/' << usedReps.size();
    if (isConsensus) {
        final->setOneComment("DotConsensus", dotFinalss.str());
    } else {
        final->setOneComment("DotFinal", dotFinalss.str());
    }

    // record fraction assigned numbers
    double meanFracAssigned = 0.0;
    if (numAnnotated > 0) {
        meanFracAssigned = sumFracAssigned / (double) numAnnotated;
        stringstream repAssignedss;
        repAssignedss.precision(3);
        repAssignedss << fixed << showpoint << meanFracAssigned;
        final->setOneComment("RepFracAssignedPeaks", repAssignedss.str());
    }

    stringstream finalAssignedss;
    finalAssignedss.precision(3);
    finalAssignedss << fixed << showpoint << (double) (pl->getNumAssignedPeaks()) / (double) (pl->getNumPeaks());
    if (isConsensus) {
        final->setOneComment("ConsFracAssignedPeaks", finalAssignedss.str());
    } else {
        final->setOneComment("FinalFracAssignedPeaks", finalAssignedss.str());
    }

    // recalculate xrea (useless, since consensus typically has horrible xrea after de-noising)
    // stringstream xreass;
    // xreass.precision(3);
    // xreass << fixed << pl->calcXrea(true);
    // final->setOneComment("Xrea", xreass.str());

    // record retention time
    if (maxRt >= 0.0) {
        stringstream rtss;
        rtss.precision(1);
        rtss << fixed << showpoint << SpectraSTReplicates::getMedian(RTs)
             << ','; // sumRt / (double)sumNumUsedWithRt << ',';
        rtss << fixed << showpoint << maxRt << ',';
        rtss << fixed << showpoint << minRt;
        final->setOneComment("RetentionTime", rtss.str());
    }
    if (maxIRT >= 0.0) {
        stringstream irtss;
        irtss.precision(1);
        irtss << fixed << showpoint << SpectraSTReplicates::getMedian(iRTs)
              << ','; // sumRt / (double)sumNumUsedWithRt << ',';
        irtss << fixed << showpoint << maxIRT << ',';
        irtss << fixed << showpoint << minIRT;
        final->setOneComment("iRT", irtss.str());
    }
    if (maxImb >= 0.0) {
        stringstream imbss;
        imbss.precision(6);
        imbss << fixed << showpoint << SpectraSTReplicates::getMedian(imbs) << ',';
        imbss << fixed << showpoint << maxImb << ',';
        imbss << fixed << showpoint << minImb;
        final->setOneComment("IonMobility", imbss.str());
    }

    if (ticCount > 0) {
        stringstream ticss;
        ticss.precision(2);
        ticss << sumTic / (double) (ticCount);
        final->setOneComment("TotalIonCurrent", ticss.str());
    }

    if (precIntCount > 0) {
        stringstream precIntss;
        precIntss.precision(2);
        precIntss << sumPrecInt / (double) (precIntCount);
        final->setOneComment("PrecursorIntensity", precIntss.str());
    }

    if (origMaxIntCount > 0) {
        stringstream origMaxIntss;
        origMaxIntss.precision(2);
        origMaxIntss << sumOrigMaxInt / (double) (origMaxIntCount);
        final->setOneComment("OrigMaxIntensity", origMaxIntss.str());
    }


    map<int, int> massdiffMap;
    for (vector<Replicate *>::iterator r = usedReps.begin(); r != usedReps.end(); r++) {
        int massdiffInt = (*r)->entry->getMassDiffInt();
        if (massdiffMap.find(massdiffInt) == massdiffMap.end()) {
            massdiffMap[massdiffInt] = 1;
        } else {
            massdiffMap[massdiffInt]++;
        }
    }
    if (!massdiffMap.empty()) {
        stringstream mdcss;
        mdcss << massdiffMap.size() << '/';
        for (map<int, int>::iterator md = massdiffMap.begin(); md != massdiffMap.end(); md++) {
            if (md != massdiffMap.begin()) mdcss << ',';
            mdcss << md->first << ':' << md->second;
        }
        final->setOneComment("MassDiffCounts", mdcss.str());
    }

    // record probability range
    // ProbRange lists 4 possible aggregate probabilites:
    // adjProb = 1 - (1-P1)(1-P2)...(1-Pn)  i.e. 1 - chance of all of replicates being false
    // maxProb - max(P1, P2,...Pn), which I think is the most reasonable estimate
    // aveProb = ave(P1, P2,...Pn)
    // minProb = min(P1, P2,...Pn)
    aveProb /= (double) (usedReps.size());
    adjProb = 1 - adjProb;
    stringstream probRangess;
    probRangess << adjProb << ',' << maxProb << ',' << aveProb << ',' << minProb;
    final->setOneComment("ProbRange", probRangess.str());

    // record maximum replicate S/N
    stringstream maxsnss;
    maxsnss.precision(1);
    maxsnss << fixed << showpoint << maxSN;
    final->setOneComment("MaxRepSN", maxsnss.str());

    // CONSENSUS-ONLY AGGREGATIONS

    if (isConsensus) {

        // prob is taken as the maximum of the probabilities of the replicates
        stringstream probss;
        probss.precision(4);
        probss << fixed << maxProb;
        final->setOneComment("Prob", probss.str());


        // The consensus's S/N ratio is taken to be the greater of 200 and the maximum S/N among all replicates
        // Typically, since the consensus creation gets rid of most of the noise peaks, the S/N calculation doesn't work
        // any more for consensus spectra. We just assume the consensus must have very good S/N, which they almost
        // always do, at least at good as 200.0 and up to the best S/N of all the replicates.
        stringstream snss;
        snss.precision(1);
        if (maxSN < 200.0) {
            snss << fixed << showpoint << 200.0;
        } else {
            snss << fixed << showpoint << maxSN;
        }
        final->setOneComment("SN", snss.str());

        // The consensus's xcorr is taken to be the maximum xcorr of all the replicates
        if (!m_missingXCorr) {
            stringstream xcorrss;
            xcorrss.precision(4);
            xcorrss << fixed << showpoint << maxXCorr;
            final->setOneComment("XCorr", xcorrss.str());
        } else {
            final->deleteOneComment("XCorr");
        }
    }


}

// plotSingle - plots a single spectrum
void SpectraSTReplicates::plotSingle(SpectraSTLibEntry *single) {

    stringstream fnss;

    fnss << m_plotPath << single->getSafeName() << "_SING";

    single->getPeakList()->plot(fnss.str(), "SINGLE");

}

// plotConsensus - plots a consensus spectrum
void SpectraSTReplicates::plotConsensus(SpectraSTLibEntry *consensus, SpectraSTPeakList *consensusPeakList,
                                        vector<Replicate *> &usedReps) {

    stringstream fnss;

    fnss << m_plotPath << consensus->getSafeName() << "_CONS";
    consensusPeakList->plot(fnss.str(), "CONSENSUS");

    unsigned int icount = 0;
    for (vector<Replicate *>::iterator rep = usedReps.begin(); rep != usedReps.end(); rep++) {
        stringstream rfnss;
        icount++;
        if (icount == 1) {
            rfnss << m_plotPath << consensus->getSafeName() << "_" << "BEST";
        } else {
            rfnss << m_plotPath << consensus->getSafeName() << "_" << icount;
        }
        stringstream rss;
        rss.precision(3);
        rss << (*rep)->entry->getPeakList()->getWeight();
        (*rep)->entry->getPeakList()->plot(rfnss.str(), rss.str());
    }


}

// sortReplicatesByWeight - comparator for sorting replicates by weight.
bool SpectraSTReplicates::sortReplicatesByWeight(Replicate a, Replicate b) {
    return (a.entry->getPeakList()->getWeight() > b.entry->getPeakList()->getWeight());
}

double SpectraSTReplicates::getMedian(vector<double> &v) {

    if (v.size() % 2 == 0) {
        vector<double>::iterator mid1 = v.begin() + v.size() / 2 - 1;
        vector<double>::iterator mid2 = v.begin() + v.size() / 2;

        std::nth_element(v.begin(), mid1, v.end());
        std::nth_element(v.begin(), mid2, v.end());

        return ((*mid1 + *mid2) / 2.0);

    }

    vector<double>::iterator mid = v.begin() + v.size() / 2;
    std::nth_element(v.begin(), mid, v.end());

    return (*mid);
}
