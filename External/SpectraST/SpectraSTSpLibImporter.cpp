#include "SpectraSTSpLibImporter.hpp"
#include "SpectraSTReplicates.hpp"
#include "SpectraSTFastaFileHandler.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>
#include <set>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

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

/* Class: SpectraSTSpLibImporter
 * 
 * Implements a library importer for the .splib file format. This takes a SpectraST processed
 * .splib file (which is already searchable), and performs certain actions to the library,
 * and writes the resulting new library to a new .splib file (with the accompanying indices).
 * 
 * These actions include:
 * - Uniquify spectra by (a) taking the best replicate; (b) taking consensus of all replicates
 * - Union (based on peptide) multiple libraries
 * - Intersect (based on peptide) multiple libraries
 * - Subtract one library from another (based on peptide)
 * - Filtering based on criteria
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTSpLibImporter::SpectraSTSpLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                               SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params),
        m_splibFins(),
        m_pepIndices(),
        m_mzIndices(),
        m_plotPath(""),
        m_QFSearchParams(NULL),
        m_QFSearchLib(NULL),
        m_ppMappings(NULL),
        m_denoisers(NULL),
        m_singletonPeptideIons() {

    if (params.outputFileName.empty()) {
        // the constructor of the parent class SpectraSTLibImporter only creates the "default" outputFileName
        // need to overwrite this using our own constructOutputFileName
        m_outputFileName = constructOutputFileName();
    }

    // if plotting is required, make a directory for it
    FileName fn;
    parseFileName(m_outputFileName, fn);
    m_plotPath = fn.path + fn.name + "_spplot/";
    if (!m_params.plotSpectra.empty()) {
        makeDir(m_plotPath);
    }

    if (params.useBayesianDenoiser) {
        m_denoisers = new vector<SpectraSTDenoiser *>(MAX_CHARGE + 1, (SpectraSTDenoiser *) NULL);
        for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
            (*m_denoisers)[charge] = new SpectraSTDenoiser();
        }
        if (!(params.trainBayesianDenoiser)) {
            // unless we are training the denoiser on the fly (from consensus building), just use defaults
            if (!(params.denoiserParamFile.empty())) {
                readDenoiserParamFile(params.denoiserParamFile);

                string dpc("denoise.params.copy");
                writeDenoiserParamFile(dpc);

            } else {

                for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
                    (*m_denoisers)[charge]->useDefault();
                }
            }
        } else {

            for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
                (*m_denoisers)[charge]->initializeTraining();
            }
        }
    }

}

// destructor
SpectraSTSpLibImporter::~SpectraSTSpLibImporter() {

    // closes the input files
    for (vector<ifstream *>::iterator i = m_splibFins.begin(); i != m_splibFins.end(); i++) {
        if (*i) delete (*i);
    }
    // deletes the peptide indices loaded into memory
    for (vector<SpectraSTPeptideLibIndex *>::iterator j = m_pepIndices.begin(); j != m_pepIndices.end(); j++) {
        if (*j) delete (*j);
    }

    // deletes the mz indices loaded into memory
    for (vector<SpectraSTMzLibIndex *>::iterator k = m_mzIndices.begin(); k != m_mzIndices.end(); k++) {
        if (*k) delete (*k);
    }

    if (m_ppMappings) {
        for (map<string, vector<pair<string, string> > *>::iterator m = m_ppMappings->begin();
             m != m_ppMappings->end(); m++) {
            if (m->second) delete (m->second);
        }
        delete (m_ppMappings);
    }

    // deletes the library and search params object used for quality filter - conflicting ID
    if (m_QFSearchParams) delete m_QFSearchParams;
    if (m_QFSearchLib) delete m_QFSearchLib;


    if (m_denoisers) {
        for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
            if ((*m_denoisers)[charge]) delete ((*m_denoisers)[charge]);
        }
        delete (m_denoisers);
    }
}


// import - reads entries from .splib files
void SpectraSTSpLibImporter::import() {
    // add by long
    cout << "[Info] inside importer" << endl;

    // check conflict of SUBTRACT_HOMOLOGS and build actions
    if (m_params.combineAction == "SUBTRACT_HOMOLOGS" && (!m_params.buildAction.empty())) {
        g_log->error("CREATE", "Cannot perform build action " + m_params.buildAction +
                               " together with combine action SUBTRACT_HOMOLOGS.");
        return;
    }
    // add by long
    cout << "[Info] inside importer 1" << endl;

    if (m_params.buildAction == "QUALITY_FILTER") {
        // quality filter is very different, so we'll go to another method to do it
        doQualityFilter();
        return;
    }

    cout << "[Info] inside importer 2" << endl;

    if (m_params.buildAction == "DECOY") {
        doGenerateDecoy();
        return;
    }
    cout << "[Info] inside importer 3" << endl;

    if (m_params.buildAction == "SORT_BY_NREPS") {
        doSortByNreps();
        return;
    }

    cout << "[Info] inside importer 4" << endl;

    if (m_params.buildAction == "USER_SPECIFIED_MODS") {
        doUserSpecifiedModifications();
        return;
    }

    cout << "[Info] inside importer 5" << endl;

    if (m_params.buildAction == "SIMILARITY_CLUSTERING") {
        doSimilarityClustering();
        return;
    }

    cout << "[Info] inside importer 6" << endl;

    if (m_params.buildAction == "SEMI_EMPIRICAL_SPLIB") {
        doBuildSemiEmpiricalSplib();
        return;
    }

    cout << "[Info] inside importer 7" << endl;

    string fileListStr = constructFileListStr();

    // print starting message to console
    if (!g_quiet) {
        if (m_params.buildAction == "BEST_REPLICATE") {
            cout << "Creating BEST REPLICATE library from " << fileListStr << endl;
        } else if (m_params.buildAction == "CONSENSUS") {
            cout << "Creating CONSENSUS library from " << fileListStr << endl;
        } else {
            cout << "Creating library from " << fileListStr << endl;
        }
    }

    cout << "[Info] inside importer 8" << endl;

    string desc = m_params.constructDescrStr(fileListStr, ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);
    cout << "[Info] inside importer 8.1" << endl;

    // for SUBTRACT_HOMOLOGS, need a slightly different order of retrieving peptides, so branch out to another method here
    if (m_params.combineAction == "SUBTRACT_HOMOLOGS") {
        doSubtractHomologs();
        return;
    }

    cout << "[Info] inside importer 9" << endl;

    openSplibs(false, 0.0, true, false, true);

    // write the preamble of the generated .splib file
    m_lib->writePreamble(m_preamble);

    if (!(m_params.refreshDatabase.empty())) {
        refresh();
    }

    // add by long
    cout << "[Info] start reading file" << endl;

    // now go through the peptide indices one by one and load the entries
    // the procedure is a bit counterintuitive, as follows: open the first file, for each peptide ion, add
    // the entries from the first file, then go look for that peptide ion in the rest of the files,
    // add them as well. then open the second file, for those peptide ions that are not already
    // added in the first pass, do the same. then open the third file... etc...

    string peptide;
    string mods;

    vector<string> subkeys;

    vector<SpectraSTPeptideLibIndex *>::size_type curPepIndex = 0;
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[curPepIndex];
    if (!pepIndex) return; // require the first file to be okay

    ProgressCount pc(!g_quiet && !g_verbose, 500, 0);
    pc.start("Importing ions");

    while (true) {
        if (curPepIndex >= m_pepIndices.size()) {
            break;
        }

        while ((pepIndex = m_pepIndices[curPepIndex]) &&
               (pepIndex->nextPeptide(peptide, subkeys))) {

            for (vector<string>::iterator k = subkeys.begin(); k != subkeys.end(); k++) {

                bool include = true;

                if (m_params.combineAction == "UNION" || m_params.combineAction == "APPEND") {
                    if (curPepIndex != 0 && (m_lib->getPeptideLibIndexPtr()->isInIndex(peptide, *k))) {
                        // already included this peptide ion, don't do it again
                        include = false;
                    }
                } else if (m_params.combineAction == "INTERSECT") {
                    // only include if this peptide ion is in every one of the file!
                    for (vector<SpectraSTPeptideLibIndex *>::iterator i = m_pepIndices.begin() + 1;
                         i != m_pepIndices.end(); i++) {
                        if ((*i) && !((*i)->isInIndex(peptide, *k))) {
                            include = false;
                            break;
                        }
                    }
                } else if (m_params.combineAction == "SUBTRACT") {
                    // only include if this peptide ion is in the first file but not in any of the rest
                    for (vector<SpectraSTPeptideLibIndex *>::iterator i = m_pepIndices.begin() + 1;
                         i != m_pepIndices.end(); i++) {
                        if ((*i) && (*i)->isInIndex(peptide, *k)) {
                            include = false;
                            break;
                        }
                    }
                }

                // now actually retrieve the entries
                vector<SpectraSTLibEntry *> entries;
                if (include) {

                    if (g_verbose) {
                        cout << "Importing ion " << peptide << '/' << (*k) << " ... ";
                    }

                    for (vector<SpectraSTPeptideLibIndex *>::iterator i = m_pepIndices.begin();
                         i != m_pepIndices.end(); i++) {
                        if (*i) {
                            (*i)->retrieve(entries, peptide, *k);
                            if (m_params.combineAction == "APPEND" && (!(entries.empty()))) {
                                break;
                            }
                        }
                        if (m_params.combineAction == "SUBTRACT") {
                            // by definition of subtraction, no need to look into other files besides the first one
                            break;
                        }
                    }

                } else {
                    continue;
                }

                m_count++;
                pc.increment();

                if (g_verbose) {
                    cout << " (" << entries.size() << " replicates)" << endl;
                }

                // perform the build actions
                doBuildAction(entries);

                // done, delete all entries
                for (vector<SpectraSTLibEntry *>::iterator den = entries.begin(); den != entries.end(); den++) {
                    delete (*den);
                }


            } // for (subkeys)

        } // while (nextPeptide)

        //   cerr << "Ave B/Y Assigned Ratio = " << m_totBYRatio / (float)(m_lib->getCount()) << endl;


        if (m_params.combineAction != "UNION" && m_params.combineAction != "APPEND") {
            // if we are not doing a UNION, one pass is enough.
            break;
        }

        // onto the next file -- necessary for UNION
        curPepIndex++;

    } // while true (looping over all peptide indices)

    pc.done();

    if (m_denoisers && m_params.trainBayesianDenoiser) {

        for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
            (*m_denoisers)[charge]->generateBayesianModel();
        }

        if (!(m_params.denoiserParamFile.empty())) {
            writeDenoiserParamFile(m_params.denoiserParamFile);
        }

        if (!(m_singletonPeptideIons.empty())) {
            // re-read the singletons and denoise them before writing to library
            reloadAndProcessSingletons();
        }

    }

}

void SpectraSTSpLibImporter::reloadAndProcessSingletons() {

    ProgressCount pc(!g_quiet, 1, (unsigned int) (m_singletonPeptideIons.size()));
    pc.start("De-noise singleton raw spectra using trained Bayesian signal/noise classifier");

    for (vector<pair<string, string> >::iterator si = m_singletonPeptideIons.begin();
         si != m_singletonPeptideIons.end(); si++) {
        vector<SpectraSTLibEntry *> entries;
        for (vector<SpectraSTPeptideLibIndex *>::iterator i = m_pepIndices.begin(); i != m_pepIndices.end(); i++) {
            if (*i) {
                (*i)->retrieve(entries, si->first, si->second);
                if (m_params.combineAction == "APPEND" && (!(entries.empty()))) {
                    break;
                }
            }
            if (m_params.combineAction == "SUBTRACT") {
                // by definition of subtraction, no need to look into other files besides the first one
                break;
            }
        }

        // perform the build actions -- although they are singletons,
        // still have to do it to deal with removeDissimilar issue and activate de-noiser
        doBuildAction(entries);

        // done, delete all entries
        for (vector<SpectraSTLibEntry *>::iterator den = entries.begin(); den != entries.end(); den++) {
            delete (*den);
        }
    }

    pc.done();

}

// doSubtractHomologs - performs the special join action SUBTRACT_HOMOLOGS
// this will include only those entries in the first .splib files that don't have a homolog in any of the other .splib files
void SpectraSTSpLibImporter::doSubtractHomologs() {

    openSplibs(true, 13.0, false, false, true);

    m_lib->writePreamble(m_preamble);

    string peptide;
    string mods;

    ProgressCount pc(!g_quiet && !g_verbose, 500, 0);
    pc.start("Importing peptide ions");

    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];
    if (!mzIndex) return;

    // loop through each entries in the m/z index of the first .splib file
    SpectraSTLibEntry *entry = NULL;
    while ((entry = mzIndex->nextEntry())) {

        Peptide *pep = entry->getPeptidePtr();

        if (!pep) continue;

        double mz = entry->getPrecursorMz();
        int charge = entry->getCharge();
        bool include = true;

        for (vector<SpectraSTMzLibIndex *>::iterator i = m_mzIndices.begin() + 1; i != m_mzIndices.end(); i++) {

            if (!(*i)) continue;

            vector<SpectraSTLibEntry *> isobaricEntries;

            // retrieve all entries within 4.5 Th of this peptide ion from all the other .splib files
            (*i)->retrieve(isobaricEntries, mz - 4.5, mz + 4.5);

            for (vector<SpectraSTLibEntry *>::iterator en = isobaricEntries.begin();
                 en != isobaricEntries.end(); en++) {

                int identity = 0;
                Peptide *thisPep = (*en)->getPeptidePtr();
                if (!thisPep) {
                    continue;
                }

                if (*pep == *thisPep || (charge == (*en)->getCharge() && pep->isHomolog(*thisPep, 0.7, identity))) {
                    stringstream logss;
                    logss << pep->interactStyleWithCharge() << " (m/z = " << mz << ") is homologous (" << identity
                          << ") to ";
                    logss << thisPep->interactStyleWithCharge() << " (m/z = " << (*en)->getPrecursorMz()
                          << "). Removed.";
                    g_log->log("CREATE", logss.str());
                    include = false;
                    break;
                }
            }

            if (!include) {
                break;
            }
        }

        if (include && insertOneEntry(entry, "CREATE")) {

            m_count++;
            pc.increment();

            if (g_verbose) {
                cout << "Importing peptide ion: " << pep->interactStyleWithCharge() << endl;
            }
        }

        delete (entry);
    }

    pc.done();


}

// doBuildAction - performs the build actions BEST_REPLICATE and CONSENSUS
void SpectraSTSpLibImporter::doBuildAction(vector<SpectraSTLibEntry *> &entries) {

    if (m_params.buildAction == "BEST_REPLICATE") {

        SpectraSTReplicates *replicates = new SpectraSTReplicates(entries, m_params);
        replicates->setPlotPath(m_plotPath);

        // pick the best replicate and only insert that one into the library
        SpectraSTLibEntry *best = replicates->findBestReplicate();
        if (best) {
            //  best->getPeakList()->rankTransformWithQuota(50, true, 180.0);
            insertOneEntry(best, "BEST_REPLICATE");
        }
        delete (replicates);

    } else if (m_params.buildAction == "CONSENSUS") {
        SpectraSTReplicates *replicates = new SpectraSTReplicates(entries, m_params, m_denoisers);
        replicates->setPlotPath(m_plotPath);
        // make a consensus spectrum of replicates and insert that into the library
        SpectraSTLibEntry *consensus = replicates->makeConsensusSpectrum();

        if (consensus) {
            if (m_denoisers && !((*m_denoisers)[0]->isFilterReady()) && consensus->getNrepsUsed() == 1) {
                // just remember this peptide ion, don't write to library yet
                pair<string, string> p;
                p.first = consensus->getPeptidePtr()->stripped;
                p.second = SpectraSTPeptideLibIndex::constructSubkey(consensus);
                m_singletonPeptideIons.push_back(p);
            } else {
                insertOneEntry(consensus, "CONSENSUS");
            }
        }
        delete (replicates);

    } else {

        // no build action. Just add everything back (except those filtered)
        for (vector<SpectraSTLibEntry *>::iterator ien = entries.begin(); ien != entries.end(); ien++) {
            insertOneEntry(*ien, "CREATE");
        }

    }


}


// doQualityFilter - performs the quality filters on all library entries
void SpectraSTSpLibImporter::doQualityFilter() {

    QFStats qfstats;
    memset(&qfstats, '\0', sizeof(qfstats));

    if (m_impFileNames.size() != 1) {
        g_log->error("QUALITY FILTER",
                     "Quality filter must be applied to one .splib file only. No filtering performed.");
        return;
    }

    if (!g_quiet) {
        cout << "Applying QUALITY FILTER to \"" << m_impFileNames[0] << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(true, 13.0, true, true, true);
    ifstream *splibFin = m_splibFins[0];
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[0];

    if (splibFin && !pepIndex) {
        // this means checkUniqueness in openSplibs failed.
        // some peptide ions are non-unique
        g_log->error("QUALITY_FILTER", "Quality filter requires unique library. Library \"" + m_impFileNames[0] +
                                       "\" is non-unique. No filtering is performed.");
        return;
    }

    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];
    if (!mzIndex) return;

    m_lib->writePreamble(m_preamble);

    // if conflicting ID's are to be detected, each library entry is searched against the same library
    // need to open the same library as a different object for searching
    if (m_params.qualityLevelMark >= 2 || m_params.qualityLevelRemove >= 2) {
        m_QFSearchParams = new SpectraSTSearchParams();
        m_QFSearchLib = new SpectraSTLib(m_impFileNames[0], m_QFSearchParams, true);
    }

    ProgressCount pc(!g_quiet && !g_verbose, 500, 0);
    pc.start("Importing peptide ions");

    SpectraSTLibEntry *entry = NULL;

    // loop through the library by precursor m/z. entry will get NULL if the end is reached.
    while ((entry = mzIndex->nextEntry())) {

        pc.increment();

        if (applyQualityFilter(entry, qfstats)) {
            insertOneEntry(entry, "QUALITY_FILTER");
        }

        delete (entry);
        entry = NULL;
    }

    pc.done();

    // if we actually go through all 5 levels, and don't remove any, then we have some good statistics
    // output them to the log file for the user to examine
    if (m_params.qualityLevelMark >= 5 && m_params.qualityLevelRemove == 0) {

        stringstream qfstatss;
        qfstatss << "immune_prob = " << qfstats.immuneProb;
        qfstatss << "; immune_engine = " << qfstats.immuneEngine << endl;

        g_log->log("QUALITY_FILTER STATS", qfstatss.str());

        // these calculates the number of spectra that will remain at each quality level
        unsigned int level1 = m_count - qfstats.Q1;
        unsigned int level2 = level1 - qfstats.Q2 + qfstats.Q1Q2;
        unsigned int level3 = level2 - qfstats.Q3 + qfstats.Q1Q3 + qfstats.Q2Q3 - qfstats.Q1Q2Q3;
        unsigned int level4 =
                level3 - qfstats.Q4 + qfstats.Q1Q4 + qfstats.Q2Q4 + qfstats.Q3Q4 - qfstats.Q1Q2Q4 - qfstats.Q1Q3Q4 -
                qfstats.Q2Q3Q4 + qfstats.Q1Q2Q3Q4;
        unsigned int level5 =
                level4 - qfstats.Q5 + qfstats.Q1Q5 + qfstats.Q2Q5 + qfstats.Q3Q5 + qfstats.Q4Q5 - qfstats.Q1Q2Q5 -
                qfstats.Q1Q3Q5 - qfstats.Q1Q4Q5 - qfstats.Q2Q3Q5 - qfstats.Q2Q4Q5 - qfstats.Q3Q4Q5 + qfstats.Q1Q2Q3Q5 +
                qfstats.Q1Q2Q4Q5 + qfstats.Q1Q3Q4Q5 + qfstats.Q2Q3Q4Q5 - qfstats.Q1Q2Q3Q4Q5;

        stringstream qfstatss2;
        qfstatss2 << "Level 0 = " << m_count;
        qfstatss2 << "; Level 1 = " << level1;
        qfstatss2 << "; Level 2 = " << level2;
        qfstatss2 << "; Level 3 = " << level3;
        qfstatss2 << "; Level 4 = " << level4;
        qfstatss2 << "; Level 5 = " << level5;
        g_log->log("QUALITY_FILTER STATS", qfstatss2.str());
    }


}


// applyQualityFilter - applies the various quality filters, sets the status and comments accordingly, and
// returns false if the entry is to be deleted
bool SpectraSTSpLibImporter::applyQualityFilter(SpectraSTLibEntry *entry, QFStats &qfstats) {

    bool inquorate = false;
    bool singleton = false;
    bool impure = false;
    bool unconfirmed = false;
    bool badConflictingID = false;

    // reset status
    entry->setStatus("Normal");
    m_count++;

    double prob = entry->getProb();

    if (prob >= m_params.qualityImmuneProbThreshold) {
        // immune by probability
        qfstats.immuneProb++;
        return (true);
    }

    unsigned int numSeqEngines = 1;
    if (m_params.qualityImmuneMultipleEngines) {
        string seqStr("");
        string::size_type dummy = 0;
        if (entry->getOneComment("Se", seqStr)) {
            numSeqEngines = atoi(nextToken(seqStr, 0, dummy, "^/", "").c_str());
            if (numSeqEngines > 1) {
                // immune by multiple search engines
                qfstats.immuneEngine++;
                return (true);
            }
        }
    }

    // check minimum numbers of replicates
    unsigned int numRepsUsed = entry->getNrepsUsed();


    // create a string to describe this entry
    stringstream tagss;
    tagss << "entry #" << entry->getLibId() << " : " << entry->getName();
    tagss << " (" << numRepsUsed << " replicates; ";
    if (numSeqEngines > 1) tagss << numSeqEngines << " engines; ";
    tagss << "P=" << prob << ") ";
    string tag = tagss.str();

    // in quality filter, the repQuorum must be at least 2 (all singletons are inquorate)
    // if it's set at 1 then the logic will fail
    unsigned int repQuorum = m_params.minimumNumReplicates;
    if (repQuorum < 2) repQuorum = 2;

    if (numRepsUsed < repQuorum) {
        inquorate = true;
    }

    // level 5 is inquorate entries (i.e. numRepsUsed < repQuorum)
    if (m_params.qualityLevelRemove >= 5 || m_params.qualityLevelMark >= 5) {
        if (inquorate) {
            qfstats.Q5++;

            if (m_params.qualityLevelRemove >= 5) {
                g_log->log("QUALITY_FILTER", "Remove INQUORATE " + tag);
                return (false);
            }
            if (m_params.qualityLevelMark >= 5) {
                g_log->log("QUALITY_FILTER", "Mark INQUORATE " + tag);
                entry->setStatus("Inquorate");
            }
        }
    }

    // level 4 is singleton entries (those with numUsedRep == 1)
    // by definition all singleton entries must be inquorate too
    if (inquorate && (m_params.qualityLevelRemove >= 4 || m_params.qualityLevelMark >= 4)) {
        if (numRepsUsed == 1) {
            singleton = true;
            qfstats.Q4++;
            qfstats.Q4Q5++;
            if (m_params.qualityLevelRemove >= 4) {
                g_log->log("QUALITY_FILTER", "Remove SINGLETON " + tag);
                return (false);
            }
            if (m_params.qualityLevelMark >= 4) {
                g_log->log("QUALITY_FILTER", "Mark SINGLETON " + tag);
                entry->setStatus("Singleton");

            }
        }
    }

    // level 3 is inquorate_unconfirmed entries (those that don't have any other peptide ion in the library with a shared sequence)
    // only apply this to inquorate entries
    if (inquorate && (m_params.qualityLevelRemove >= 3 || m_params.qualityLevelMark >= 3)) {
        if (!hasSharedSequence(entry)) {
            unconfirmed = true;
            qfstats.Q3++;
            qfstats.Q3Q5++;
            if (singleton) {
                qfstats.Q3Q4++;
                qfstats.Q3Q4Q5++;
            }
            stringstream ss;
            if (m_params.qualityLevelRemove >= 3) {
                g_log->log("QUALITY_FILTER", "Remove INQUORATE_UNCONFIRMED " + tag);
                return (false);
            }
            if (m_params.qualityLevelMark >= 3) {
                g_log->log("QUALITY_FILTER", "Mark INQUORATE_UNCONFIRMED " + tag);
                entry->setStatus("Inquorate_Unconfirmed");
            }
        }
    }

    // level 2 is conflicting IDs (those that have a spectrally similar counterpart in the library with a different ID.)
    if (m_params.qualityLevelRemove >= 2 || m_params.qualityLevelMark >= 2) {

        string msg("");
        if (isBadConflictingID(entry, numRepsUsed, m_params.qualityPenalizeSingletons, msg)) {
            badConflictingID = true;
            qfstats.Q2++;
            if (inquorate) qfstats.Q2Q5++;
            if (singleton) {
                qfstats.Q2Q4++;
                qfstats.Q2Q4Q5++;
            }
            if (unconfirmed) {
                qfstats.Q2Q3++;
                qfstats.Q2Q3Q5++;
            }
            if (singleton && unconfirmed) {
                qfstats.Q2Q3Q4++;
                qfstats.Q2Q3Q4Q5++;
            }
            if (m_params.qualityLevelRemove >= 2) {
                g_log->log("QUALITY_FILTER", "Remove CONFLICTING_ID " + tag + msg);
                return (false);
            }
            if (m_params.qualityLevelMark >= 2) {
                g_log->log("QUALITY_FILTER", "Mark CONFLICTING_ID " + tag + msg);
                entry->setStatus("Conflicting_ID");
            }
        } else if (!msg.empty()) {
            g_log->log("QUALITY_FILTER", "Keep CONFLICTING_ID " + tag + msg);
        }

    }

    // level 1 is impure entries (those entries with impure spectra)
    if (m_params.qualityLevelRemove >= 1 || m_params.qualityLevelMark >= 1) {
        string msg("");
        if (isImpure(entry, numRepsUsed, m_params.qualityPenalizeSingletons, msg)) {
            impure = true;
            qfstats.Q1++;
            if (inquorate) qfstats.Q1Q5++;
            if (singleton) {
                qfstats.Q1Q4++;
                qfstats.Q1Q4Q5++;
            }
            if (unconfirmed) {
                qfstats.Q1Q3++;
                qfstats.Q1Q3Q5++;
            }
            if (singleton && unconfirmed) {
                qfstats.Q1Q3Q4++;
                qfstats.Q1Q3Q4Q5++;
            }
            if (badConflictingID) qfstats.Q1Q2++;
            if (badConflictingID && inquorate) qfstats.Q1Q2Q5++;
            if (badConflictingID && singleton) {
                qfstats.Q1Q2Q4++;
                qfstats.Q1Q2Q4Q5++;
            }
            if (badConflictingID && unconfirmed) {
                qfstats.Q1Q2Q3++;
                qfstats.Q1Q2Q3Q5++;
            }
            if (badConflictingID && singleton && unconfirmed) {
                qfstats.Q1Q2Q3Q4++;
                qfstats.Q1Q2Q3Q4Q5++;
            }


            if (m_params.qualityLevelRemove >= 1) {
                g_log->log("QUALITY_FILTER", "Remove IMPURE " + tag + msg);
                return (false);
            }
            if (m_params.qualityLevelMark >= 1) {
                g_log->log("QUALITY_FILTER", "Mark IMPURE " + tag + msg);
                entry->setStatus("Impure");
            }
        }
    }

    // enforce peak quorum
    unsigned int minNumRepWithPeak = (unsigned int) (numRepsUsed * m_params.peakQuorum - 0.00001) + 1;
    if (minNumRepWithPeak < 1) minNumRepWithPeak = 1;
    entry->getPeakList()->removeInquoratePeaks(minNumRepWithPeak);

    return (true);

}

// isImpure - checks the purity of the spectrum
bool SpectraSTSpLibImporter::isImpure(SpectraSTLibEntry *entry, unsigned int numRepsUsed, bool penalizeSingletons,
                                      string &msg) {

    if (!entry->getPeptidePtr()) {
        return (false);
    }

    entry->annotatePeaks();

    unsigned int numUnassignedAll = 0;
    unsigned int numAssignedAll = 0;
    unsigned int numUnassignedTop20 = 0;
    unsigned int numAssignedTop20 = 0;
    unsigned int numUnassignedTop5 = 0;
    unsigned int numAssignedTop5 = 0;
    double fracUnassignedAll = 0.0;
    double fracUnassignedTop20 = 0.0;
    double fracUnassignedTop5 = 0.0;

    // check fraction unassigned
    string fracUnassignedStr("");
    if (entry->getOneComment("FracUnassigned", fracUnassignedStr)) {
        // already calculated, just parse out the information

        string::size_type semicolonPos = 0;
        string::size_type commaPos = 0;
        string top5Str = nextToken(fracUnassignedStr, semicolonPos, semicolonPos, ";\r\t\n");
        string top20Str = nextToken(fracUnassignedStr, semicolonPos + 1, semicolonPos, ";\r\t\n");
        string allStr = nextToken(fracUnassignedStr, semicolonPos + 1, semicolonPos, ";\r\t\n");

        fracUnassignedTop20 = atof(nextToken(top20Str, 0, commaPos, ",\t\r\n").c_str());
        numUnassignedTop20 = atoi(nextToken(top20Str, commaPos + 1, commaPos, "/\t\r\n").c_str());
        numAssignedTop20 = atoi(nextToken(top20Str, commaPos + 1, commaPos, "/\t\r\n").c_str()) - numUnassignedTop20;

        fracUnassignedTop5 = atof(nextToken(top5Str, 0, commaPos, ",\t\r\n").c_str());
        numUnassignedTop5 = atoi(nextToken(top5Str, commaPos + 1, commaPos, "/\t\r\n").c_str());
        numAssignedTop5 = atoi(nextToken(top5Str, commaPos + 1, commaPos, "/\t\r\n").c_str()) - numUnassignedTop5;


        // don't need the other unassigned info, so don't parse

    } else {
        // not calculated, need to calculate anew

        unsigned int NAA = entry->getPeptidePtr()->NAA();

        fracUnassignedAll = entry->getPeakList()->calcFractionUnassigned(999999, numUnassignedAll, numAssignedAll);
        fracUnassignedTop20 = entry->getPeakList()->calcFractionUnassigned(20, numUnassignedTop20, numAssignedTop20,
                                                                           true, true);
        fracUnassignedTop5 = entry->getPeakList()->calcFractionUnassigned(5, numUnassignedTop5, numAssignedTop5, true,
                                                                          true);

        // as a bonus, stick this information into the Comment for future use
        stringstream fracUnassignedss;
        fracUnassignedss.precision(2);
        fracUnassignedss << fixed << showpoint << fracUnassignedTop5 << ',' << numUnassignedTop5 << '/'
                         << numUnassignedTop5 + numAssignedTop5 << ';';
        fracUnassignedss << fixed << showpoint << fracUnassignedTop20 << ',' << numUnassignedTop20 << '/'
                         << numUnassignedTop20 + numAssignedTop20 << ';';
        fracUnassignedss << fixed << showpoint << fracUnassignedAll << ',' << numUnassignedAll << '/'
                         << numUnassignedAll + numAssignedAll;
        entry->setOneComment("FracUnassigned", fracUnassignedss.str());

    }

    // logic: charge +1 spectra are immune (they often appear to be impure -- many unassigned peaks -- anyway)
    // otherwise, impure if fracUnassigned of the top 20 peaks is above 0.4.
    // if penalizeSingletons is TRUE, then singleton spectra are considered impure if the fracUnassigned of the top 20 peaks
    // is above 0.4 OR the number of unassigned peaks in the top 20 is more than 9.
    if (entry->getCharge() == 1 ||
        ((!penalizeSingletons || numRepsUsed > 1) && fracUnassignedTop20 < 0.4) ||
        (penalizeSingletons && numRepsUsed == 1 &&
         (fracUnassignedTop20 < 0.4 && numUnassignedTop20 < numAssignedTop20 - 2))) {

        msg = "";
        return (false);

    } else {
        stringstream ss;
        ss << "| FRAC UNASSIGNED " << fracUnassignedTop20 << ";" << numUnassignedTop20 << '/'
           << numUnassignedTop20 + numAssignedTop20;
        msg = ss.str();
        return (true);
    }
}

// hasSharedSequence - checks if there is another peptide ion with a shared sequence with this entry. If so,
// returns true. (This means this entry is less likely to be a false positive)
bool SpectraSTSpLibImporter::hasSharedSequence(SpectraSTLibEntry *entry) {

    if (!entry->getPeptidePtr()) {
        return (false);
    }

    string foundPeptide("");
    string foundSubkey("");
    if (m_QFSearchLib->getPeptideLibIndexPtr()->hasSharedSequence(*(entry->getPeptidePtr()), foundPeptide,
                                                                  entry->getFragType())) {
        return (true);
    } else {
        return (false);
    }

}

// isBadConflictingID - searches entry against the library; if there's a highly similar spectrum at similar precursor m/z,
// returns true if entry is "worse" than the similar-looking library spectrum.
bool
SpectraSTSpLibImporter::isBadConflictingID(SpectraSTLibEntry *entry, unsigned int numRepsUsed, bool penalizeSingletons,
                                           string &msg) {

    if (!entry->getPeptidePtr()) {
        return (false);
    }

    double precursorMz = entry->getPrecursorMz();

    double prob = entry->getProb();

    // retrieves all entries from the library within the tolerable m/z range
    vector<SpectraSTLibEntry *> entries;
    m_QFSearchLib->retrieve(entries, precursorMz - 4.5, precursorMz + 4.5);

    stringstream ss;

//  SpectraSTPeakList* copiedPeakList = new SpectraSTPeakList(*(entry->getPeakList()));

    for (vector<SpectraSTLibEntry *>::iterator i = entries.begin(); i != entries.end(); i++) {
        if (entry->getLibId() != (*i)->getLibId()) {

            double dot = entry->getPeakList()->compare((*i)->getPeakList());

            if (dot >= 0.70 || (penalizeSingletons && numRepsUsed == 1 && dot >= 0.65)) {
                // similar spectra!

                unsigned int matchNumRepsUsed = (*i)->getNrepsUsed();
                double matchProb = (*i)->getProb();

                ss << "| SIMILAR (" << dot << ") to " << (*i)->getLibId() << " : "
                   << (*i)->getPeptidePtr()->interactStyleWithCharge();
                ss << " (" << matchNumRepsUsed << " replicates; P=" << matchProb << ")";

                // check homology, if it is homologous, then we don't apply the conflicting ID filter
                int identity = 0;
                if (entry->getPeptidePtr()->isHomolog(*((*i)->getPeptidePtr()), 0.6, identity)) {
                    // if the spectrum similar to this one belongs to a homologous sequence, then spare it from destruction
                    // Note that this means any search against this library will need the detectHomolog option turned on!
                    ss << " HOMOLOG ";
                    continue;
                }

                // apply filter. always keep the one with more replicates; in case of a tie, the probabilities are the tie-breakers.
                // if even the probs are the same, then keep both.
                if (matchNumRepsUsed > numRepsUsed) {
                    msg = ss.str();
                    return (true);
                } else if (matchNumRepsUsed == numRepsUsed) {
                    if (matchProb > prob) {
                        msg = ss.str();
                        return (true);
                    }
                }
            }
        }
    }


    msg = ss.str();
    return (false);

}

// parsePreamble - parses out the preamble of each of the input .splib files
// then put them back into the final product to leave a trace of what has been done
void SpectraSTSpLibImporter::parsePreamble(ifstream &splibFin, bool binary) {

    if (!binary) {
        char firstChar = (char) (splibFin.peek());
        if (firstChar != '#') return;

        string line;
        string fileName("");
        string firstLine("");
        string versionLine("");
        string::size_type pos = 0;
        while (nextLine(splibFin, line, "### ===", "")) {
            string prefix = line.substr(0, 3);
            if (prefix != "###") {
                return;
            }
            string rest = nextToken(line, 3, pos, "\r\n", " \t");

            if (!rest.empty()) {
                if (fileName.empty()) {
                    fileName = rest;
                } else if (versionLine.empty() && rest.compare(0, 18, "SpectraST (version") == 0) {
                    versionLine = rest;
                } else if (firstLine.empty()) {
                    firstLine = fileName + " : " + rest;
                    m_preamble.push_back("> " + firstLine);
                } else {
                    m_preamble.push_back("> " + rest);
                }
            }
        }

    } else {
        // binary format
        int spectrastVersion = 0;
        int spectrastSubVersion = 0;
        unsigned int numLines = 0;
        string line("");

        splibFin.read((char *) (&spectrastVersion), sizeof(int));
        splibFin.read((char *) (&spectrastSubVersion), sizeof(int));

        if (!nextLine(splibFin, line)) {
            g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
            g_log->crash();
        }
        string fileName(line);

        string firstLine("");

        splibFin.read((char *) (&numLines), sizeof(unsigned int));
        for (unsigned int i = 0; i < numLines; i++) {

            if (!nextLine(splibFin, line)) {
                g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
                g_log->crash();
            }

            if (firstLine.empty()) {
                firstLine = fileName + " : " + line;
                m_preamble.push_back("> " + firstLine);
            } else {
                m_preamble.push_back("> " + line);
            }
        }

    }

}

// constructOutputName - if the user doesn't specify the output file name, try to construct one that makes sense
string SpectraSTSpLibImporter::constructOutputFileName() {

    stringstream ss;

    FileName fn0;
    parseFileName(m_impFileNames[0], fn0);
    ss << fn0.path << fn0.name;

    char oper = 'U';

    if (m_params.combineAction == "UNION") {
        oper = 'U';
    } else if (m_params.combineAction == "INTERSECT") {
        oper = 'I';
    } else if (m_params.combineAction == "SUBTRACT") {
        oper = 'S';
    } else if (m_params.combineAction == "SUBTRACT_HOMOLOG") {
        oper = 'H';
    } else if (m_params.combineAction == "APPEND") {
        oper = 'A';
    }

    if (m_impFileNames.size() < 4) {
        // build a chain name
        for (vector<string>::size_type i = 1; i < m_impFileNames.size(); i++) {
            FileName fn;
            parseFileName(m_impFileNames[i], fn);
            ss << '_' << oper << '_' << fn.name;
        }

    } else {
        // use the first name + operator + "plus"
        ss << '_' << oper << "_plus";
    }

    if (m_params.buildAction == "BEST_REPLICATE") {
        ss << "_best";
    } else if (m_params.buildAction == "CONSENSUS") {
        ss << "_consensus";
    } else if (m_params.buildAction == "QUALITY_FILTER") {
        ss << "_quality";
    } else if (m_params.buildAction == "DECOY") {
        ss << "_decoy";
    } else if (m_params.buildAction == "SORT_BY_NREPS") {
        ss << "_sorted";
    } else if (m_params.buildAction == "USER_SPECIFIED_MODS") {
        ss << "_mods";
    } else if (m_params.buildAction == "SIMILARITY_CLUSTERING") {
        ss << "_clustered";
    } else if (m_params.buildAction == "SEMI_EMPIRICAL_SPLIB") {
        ss << "_semi";
    } else {
        ss << "_new";
    }


    ss << ".splib";

    return (ss.str());
}

// plot - plots an entry
void SpectraSTSpLibImporter::plot(SpectraSTLibEntry *entry) {

    stringstream fnss;

    fnss << m_plotPath << entry->getSafeName();
    entry->getPeakList()->plot(fnss.str(), entry->getStatus());


}

void SpectraSTSpLibImporter::doGenerateDecoy() {

    if (m_params.decoyPrecursorSwap) {
        doGeneratePrecursorSwapDecoy();
    } else {
        doGenerateShuffleAndRepositionDecoy();
    }

}

void SpectraSTSpLibImporter::doGeneratePrecursorSwapDecoy() {

    if (m_impFileNames.size() != 1) {
        g_log->error("DECOY", "Decoy generation must be applied to one .splib file only. No decoy library created.");
        return;
    }

    // don't do reduce in decoy generation - unclear how to reduce both
    // real and corresponding decoys
    m_params.reduceSpectrum = 0;

    if (!g_quiet) {
        cout << "Generating DECOY to \"" << m_impFileNames[0] << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(true, 200.0, true, true, true);

    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[0];
    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];

    if (!pepIndex) {
        // this means checkUniqueness in openSplibs failed.
        // non-unique library, i.e. some peptide ions have multiple spectra.
        g_log->error("DECOY", "Decoy generation requires unique library. Library \"" + m_impFileNames[0] +
                              "\" is non-unique. No decoy is generated.");
        return;
    }

    if (!mzIndex) {
        g_log->error("DECOY", "Cannot open the index of library \"" + m_impFileNames[0] + "\". Aborted.");
        return;
    }

    unsigned int count = mzIndex->getEntryCount();

    m_lib->writePreamble(m_preamble);

    ProgressCount pc(!g_quiet, 1, count);
    pc.start("Generating decoy spectra (precursor swap method)");

    srand(time(NULL));

    set<fstream::off_type> usedEntries;
    vector<SpectraSTLibEntry *> entries;
    fstream::off_type currentOffset;
    double minMz = -1.0;

    while (mzIndex->nextSortedFileOffset(currentOffset)) {

        if (usedEntries.find(currentOffset) != usedEntries.end()) {
            // already used this one, don't read it
            continue;
        }

        SpectraSTLibEntry *entry = mzIndex->thisEntry();

        double precursorMz = entry->getPrecursorMz();
        if (minMz < 0.0) minMz = precursorMz;

        pc.increment();
        pc.increment(); // increment twice -- expect to swap 2 entries at once in most cases

        if (precursorMz > minMz + 3.0) {
            // done with one batch, go find swapping candidates
            // cerr << "BATCH" << endl;
            generatePrecursorSwapDecoys(mzIndex, entries, usedEntries, minMz + 12.0, precursorMz + 12.0);

            // here, the entries list still contain failed cases, will try to handle these in the next round
            minMz = precursorMz;

            if (usedEntries.find(entry->getLibFileOffset()) ==
                usedEntries.end()) { // check to see if this entry got picked as swap candidate in last batch
                entries.push_back(entry);
            }

        } else {

            entries.push_back(entry);

        }


    }

    // last batch, contains all the orphaned cases

    // cerr << "ORPHANS " << entries.size() << " left." << endl;
    for (vector<SpectraSTLibEntry *>::iterator i = entries.begin(); i != entries.end(); i++) {

        SpectraSTLibEntry *e1 = (*i);
        if (!e1) continue;

        // cerr << "Finding partner for " << e1->getName() << endl;

        for (vector<SpectraSTLibEntry *>::iterator j = entries.begin(); j != entries.end(); j++) {

            SpectraSTLibEntry *e2 = (*j);
            if ((!e2) || e1 == e2) continue;

            if (fabs(e1->getPrecursorMz() - e2->getPrecursorMz()) > 9.0 &&
                e1->getCharge() == e2->getCharge()) { // don't check homolog

                // cerr << "DONE: Swapping with " << e2->getName() << endl;

                swapPrecursors(e1, e2);

                usedEntries.insert(e1->getLibFileOffset());
                usedEntries.insert(e2->getLibFileOffset());

                delete (e1);
                delete (e2);

                e1 = NULL;
                e2 = NULL;
                (*i) = NULL;
                (*j) = NULL;
                break;
            }
        }

        if (e1) {
            // still can't find partner for e1
            // forget about them for now -- no decoy generated
            stringstream msg;
            msg << "Cannot find swapping partner for " << e1->getName() << " (";
            msg.precision(4);
            msg << e1->getPrecursorMz() << "). No decoy is generated for this spectrum.";
            g_log->log("DECOY", msg.str());

            if (m_params.decoyConcatenate) {
                insertOneEntry(e1, "DECOY");
            }

            delete (e1);
            e1 = NULL;
        }
    }


}

void SpectraSTSpLibImporter::swapPrecursors(SpectraSTLibEntry *e1, SpectraSTLibEntry *e2) {

    stringstream msg;
    msg << "Swapping precursors " << e1->getName() << " (";
    msg.precision(4);
    msg << fixed << e1->getPrecursorMz() << ") and " << e2->getName() << " (";
    msg.precision(4);
    msg << fixed << e2->getPrecursorMz() << ")";
    g_log->log("DECOY", msg.str());

    if (m_params.decoyConcatenate) {
        insertOneEntry(e1, "DECOY");
        insertOneEntry(e2, "DECOY");
    }

    SpectraSTLibEntry *decoy1 = new SpectraSTLibEntry(*e2);
    decoy1->setPrecursor(e1->getPrecursorMz(), e1->getCharge());
    decoy1->labelAsDecoy();
    // decoy1->getPeakList()->shiftAllPeaks(e1->getPrecursorMz() - e2->getPrecursorMz());
    insertOneEntry(decoy1, "DECOY");
    delete (decoy1);

    SpectraSTLibEntry *decoy2 = new SpectraSTLibEntry(*e1);
    decoy2->setPrecursor(e2->getPrecursorMz(), e2->getCharge());
    decoy2->labelAsDecoy();
    // decoy2->getPeakList()->shiftAllPeaks(e2->getPrecursorMz() - e1->getPrecursorMz());
    insertOneEntry(decoy2, "DECOY");
    delete (decoy2);

}

void
SpectraSTSpLibImporter::generatePrecursorSwapDecoys(SpectraSTMzLibIndex *mzIndex, vector<SpectraSTLibEntry *> &entries,
                                                    set<fstream::off_type> &usedEntries, double minMz, double maxMz) {

    // retrieve swapping candidates
    vector<SpectraSTLibEntry *> candidateEntries;
    mzIndex->retrieve(candidateEntries, minMz, maxMz, false);

    // cerr << "Found " << candidateEntries.size() << " candidates for range " << minMz << " to " << maxMz << "." << endl;

    // sort the candidates by charge
    map<int, list<SpectraSTLibEntry *> > candMap;
    for (vector<SpectraSTLibEntry *>::iterator ca = candidateEntries.begin(); ca != candidateEntries.end(); ca++) {

        if (usedEntries.find((*ca)->getLibFileOffset()) != usedEntries.end()) continue;
        if (!passAllFilters(*ca)) continue;

        int candidateCharge = (*ca)->getCharge();
        if (candidateCharge == 0) candidateCharge = 2;

        candMap[candidateCharge].push_back(*ca);
    }


    for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {

        SpectraSTLibEntry *entry = (*en);

        if (!passAllFilters(entry)) {
            // cerr << "FAIL FILTER" << endl;
            delete (entry);
            (*en) = NULL;
            continue;
        }

        // cerr << "Finding partner for " << entry->getName() << endl;

        Peptide *pep = entry->getPeptidePtr();
        int charge = entry->getCharge();


        if (candMap[charge].empty()) {
            // no candidate with a suitable charge...
            // cerr << "Failed " << entry->getName() << " (" << entry->getPrecursorMz() << ") No candidate for that charge." << endl;
            continue;
        }

        SpectraSTLibEntry *swap = NULL;
        unsigned int numTries = 0;

        do {

            numTries++;
            swap = candMap[charge].front();
            candMap[charge].pop_front();

            Peptide *swapPep = swap->getPeptidePtr();
            int dummy = 0;
            if (pep && swapPep && pep->isHomolog(*swapPep, 0.6, dummy)) {
                // cerr << "Homolog: " << pep->interactStyleWithCharge() << " | " << swapPep->interactStyleWithCharge() << endl;
                candMap[charge].push_back(swap);
                swap = NULL;
            }

        } while (!swap && numTries < candMap[charge].size());

        if (!swap) {
            // no candidate is not a homolog
            // cerr << "Failed " << entry->getName() << " (" << entry->getPrecursorMz() << ") No candidate not homolog." << endl;
            continue;
        }

        // cerr << "Done, swapping with " << swap->getName() << endl;

        swapPrecursors(entry, swap);

        usedEntries.insert(entry->getLibFileOffset());
        usedEntries.insert(swap->getLibFileOffset());

        delete (entry);
        (*en) = NULL;
    }

    // cerr << "Batch size = " << entries.size() << endl;

    vector<SpectraSTLibEntry *>::iterator i = entries.begin();
    while (i != entries.end()) {
        if (!(*i)) {
            // already done with this one, erase from vector
            i = entries.erase(i);
        } else {
            // cerr << "Failed: " << (*i)->getName() << endl;
            i++;
        }
    }

    // cerr << "Remaining size = " << entries.size() << endl;
    // at return time, the entries vector still contains failed cases (with non-null SpectraSTLibEntry*)

}


void SpectraSTSpLibImporter::doGenerateShuffleAndRepositionDecoy() {

    if (m_impFileNames.size() != 1) {
        g_log->error("DECOY", "Decoy generation must be applied to one .splib file only. No decoy library created.");
        return;
    }

    // don't do reduce in decoy generation - unclear how to reduce both
    // real and corresponding decoys
    m_params.reduceSpectrum = 0;

    if (!g_quiet) {
        cout << "Generating DECOY to \"" << m_impFileNames[0] << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(true, 20.0, true, true, true);

    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[0];
    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];

    if (!pepIndex) {
        // this means checkUniqueness in openSplibs failed.
        // non-unique library, i.e. some peptide ions have multiple spectra.
        g_log->error("DECOY", "Decoy generation requires unique library. Library \"" + m_impFileNames[0] +
                              "\" is non-unique. No decoy is generated.");
        return;
    }

    if (!mzIndex) {
        g_log->error("DECOY", "Cannot open the index of library \"" + m_impFileNames[0] + "\". Aborted.");
        return;
    }

    unsigned int count = mzIndex->getEntryCount();

    m_lib->writePreamble(m_preamble);

    ProgressCount pc(!g_quiet, 1, count);
    pc.start("Generating decoy spectra (shuffle and reposition method)");

    srand(time(NULL));

    // keep track of all sequences in library so that we won't shuffle to one that collides with a sequence already present
    map<int, map<string, vector<string> *> > allSequences; // length (NAA) -> set of all peptide sequences of that length

    // add all real sequences
    vector<string> allRealSequences;
    pepIndex->getAllSequences(allRealSequences);
    // cerr << "Loading real sequences " << allRealSequences.size();
    for (vector<string>::iterator i = allRealSequences.begin(); i != allRealSequences.end(); i++) {
        allSequences[(int) (i->length())][*i] = NULL;
    }

    set<string> allDecoySequences;

    SpectraSTLibEntry *entry;
    while (entry = mzIndex->nextEntry()) {

        Peptide *pep = entry->getPeptidePtr();

        if (!pep) {
            // not a peptide, no idea how to create decoy
            delete (entry);
            continue;
        }


        if (!passAllFilters(entry)) {
            delete (entry);
            continue;
        }

        pc.increment();

        if (m_params.decoyConcatenate) {
            insertOneEntry(entry, "DECOY");
        }

        map<string, vector<string> *>::iterator foundSeq = allSequences[pep->NAA()].find(pep->stripped);

        if (foundSeq == allSequences[pep->NAA()].end()) {
            // something is wrong, can't even find the sequence just put in...
            exit(1);
        }

        if (!(foundSeq->second)) {

            //      cerr << "Preparing to shuffle " << pep->interactStyleWithCharge() << endl;

            // first time this sequence is encountered

            Peptide *dummyPep = new Peptide(
                    *pep); // copy real peptide, this is just as a holder to be used in shuffling

            vector<SpectraSTLibEntry *> allIons;
            pepIndex->retrieve(allIons, pep->stripped, 0, "", entry->getFragType());

            for (vector<SpectraSTLibEntry *>::iterator i = allIons.begin(); i != allIons.end(); i++) {
                // this marks all amino acid that has been observed modified in this sequence
                // so that we won't shuffle those amino acid below
                //	cerr << " Adding mods of " << (*i)->getPeptidePtr()->interactStyleWithCharge() << endl;
                bool dummy = dummyPep->parseMspModStr((*i)->getPeptidePtr()->mspMods(), true);

            }

            //      cerr << " Dummy is " << dummyPep->interactStyleWithCharge() << endl;

            foundSeq->second = new vector<string>;

            for (int fold = 0; fold < m_params.decoySizeRatio; fold++) {

                bool okay = true;
                unsigned int numTries = 0;

                SpectraSTLibEntry *decoyEntry = NULL;
                Peptide *decoyPep = NULL;

                unsigned int numAAToAdd = 0;

                do {

                    Peptide *p = NULL;

                    numTries++;

                    p = dummyPep->shuffle(0.6, numAAToAdd);

                    if (numTries > 3) numAAToAdd = 1;
                    if (numTries > 6) numAAToAdd = 2;

                    okay = true;

                    //	  cerr << " ---> " << p->interactStyleWithCharge() << "  |" << numTries << endl;

                    if (allSequences[p->NAA()].find(p->stripped) != allSequences[p->NAA()].end()) {
                        // this shuffle already exists as a real sequence
                        //   cerr << " Already exists as real " << endl;
                        okay = false;
                    }
                    if (okay && allDecoySequences.find(p->stripped) != allDecoySequences.end()) {
                        // this shuffle already exists as a decoy sequence
                        //   cerr << " Already exists as decoy " << endl;
                        okay = false;
                    }

                    if (okay) {

                        decoyPep = new Peptide(p->stripped, pep->charge, pep->mspMods());
                        decoyPep->prevAA = pep->prevAA;
                        decoyPep->nextAA = pep->nextAA;
                        decoyEntry = new SpectraSTLibEntry(*entry); // copy
                        decoyEntry->makeDecoy(decoyPep, fold);
                        double precursorMz = decoyEntry->getPrecursorMz();

                        vector<SpectraSTLibEntry *> isobaricEntries;
                        mzIndex->retrieve(isobaricEntries, precursorMz - 2.0, precursorMz + 2.0, true);
                        // cerr << " Comparing to library " << isobaricEntries.size() << endl;;

                        for (vector<SpectraSTLibEntry *>::iterator ien = isobaricEntries.begin();
                             ien != isobaricEntries.end(); ien++) {
                            if (decoyEntry->getPeakList()->compare((*ien)->getPeakList()) > 0.7) {
                                // cerr << " Too similar to library spectrum " << (*ien)->getName() << endl;
                                okay = false;
                                break;
                            }
                        }
                    }

                    //	  numTries++;

                    delete (p);

                } while (!okay && numTries < 10);

                if (!okay) {
                    // change precursor m/z to the decoy sequence's (since AAs must have been added)
                    decoyEntry->synchWithPep();
                    // still not ok!
                    // delete (decoyEntry); // this will delete decoyPep also
                    // decoyEntry = NULL;
                    // cerr << "Shifting precursor m/z" << endl;
                    // cout << "Failed to generate decoy for " << pep->interactStyleWithCharge() << endl;
                    // break; // something is wrong -- probably still can't find one after 20 tries, no need to keep trying
                }

                allDecoySequences.insert(decoyPep->stripped);
                foundSeq->second->push_back(decoyPep->stripped);

                stringstream dss;
                dss << "Shuffle " << pep->stripped << " to ";
                dss << decoyPep->stripped << " .";
                if (pep->NAA() != decoyPep->NAA()) {
                    dss << " (" << decoyPep->NAA() - pep->NAA() << " AAs added randomly.)";
                }
                if (fabs(entry->getPrecursorMz() - decoyEntry->getPrecursorMz()) > 0.1) {
                    dss << " (Shift " << decoyEntry->getPrecursorMz() - entry->getPrecursorMz()
                        << " Th to decoy peptide's precursor m/z.)";
                }
                g_log->log("DECOY", dss.str());

                insertOneEntry(decoyEntry, "DECOY");

                delete (decoyEntry);
            } // END for (int fold = 0; fold < m_params.decoySizeRatio; fold++)

            delete (dummyPep);

        } else { // END if (!(*foundSeq))

            // this sequence has been encountered before and good shuffles were found
            vector<string> *goodShuffledSeqs = foundSeq->second;

            for (int fold = 0; fold < goodShuffledSeqs->size(); fold++) {

                Peptide *decoyPep = new Peptide((*goodShuffledSeqs)[fold], pep->charge, pep->mspMods());
                decoyPep->prevAA = pep->prevAA;
                decoyPep->nextAA = pep->nextAA;
                SpectraSTLibEntry *decoyEntry = new SpectraSTLibEntry(*entry); // copy
                decoyEntry->makeDecoy(decoyPep, fold);

                insertOneEntry(decoyEntry, "DECOY");
                delete (decoyEntry);

            }

        }

        delete (entry);

    }


    for (map<int, map<string, vector<string> *> >::iterator ii = allSequences.begin(); ii != allSequences.end(); ii++) {
        for (map<string, vector<string> *>::iterator jj = ii->second.begin(); jj != ii->second.end(); jj++) {
            if (jj->second) delete (jj->second);
        }
    }

    pc.done();


}

/*

void SpectraSTSpLibImporter::doGenerateDecoy() {
  
  if (m_impFileNames.size() != 1) {
    g_log->error("DECOY", "Decoy generation must be applied to one .splib file only. No decoy library created.");
    return;
  }

  // don't do reduce in decoy generation - unclear how to reduce both
  // real and corresponding decoys
  m_params.reduceSpectrum = 0;

  if (!g_quiet) {
    cout << "Generating DECOY to \"" << m_impFileNames[0] << "\"." << endl; 
  }

  // put the description in the preamble
  string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
  m_preamble.push_back(desc);
  g_log->log("CREATE", desc);

  openSplibs(false, 13.0, true, true, true);

  ifstream* splibFin = m_splibFins[0];
  if (!splibFin) return;
  SpectraSTPeptideLibIndex* pepIndex = m_pepIndices[0];

  if (splibFin && !pepIndex) {
    // this means checkUniqueness in openSplibs failed.
    // non-unique library, i.e. some peptide ions have multiple spectra.
    g_log->error("DECOY", "Decoy generation requires unique library. Library \"" + m_impFileNames[0] + "\" is non-unique. No decoy is generated.");
    return;
  }

  unsigned int count = pepIndex->getEntryCount();
  
  m_lib->writePreamble(m_preamble);
  
  ProgressCount pc(!g_quiet, 1, count);
  pc.start("Generating decoy spectra");

  srand(time(NULL));

  // keep track of all sequences in library so that we won't shuffle to one that collides with a sequence already present
  map<int, set<string> > allSequences; // length (NAA) -> set of all peptide sequences of that length
  // add all real sequences 
  vector<string> allRealSequences;
  pepIndex->getAllSequences(allRealSequences);
  for (vector<string>::iterator i = allRealSequences.begin(); i != allRealSequences.end(); i++) {
    allSequences[(int)(i->length())].insert(*i);
  }
  
  string origPeptide("");
  vector<string> subkeys;
  while (pepIndex->nextPeptide(origPeptide, subkeys)) {

    if (origPeptide[0] == '_') {
      // not a peptide, no idea how to create decoy
      continue;
    }
    
    Peptide p(origPeptide, 2, ""); // the stripped peptide

    vector<SpectraSTLibEntry*> entries(subkeys.size(), NULL);

    for (unsigned int sk = 0; sk < (unsigned int)(subkeys.size()); sk++) {

      vector<SpectraSTLibEntry*> entryHolder;
      pepIndex->retrieve(entryHolder, origPeptide, subkeys[sk]);

      if (entryHolder.size() != 1) {
        // unique lib; should only retrieve one entry -- so this shouldn't happen...
        continue;
      }

      if (!passAllFilters(entryHolder[0])) {
        delete (entryHolder[0]);
        continue;
      }

      pc.increment();

      entries[sk] = entryHolder[0];

      if (m_params.decoyConcatenate) {
        insertOneEntry(entries[sk], "DECOY");
      }

      // this marks all amino acid that has been observed modified in this stripped peptide
      // so that we won't shuffle those amino acid below
      int dummyCharge = 0;
      string dummyMods("");
      string dummyFrag("");
      SpectraSTPeptideLibIndex::parseSubkey(subkeys[sk], dummyCharge, dummyMods, dummyFrag);
      bool dummy = p.parseMspModStr(dummyMods, true);
    }

 
    // decoys[origPeptide] = 1;

    for (int fold = 0; fold < m_params.decoySizeRatio; fold++) {

      Peptide* decoyp = p.shufflePeptideSequence(allSequences);
      //  string decoyPeptide = p.reverse();

      allSequences[decoyp->NAA()].insert(decoyp->stripped);

      for (unsigned int sk = 0; sk < (unsigned int)(subkeys.size()); sk++) {

        if (!entries[sk]) {
          // the original entry doesn't pass filter, so won't create decoy either
          continue;
        }

	int decoyCharge = 0;
	string decoyMods("");
	string decoyFrag("");
	SpectraSTPeptideLibIndex::parseSubkey(subkeys[sk], decoyCharge, decoyMods, decoyFrag);

        Peptide* origPep = entries[sk]->getPeptidePtr();

	// SHUFFLE
	Peptide* decoyPep = new Peptide(decoyp->stripped, decoyCharge, decoyMods);

	// REVERSE
        // Peptide* decoyPep = new Peptide(*origPep);
        // decoyPeptide = decoyPep->reverse_sequence();

	// SHIFT
	// Peptide* decoyPep = new Peptide(*origPep);
	// 

        decoyPep->prevAA = origPep->prevAA;
        decoyPep->nextAA = origPep->nextAA;

        SpectraSTLibEntry* decoyEntry = new SpectraSTLibEntry(*(entries[sk]));
        decoyEntry->makeDecoy(decoyPep, (unsigned int)fold);

        stringstream dss;
        dss << "Shuffle " << origPep->interactStyleWithCharge() << " to ";
        dss << decoyPep->interactStyleWithCharge() << " .";
        if (origPep->NAA() != decoyPep->NAA()) {
          dss << " Two AAs added randomly.";
        }
        g_log->log("DECOY", dss.str());

	insertOneEntry(decoyEntry, "DECOY");
	
	delete (decoyEntry);
	
        if (fold >= m_params.decoySizeRatio - 1) {
          delete (entries[sk]);
        }

      }

      delete (decoyp);
    }

  }

  pc.done();
  
}

*/

// doSortByNreps - performs the build action SORT_BY_NREPS. Uses the mzIndex object to sort entries by Nreps and spit
// them up in descending order of Nreps.
void SpectraSTSpLibImporter::doSortByNreps() {

    if (m_impFileNames.size() != 1) {
        g_log->error("SORT_BY_NREPS", "Sorting by Nreps must be applied to one .splib file only. No library created.");
        return;
    }

    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(true, 13.0, false, false, true);
    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];
    if (!mzIndex) return;

    m_lib->writePreamble(m_preamble);

    if (!g_quiet) {
        cout << "SORT entries in \"" << m_impFileNames[0] << "\" by descending number of replicates...";
        cout.flush();
    }

    mzIndex->sortEntriesByNreps();

    if (!g_quiet) {
        cout << "DONE!" << endl;
    }

    ProgressCount pc(!g_quiet, 1, (int) (mzIndex->getEntryCount()));

    pc.start("Rewriting entries in the order of descending number of replicates");

    SpectraSTLibEntry *entry;
    while ((entry = mzIndex->nextSortedEntry())) {

        if (insertOneEntry(entry, "SORT_BY_NREPS")) {
            pc.increment();
        }
        delete (entry);
    }

    pc.done();

//  printProteinList();


}

// doUserSpecifiedModifications - create the semi-empirical spectra based on user-specified modifications
void SpectraSTSpLibImporter::doUserSpecifiedModifications() {

    if (m_impFileNames.size() != 1) {
        g_log->error("SEMI-EMPIRICAL",
                     "Semi-empirical spectrum generation must be applied to one .splib file only. No UserSpMods library created.");
        return;
    }

    if (m_params.allowableModTokens.empty()) {
        if (!m_params.predictionOrigTargetFileName.empty()) {
            buildSemiempiricalSpectraByOrigTargetList();
            return;
        }
        g_log->error("SEMI-EMPIRICAL",
                     "No user-specified modifications specified. Please use -cx option to list all allowable mod tokens. No action performed.");
        return;
    }

    m_params.reduceSpectrum = 0;

    if (!g_quiet) {
        cout << "Generating semi-empirical spectra for user-specified modifications for \"" << m_impFileNames[0]
             << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(false, 13.0, true, true, true);

    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[0];

    if (splibFin && !pepIndex) {
        // this means checkUniqueness in openSplibs failed.
        // non-unique library, i.e. some peptide ions have multiple spectra.
        g_log->error("SEMI-EMPIRICAL",
                     "Semi-empirical spectrum generation requires unique library. Library \"" + m_impFileNames[0] +
                     "\" is non-unique. No semi-empirical spectrum is generated.");
        return;
    }

    unsigned int count = pepIndex->getEntryCount();

    m_lib->writePreamble(m_preamble);

    // vector of several (aa => {set of mod tokens for that aa}) sets
    vector<map<char, set<string> > > allowableTokenSets;

    parseAllowableTokensStr(m_params.allowableModTokens, allowableTokenSets);

    if (!g_quiet) {
        unsigned int setCount = 0;
        for (vector<map<char, set<string> > >::iterator se = allowableTokenSets.begin();
             se != allowableTokenSets.end(); se++) {
            cout << "Allowable tokens (Set #" << ++setCount << "): ";
            for (map<char, set<string> >::iterator at = se->begin(); at != se->end(); at++) {
                string tok("");
                for (set<string>::iterator to = at->second.begin(); to != at->second.end(); to++) {
                    if (!tok.empty()) tok += "/";
                    tok += *to;
                }
                cout << " " << tok;
            }
            cout << "." << endl;
        }
    }

    ProgressCount pc2(!g_quiet, 1, count);
    pc2.start("Generating spectra");
    srand(time(NULL));

    // retrieve peptides from the splib and process the entries
    string origPeptide("");
    vector<string> subkeys;
    while (pepIndex->nextPeptide(origPeptide, subkeys)) {

        if (origPeptide[0] == '_') {
            // not a peptide
            pc2.increment();
            continue;
        }

        vector<SpectraSTLibEntry *> origEntries; // save pointers for deleting purposes
        map<string, pair<int, SpectraSTLibEntry *> > newIons; // peptide ion => (# AA changes, closest existing ion)

        for (vector<string>::iterator sk = subkeys.begin(); sk != subkeys.end(); sk++) {
            int charge = 0;
            string mod = "";
            string frag = "";
            SpectraSTPeptideLibIndex::parseSubkey(*sk, charge, mod, frag);

            //     if (!frag.empty() && frag != "CID") { // consider CID only for now
            //       continue;
            //     }

            vector<SpectraSTLibEntry *> entries;
            pepIndex->retrieve(entries, origPeptide, *sk);

            // should always retrieve exactly one entry, since this is a unique library
            SpectraSTLibEntry *entry = entries[0];

            origEntries.push_back(entry);

            vector<pair<string, int> > permutations;
            entry->getPeptidePtr()->permuteModTokens(allowableTokenSets, permutations);

            for (vector<pair<string, int> >::iterator pe = permutations.begin(); pe != permutations.end(); pe++) {

                map<string, pair<int, SpectraSTLibEntry *> >::iterator foundNewIon = newIons.find(pe->first);

                if (foundNewIon == newIons.end()) {
                    pair<int, SpectraSTLibEntry *> p;
                    p.first = pe->second;
                    p.second = entry;
                    newIons[pe->first] = p;

                } else {

                    if (foundNewIon->second.first > pe->second) {
                        // closer match
                        foundNewIon->second.first = pe->second;
                        foundNewIon->second.second = entry;
                    }
                }

            }

        }

        // done finding all permutations, now create semi-empirical spectra where needed
        for (map<string, pair<int, SpectraSTLibEntry *> >::iterator ni = newIons.begin(); ni != newIons.end(); ni++) {

            if (ni->second.first == 0) {
                // no change required. just insert original entry
                insertOneEntry(ni->second.second, "SEMI-EMPIRICAL");

            } else {

                Peptide *pep = new Peptide(ni->first, 0, "");

                // create semi-empirical spectrum
                SpectraSTLibEntry *closest = ni->second.second;

                Peptide *origPep = closest->getPeptidePtr();
                pep->prevAA = origPep->prevAA;
                pep->nextAA = origPep->nextAA;

                SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*closest);
                newEntry->makeSemiempiricalSpectrum(pep);

                stringstream dss;
                dss << "Perturb " << origPep->interactStyleWithCharge() << " to ";
                dss << pep->interactStyleWithCharge() << " .";
                g_log->log("SEMI-EMPIRICAL", dss.str());

                insertOneEntry(newEntry, "SEMI-EMPIRICAL");
                delete (newEntry);

            }
        } // for new ions of this peptide sequence

        // delete all old retrieved entries
        for (vector<SpectraSTLibEntry *>::iterator i = origEntries.begin(); i != origEntries.end(); i++) {
            delete (*i);
        }

        pc2.increment();
    } // while there are new peptide sequences

    pc2.done();

}

// parseAllowableTokenStr -- parse the user-specified string (option -cx) and put them in a data structure
void SpectraSTSpLibImporter::parseAllowableTokensStr(string &allowableTokensStr,
                                                     vector<map<char, set<string> > > &allowableTokenSets) {

    string tokenSetStr("");
    string::size_type bracePos = 0;
    while (!((tokenSetStr = nextToken(allowableTokensStr, bracePos, bracePos, "}\t\r\n", " {\t\r\n")).empty())) {

        bracePos++;

        string::size_type pos = 0;

        map<char, set<string> > tokenSet;

        while (pos != string::npos) {

            string token = Peptide::nextAAToken(tokenSetStr, pos, pos);

            if (token.empty() || (token[0] != 'n' && token[0] != 'c' && (token[0] < 'A' || token[0] > 'Z'))) break;

            if (token.length() > 1 && Peptide::modTokenTable->find(token) == Peptide::modTokenTable->end()) {
                g_log->error("SEMI-EMPIRICAL", "User-specified mod token \"" + token + "\" not recognized. Ignored.");
            } else {
                tokenSet[token[0]].insert(token);
            }
        }

        allowableTokenSets.push_back(tokenSet);
    }

}


// refresh - refresh protein mappings of peptide sequences in the library
void SpectraSTSpLibImporter::refresh() {

    if (!m_ppMappings) return;

    if (!g_quiet) {
        cout << "REFRESHING protein mappings...";
        cout.flush();
    }

    SpectraSTFastaFileHandler fasta(m_params.refreshDatabase);

    fasta.refresh(*m_ppMappings, m_params.refreshTrypticOnly);

    if (!g_quiet) {
        cout << "DONE!" << endl;
    }
}

// addSequencesForRefresh - add new sequences to be refreshed (into m_ppMappings)
void SpectraSTSpLibImporter::addSequencesForRefresh(vector<string> &seqs) {
    if (!m_ppMappings) {
        m_ppMappings = new map<string, vector<pair<string, string> > *>;
    }
    for (vector<string>::iterator i = seqs.begin(); i != seqs.end(); i++) {
        (*m_ppMappings)[*i] = NULL; // no protein mapped for now
    }
}

/*
// passAllFilters - check that an entry passes all the filters (-cf, -cT, -cO, -cd/-cu options)
bool SpectraSTSpLibImporter::passAllFilters(SpectraSTLibEntry* entry) {

  if (!m_ppMappings) {
    return (SpectraSTLibImporter::passAllFilters(entry));
  }
  
  if (!(entry->getPeptidePtr())) {
    return (SpectraSTLibImporter::passAllFilters(entry));
  }
  
  if (SpectraSTLibImporter::passAllFilters(entry)) {

    map<string, vector<pair<string, string> >* >::iterator found = m_ppMappings->find(entry->getPeptidePtr()->stripped);
  
    if (found == m_ppMappings->end()) {
      return (false);
    } else if (m_params.refreshDeleteUnmapped && !(found->second)) {
      return (false);
    } else if (m_params.refreshDeleteMultimapped && found->second && found->second->size() != 1) {
      return (false);
    }
    
    return (true);
    
  }

  return (false);
}
*/

bool SpectraSTSpLibImporter::insertOneEntry(SpectraSTLibEntry *entry, string mode) {

    if (!(SpectraSTLibImporter::passAllFilters(entry))) return (false);

    if (entry->getPeakList()->getNumPeaks() < m_params.minimumNumPeaksToInclude) {
        g_log->log(mode, "Too few peaks in spectrum. Skipped entry \"" + entry->getName() + "\".");
        return (false);
    }

    Peptide *pep = entry->getPeptidePtr();

    if (pep && pep->NAA() < m_params.minimumNumAAToInclude) {
        g_log->log(mode, "Peptide ID too short. Skipped entry \"" + entry->getName() + "\".");
        return (false);
    }

    entry->setOneComment("Spec", mode);

    if (pep && m_ppMappings) {

        map<string, vector<pair<string, string> > *>::iterator found = m_ppMappings->find(pep->stripped);

        if (found != m_ppMappings->end()) {

            stringstream protss;
            stringstream contextss;

            if (found->second) {
                // found mapping

                if (m_params.refreshDeleteMultimapped && found->second->size() != 1) {
                    return (false);
                }

                protss << found->second->size();
                contextss << found->second->size();


                char origPrevAA = pep->prevAA;
                char origNextAA = pep->nextAA;
                char bestPrevAA = origPrevAA;
                char bestNextAA = origNextAA;
                int origNTT = (int) (pep->NTT());
                int highestNTT = -1;
                bool foundOrig = false;

                string proteins("");
                string contexts("");

                for (vector<pair<string, string> >::iterator prot = found->second->begin();
                     prot != found->second->end(); prot++) {

                    if (proteins.empty()) {
                        proteins = prot->first;
                        contexts = prot->second;
                    } else {
                        if (prot->first.compare(0, 5, "DECOY") == 0 || prot->first.compare(0, 3, "REV") == 0 ||
                            prot->first.compare(0, 3, "rev") == 0) {
                            proteins = proteins + "/" + prot->first;
                            contexts = contexts + "/" + prot->second;
                        } else {
                            proteins = prot->first + "/" + proteins;
                            contexts = prot->second + "/" + contexts;
                        }
                    }

                    if (origPrevAA == prot->second[3] && origNextAA == prot->second[5]) {
                        foundOrig = true;
                        if (origNTT > highestNTT) highestNTT = origNTT;
                        continue;
                    }

                    if (foundOrig && origNTT == 2) {
                        continue;
                    }

                    Peptide testPep(pep->stripped, 1);
                    testPep.prevAA = prot->second[3];
                    testPep.nextAA = prot->second[5];
                    int testNTT = (int) (testPep.NTT());

                    if (testNTT > highestNTT) {
                        highestNTT = testNTT;
                        bestPrevAA = prot->second[3];
                        bestNextAA = prot->second[5];
                    }

                }

                protss << "/" << proteins;
                contextss << "/" << contexts;

                if (!foundOrig || highestNTT > origNTT) {
                    pep->prevAA = bestPrevAA;
                    pep->nextAA = bestNextAA;
                    entry->synchWithPep();
                }

            } else {

                if (m_params.refreshDeleteUnmapped) {
                    return (false);
                }

                protss << "0/UNMAPPED";
                string origProtein("");
                if (entry->getOneComment("Protein", origProtein)) {
                    entry->setOneComment("OrigProtein", origProtein);
                }
                contextss << "0/UNMAPPED";

            }

            entry->setOneComment("Protein", protss.str());
            entry->setOneComment("PepContext", contextss.str());

        } else {
            // this peptide can't be found in the ppMappings, maybe it's a decoy or semi-empirical, don't touch it
        }

    }

    if (!(m_params.setFragmentation.empty())) {
        entry->setFragType(m_params.setFragmentation);
    }

    if (m_params.reannotatePeaks) {
        entry->annotatePeaks(true);
    }

    string naaStr;
    if (pep && !(entry->getOneComment("NAA", naaStr))) {
        stringstream naass;
        naass << entry->getPeptidePtr()->NAA();
        entry->setOneComment("NAA", naass.str());
    }

    string specType("");
    if (m_params.plotSpectra == "ALL" || m_params.plotSpectra == entry->getStatus() ||
        (entry->getOneComment("Spec", specType) && !m_params.plotSpectra.empty() && m_params.plotSpectra == specType)) {
        plot(entry);
    }

    if (m_params.reduceSpectrum > 0) {
        double reducedFraction = entry->getPeakList()->reduce(m_params.reduceSpectrum, m_params.minimumMRMQ3MZ,
                                                              m_params.maximumMRMQ3MZ, entry->getNrepsUsed());
        stringstream ricss;
        ricss.precision(3);
        ricss << fixed << showpoint << reducedFraction;
        entry->setOneComment("ReducedFracIonCurrent", ricss.str());

    }

    m_lib->insertEntry(entry);
    return (true);

}

/*
// processEntry - process an entry before it is inserted into a new library
void SpectraSTSpLibImporter::processEntry(SpectraSTLibEntry* entry) {

  if (!(m_params.setFragmentation.empty())) {
    entry->setFragType(m_params.setFragmentation);
  }

  if (m_params.reannotatePeaks) {
    entry->annotatePeaks(true);
  }
  
  string naaStr;
  if (!(entry->getOneComment("NAA", naaStr)) && entry->getPeptidePtr()) {
    stringstream naass;
    naass << entry->getPeptidePtr()->NAA();
    entry->setOneComment("NAA", naass.str());
  }

  string specType("");
  if (m_params.plotSpectra == "ALL" || m_params.plotSpectra == entry->getStatus() || 
      (entry->getOneComment("Spec", specType) && !m_params.plotSpectra.empty() && m_params.plotSpectra == specType)) {
    plot(entry);
  }
	
  if (m_params.reduceSpectrum > 0) {      
    double reducedFraction = entry->getPeakList()->reduce(m_params.reduceSpectrum, m_params.minimumMRMQ3MZ, m_params.maximumMRMQ3MZ, entry->getNrepsUsed());
    stringstream ricss;
    ricss.precision(3);
    ricss << fixed << showpoint << reducedFraction;
    entry->setOneComment("ReducedFracIonCurrent", ricss.str());
    
  }

  // do the protein remapping here
  if (m_ppMappings && entry->getPeptidePtr()) {

    // hijack - remove RawSpectrum info
    // entry->deleteOneComment("RawSpectrum");
    // entry->deleteOneComment("BestRawSpectrum");
    // END hijack
    
    map<string, vector<pair<string, string> >* >::iterator found = m_ppMappings->find(entry->getPeptidePtr()->stripped);

    if (found != m_ppMappings->end()) {

      stringstream protss;
      stringstream contextss;

      if (found->second) {
	// found mapping
	protss << found->second->size();
	contextss << found->second->size();


	char origPrevAA = entry->getPeptidePtr()->prevAA;
	char origNextAA = entry->getPeptidePtr()->nextAA;
	char bestPrevAA = origPrevAA;
	char bestNextAA = origNextAA;
	int origNTT = (int)(entry->getPeptidePtr()->NTT());
	int highestNTT = -1;
	bool foundOrig = false;

	string proteins("");
	string contexts("");
	
	for (vector<pair<string, string> >::iterator prot = found->second->begin(); prot != found->second->end(); prot++) {
	  
	  if (proteins.empty()) {
	    proteins = prot->first;
	    contexts = prot->second;
	  } else {
	    if (prot->first.compare(0, 5, "DECOY") == 0 || prot->first.compare(0, 3, "REV") == 0 || prot->first.compare(0, 3, "rev") == 0) {
	      proteins = proteins + "/" + prot->first;
	      contexts = contexts + "/" + prot->second;
	    } else {	
	      proteins = prot->first + "/" + proteins;
	      contexts = prot->second + "/" + contexts;
	    }
	  }
	  
	  if (origPrevAA == prot->second[3] && origNextAA == prot->second[5]) {
	    foundOrig = true;
	    if (origNTT > highestNTT) highestNTT = origNTT;
	    continue;
	  }

	  if (foundOrig && origNTT == 2) {
	    continue;
	  }
	   
	  Peptide testPep(entry->getPeptidePtr()->stripped, 1);
	  testPep.prevAA = prot->second[3];
	  testPep.nextAA = prot->second[5];
	  int testNTT = (int)(testPep.NTT());
	  
	  if (testNTT > highestNTT) {
	    highestNTT = testNTT;
	    bestPrevAA = prot->second[3];
	    bestNextAA = prot->second[5];
	  }

	}
	
	protss << "/" << proteins;
	contextss << "/" << contexts;
	
	if (!foundOrig || highestNTT > origNTT) {
	  entry->getPeptidePtr()->prevAA = bestPrevAA;
	  entry->getPeptidePtr()->nextAA = bestNextAA;
	  entry->synchWithPep();
	}

      } else {
	protss << "0/UNMAPPED";
	string origProtein("");
	if (entry->getOneComment("Protein", origProtein)) {
	  entry->setOneComment("OrigProtein", origProtein);
	}
	contextss << "0/UNMAPPED";

      }
    
      entry->setOneComment("Protein", protss.str());
      entry->setOneComment("PepContext", contextss.str());
    }
  }
}
*/

// constructFileListStr - make a string of files processed to be printed to the preamble and to the log file
string SpectraSTSpLibImporter::constructFileListStr() {

    stringstream ss;
    string fullName(m_impFileNames[0]);
    makeFullPath(fullName);
    ss << "\"" << fullName << "\" ";

    if (m_impFileNames.size() > 9) {
        fullName = (*(m_impFileNames.end() - 1));
        makeFullPath(fullName);
        ss << m_params.combineAction << " ... \"" << fullName << "\" ";
    } else {
        for (vector<string>::iterator i = m_impFileNames.begin() + 1; i != m_impFileNames.end(); i++) {
            fullName = *i;
            makeFullPath(fullName);
            ss << m_params.combineAction << " \"" << fullName << "\" ";
        }
    }
    return (ss.str());

}

// openSplibs - helper method to handle the opening of splib libraries (and their indices) to be imported. Opens the
// ifstream on each .splib file, parse the preamble, and optionally open also the associated mz index and/or peptide index.
// If checkUniqueness and/or refresh is TRUE, the peptide index is needed anyway, so it 
// will be opened whether or not openPepIndex is TRUE.
// The opened objects will be stored in m_splibFins, m_mzIndices and m_pepIndices (the objects of the same library 
// will occupy the same index as that of the file name in m_impFileNames), to be deleted in the destructor.

void SpectraSTSpLibImporter::openSplibs(bool openMzIndex, double mzIndexCacheRange,
                                        bool openPepIndex, bool checkUniqueness, bool refresh) {

    // just to make sure this method is only called once per instance
    if (!(m_splibFins.empty())) return;

    bool hasFailure = false;

    for (vector<string>::iterator f = m_impFileNames.begin(); f != m_impFileNames.end(); f++) {

        FileName fn;
        parseFileName((*f), fn);

        ifstream *splibFin = new ifstream;

        if (!myFileOpen(*splibFin, (*f), true)) {
            g_log->error("CREATE", "Cannot open SPLIB file \"" + *f + "\" for reading. File skipped.");
            delete (splibFin);
            m_splibFins.push_back(NULL);
            m_pepIndices.push_back(NULL);
            m_mzIndices.push_back(NULL);
            continue;
        }

        // peek to see if it's a binary file or not
        bool binary = true;
        char firstChar = splibFin->peek();
        if (firstChar == '#' || firstChar == 'N') {
            binary = false;
        }

        m_splibFins.push_back(splibFin);

        if (openPepIndex || checkUniqueness || (refresh && !(m_params.refreshDatabase.empty()))) {
            SpectraSTPeptideLibIndex *pepIndex = new SpectraSTPeptideLibIndex(fn.path + fn.name + ".pepidx", splibFin,
                                                                              binary);
            if (checkUniqueness && !pepIndex->isUniqueLibrary()) {
                // non-unique library, i.e. some peptide ions have multiple spectra.
                delete (pepIndex);
                m_pepIndices.push_back(NULL);

            } else {

                m_pepIndices.push_back(pepIndex);

                if (refresh && !(m_params.refreshDatabase.empty())) {
                    vector<string> seqs;
                    pepIndex->getAllSequences(seqs);
                    addSequencesForRefresh(seqs);
                }
            }

        } else {
            m_pepIndices.push_back(NULL);
        }

        if (openMzIndex) {
            SpectraSTMzLibIndex *mzIndex = new SpectraSTMzLibIndex(fn.path + fn.name + ".spidx", splibFin,
                                                                   mzIndexCacheRange, binary);
            m_mzIndices.push_back(mzIndex);
        } else {
            m_mzIndices.push_back(NULL);
        }

        // parse the preambles of the imported .splib files - these will be appended to the
        // preamble of the generated .splib file
        parsePreamble(*splibFin, binary);


    }

}

void SpectraSTSpLibImporter::doSimilarityClustering() {

    if (m_impFileNames.size() != 1) {
        g_log->error("SIMILARITY_CLUSTERING",
                     "Similarity clustering must be applied to one .splib file only. No library created.");
        return;
    }

    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(true, 13.0, false, false, false);
    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTMzLibIndex *mzIndex = m_mzIndices[0];
    if (!mzIndex) return;

    m_lib->writePreamble(m_preamble);

    if (!g_quiet) {
        cout << "CLUSTER entries in \"" << m_impFileNames[0] << "\" by spectral similarity." << endl;
        cout.flush();
    }

    map<fstream::off_type, int> clusteredEntries;
    vector<set<fstream::off_type> *> multiclusters;

    ProgressCount pc(!g_quiet && !g_verbose, 1, (int) (mzIndex->getEntryCount()));
    pc.start("Clustering");

    mzIndex->sortEntriesBySN();

    fstream::off_type offset = 0;
    while (mzIndex->nextSortedFileOffset(offset)) {

        pc.increment();

        if (clusteredEntries.find(offset) != clusteredEntries.end()) {
            // already in a cluster
            continue;
        }

        SpectraSTLibEntry *entry = mzIndex->thisSortedEntry();

        if (entry->getLibFileOffset() != offset) {
            cerr << "Something is wrong!!" << endl;
            exit(1);
        }

        // if (!(entry->getPeakList()->passFilterUnidentified())) {
        //  // this is just double-checking, should have been removed in MzXMLLibImporter
        //  delete (entry);
        //  continue;
        // }

        set<fstream::off_type> *cluster = new set<fstream::off_type>;
        cluster->insert(entry->getLibFileOffset());

        vector<SpectraSTLibEntry *> isobaricEntries;
        double rootPrecursorMz = entry->getPrecursorMz();
        mzIndex->retrieve(isobaricEntries, rootPrecursorMz - 2.5, rootPrecursorMz + 2.5, true);

        // remove spectra that are already members of other clusters
        for (vector<SpectraSTLibEntry *>::iterator ie = isobaricEntries.begin(); ie != isobaricEntries.end(); ie++) {
            if (clusteredEntries.find((*ie)->getLibFileOffset()) != clusteredEntries.end()) {
                // already member of other clusters
                (*ie) = NULL;
            }
        }

        findSpectralNeighbors(entry, rootPrecursorMz, 0, isobaricEntries, cluster);

        if (cluster->size() == 1) {
            // singleton, just copy this entry and be done with it

            string xreaStr("");
            double xrea = 0.0;
            if (entry->getOneComment("Xrea", xreaStr)) {
                xrea = atof(xreaStr.c_str());
            } else {
                xrea = entry->getPeakList()->calcXrea(true);
                stringstream xreass;
                xreass.precision(3);
                xreass << fixed << xrea;
                entry->setOneComment("Xrea", xreass.str());
            }

            //string snStr("");
            //double sn = 0.0;
            //if (entry->getOneComment("SN", snStr)) {
            //  sn = atof(snStr.c_str());
            //}

            if (entry->getNrepsUsed() > 1 || xrea >= m_params.unidentifiedSingletonXreaThreshold) {
                // re-read unprocessed entry from library
                splibFin->seekg(entry->getLibFileOffset());
                SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*splibFin, true, false);
                insertOneEntry(newEntry, "SIMILARITY_CLUSTERING");
                delete (newEntry);
            }

            clusteredEntries[entry->getLibFileOffset()] = -1; // don't belong to any multicluster
            delete (cluster);

        } else {
            int clusterIndex = (int) (multiclusters.size());
            for (set<fstream::off_type>::iterator i = cluster->begin(); i != cluster->end(); i++) {
                clusteredEntries[*i] = clusterIndex;
            }
            multiclusters.push_back(cluster);
        }

        delete (entry);

    }

    pc.done();

    cout << "Found " << multiclusters.size() << " clusters of 2+ members." << endl;

    ProgressCount pc2(!g_quiet && !g_verbose, 1, (int) (multiclusters.size()));
    pc2.start("Generating merged spectra from clusters");

    // now deal with true clusters of more than 1 members
    for (vector<set<fstream::off_type> *>::iterator cl = multiclusters.begin(); cl != multiclusters.end(); cl++) {

        pc2.increment();

        vector<SpectraSTLibEntry *> entries;
        for (set<fstream::off_type>::iterator os = (*cl)->begin(); os != (*cl)->end(); os++) {
            splibFin->seekg(*os);
            SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*splibFin, true, false);
            entries.push_back(newEntry);
        }

        SpectraSTReplicates *replicates = new SpectraSTReplicates(entries, m_params);

        SpectraSTLibEntry *consensus = replicates->makeConsensusSpectrum();

        if (consensus) {
            insertOneEntry(consensus, "SIMILARITY_CLUSTERING");
        }

        delete (replicates);

        for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {
            delete (*en);
        }

        delete (*cl);

    }

    pc2.done();


}

void SpectraSTSpLibImporter::findSpectralNeighbors(SpectraSTLibEntry *entry, double rootPrecursorMz, unsigned int round,
                                                   vector<SpectraSTLibEntry *> &isobaricEntries,
                                                   set<fstream::off_type> *cluster) {

    SpectraSTPeakList *pl = entry->getPeakList();
    pl->quickSimplify(50, 99999, true, 0.0);
//    int level = cluster->size();

    double lowMz = rootPrecursorMz - 2.5 + (double) round;
    double highMz = rootPrecursorMz + 2.5 - (double) round;

    unsigned int numInCluster = (unsigned int) (cluster->size());
    double sumMzInCluster = rootPrecursorMz * (double) (cluster->size());
    vector<SpectraSTLibEntry *> hits;
    //  double sumDot = 0.0;
    //  double sumDotSquare = 0.0;
    // multimap<double, SpectraSTLibEntry*> hits;

    for (vector<SpectraSTLibEntry *>::iterator en = isobaricEntries.begin(); en != isobaricEntries.end(); en++) {

        if (!(*en)) {
            // do not consider, set to NULL because it is already member of another cluster, or
            // already matches one of the spectrum in this cluster with dot < 0.3
            continue;
        }

        if ((cluster->find((*en)->getLibFileOffset()) != cluster->end())) {
            // already in this cluster
            (*en) = NULL;
            continue;
        }

        double precursorMz = (*en)->getPrecursorMz();

        if (precursorMz < lowMz || precursorMz > highMz) {
            continue;
        }

        (*en)->getPeakList()->quickSimplify(50, 99999, true, 0.0);
        double dot = pl->compare((*en)->getPeakList());

        if (dot >= m_params.unidentifiedClusterMinimumDot - (double) round * 0.05) {

            cluster->insert((*en)->getLibFileOffset());
            hits.push_back(*en);

            numInCluster++;
            sumMzInCluster += precursorMz;

        } else if (dot < 0.3) {
            // hopeless, remove these from consideration in subsequent rounds
            (*en) = NULL;

        }

    }

    if (round >= 2) return;

    double meanMzInCluster = sumMzInCluster / (double) numInCluster;

    for (vector<SpectraSTLibEntry *>::iterator h = hits.begin(); h != hits.end(); h++) {

        findSpectralNeighbors(*h, meanMzInCluster, round + 1, isobaricEntries, cluster);

    }

}

bool SpectraSTSpLibImporter::hackDeamidation(SpectraSTLibEntry *entry) {

    int numDeamidation = entry->getMassDiffInt();
    cerr << "Processing " << entry->getName() << ":" << numDeamidation << endl;
    if (numDeamidation > 0) {
        Peptide *p = entry->getPeptidePtr();
        int numPossibleSites = 0;
        for (int pos = 0; pos < (int) (p->NAA()); pos++) {
            if (p->mods.find(pos) == p->mods.end() && (p->stripped[pos] == 'N' || p->stripped[pos] == 'Q')) {
                // could be deamidated
                numPossibleSites++;
            }
        }
        if (numDeamidation >= numPossibleSites) {
            for (int pos = 0; pos < (int) (p->NAA()); pos++) {
                if (p->mods.find(pos) == p->mods.end() && (p->stripped[pos] == 'N' || p->stripped[pos] == 'Q')) {
                    p->mods[pos] = "Deamidated";
                }
            }
            if (numDeamidation == numPossibleSites) {
                cerr << "UNAMBIGUOUS SITES -- " << entry->getName() << endl;
            }
            if (numDeamidation > numPossibleSites) {
                // something's wrong, can't explain by deamidation. Just flag this
                cerr << "NOT ENOUGH SITES -- " << entry->getName() << ":" << numDeamidation << "," << numPossibleSites
                     << endl;
            }
        } else {
            // uncertain sites. this time we will just throw this spectrum away
            cerr << "TOO MANY SITES -- " << entry->getName() << ":" << numDeamidation << "," << numPossibleSites
                 << endl;
            return (false);
        }

        entry->synchWithPep();
        entry->getPeakList()->setPeptidePtr(p);
        entry->getPeakList()->annotate(true);
    }

    return (true);
}

void SpectraSTSpLibImporter::readDenoiserParamFile(string &fileName) {

    ifstream fin;

    if (!(myFileOpen(fin, fileName, false))) return;

    string line("");
    while (!(nextLineCropped(fin, line, "CHARGE", ""))) {
        cout << "Reading: " << line << endl;
        if (line == "_EOF_") break;

        int charge = 0;
        if (line.length() > 6) {
            charge = atoi(line.substr(6).c_str());
            cout << "Got charge: " << charge << endl;
            (*m_denoisers)[charge]->readModel(fin);
        }
    }

}

void SpectraSTSpLibImporter::writeDenoiserParamFile(string &fileName) {

    ofstream fout;

    if (!(myFileOpen(fout, fileName, false))) return;

    for (unsigned int charge = 0; charge <= MAX_CHARGE; charge++) {
        fout << "CHARGE" << charge << endl;
        (*m_denoisers)[charge]->printModel(fout);
        fout << "END CHARGE" << endl;
    }

}

void SpectraSTSpLibImporter::buildSemiempiricalSpectraByOrigTargetList() {

    if (!g_quiet) {
        cout << "Generating semi-empirical spectra for user-specified modifications for \"" << m_impFileNames[0]
             << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(false, 13.0, true, true, true);

    ifstream *splibFin = m_splibFins[0];
    if (!splibFin) return;
    SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[0];

    if (splibFin && !pepIndex) {
        // this means checkUniqueness in openSplibs failed.
        // non-unique library, i.e. some peptide ions have multiple spectra.
        g_log->error("SEMI-EMPIRICAL",
                     "Semi-empirical spectrum generation requires unique library. Library \"" + m_impFileNames[0] +
                     "\" is non-unique. No semi-empirical spectrum is generated.");
        return;
    }

    unsigned int count = pepIndex->getEntryCount();

    m_lib->writePreamble(m_preamble);

    ProgressCount pc2(!g_quiet, 1, count);
    pc2.start("Generating spectra");
    srand(time(NULL));

    ifstream fin;
    string line;
    myFileOpen(fin, m_params.predictionOrigTargetFileName);
    while (nextLine(fin, line)) {

        string::size_type pos = 0;
        string opep = nextToken(line, pos, pos);
        string tpep = nextToken(line, pos, pos);

        Peptide *origPeptide = new Peptide(opep, 0);
        Peptide *tarPeptide = new Peptide(tpep, 0);

        vector<SpectraSTLibEntry *> entries;
        pepIndex->retrieve(entries, origPeptide->stripped, origPeptide->charge, origPeptide->mspMods());

        SpectraSTLibEntry *origEntry = entries[0];
        SpectraSTLibEntry *tarEntry = new SpectraSTLibEntry(*origEntry);

        tarEntry->makeSemiempiricalSpectrum(tarPeptide);
        insertOneEntry(tarEntry, "SEMI-EMPIRICAL");

        delete (tarEntry);
        delete origPeptide;
    }

    fin.close();
    return;
}

void SpectraSTSpLibImporter::doBuildSemiEmpiricalSplib() {

    if (m_impFileNames.size() != 1) {
        g_log->error("SEMI-EMPIRICAL",
                     "Semi-empirical spectrum generation must be applied to one .splib file only. No library created.");
        return;
    }

    if (!g_quiet) {
        cout << "Generating semi-empirical spectra for user-specified modifications for \"" << m_impFileNames[0]
             << "\"." << endl;
    }

    // put the description in the preamble
    string desc = m_params.constructDescrStr(constructFileListStr(), ".splib");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    openSplibs(false, 0.0, true, false, true);

    m_lib->writePreamble(m_preamble);

    // read the target peptides
    ifstream fin;
    string line;
    string fn = m_params.predictionTargetPeptidesFileName;

    myFileOpen(fin, fn);
    ProgressCount pc(!g_quiet && !g_verbose, 500, 0);
    pc.start("Build Semi-empirical Spectral Library");

    while (nextLine(fin, line)) {

        pc.increment();
        if ((line[0] > 'Z') || (line[0] < 'A')) {
            continue;
        }
        string targetPeptide;
        string targetSubkey;
        int targetCharge;
        string targetMods;
        parseLineOfPredictionTargetPeptide(line, targetPeptide, targetSubkey, targetCharge, targetMods);

        //cout << "targetPeptide=" << targetPeptide << endl;
        //cout << "targetSubkey=" << targetSubkey << endl;
        //cout << "targetCharge=" << targetCharge << endl;
        //cout << "targetMods=" << targetMods << endl;

        string peptide;
        string mods;
        vector<string> subkeys;

        vector<SpectraSTLibEntry *> hits;

        int isFound = 0;// find the exact matched peptide entry
        vector<SpectraSTPeptideLibIndex *>::size_type curPepIndex = 0;
        SpectraSTPeptideLibIndex *pepIndex = m_pepIndices[curPepIndex];

        while (isFound == 0) {
            if (curPepIndex >= m_pepIndices.size()) {
                break;
            }

            m_pepIndices[curPepIndex]->reset();
            while ((pepIndex = m_pepIndices[curPepIndex]) &&
                   (pepIndex->nextPeptide(peptide, subkeys)) && (!isFound)) {
                for (vector<string>::iterator s = subkeys.begin(); s != subkeys.end(); s++) {
                    //cout << "nextpeptide=" << peptide << "\t" << *s << endl;
                    int distance = calcPredictionDistance(targetPeptide, targetSubkey, peptide, *s);
                    //cout << "distance=" << distance << endl;
                    if (distance == 0) {
                        // already find the correct answer, no need further prediction
                        if (hits.size() > 0) {
                            // cout << "hits.size=" << hits.size() << endl;
                            for (vector<SpectraSTLibEntry *>::iterator h = hits.begin(); h != hits.end(); h++) {
                                delete (*h);
                            }
                        }
                        //cout << "hits.size()=" << hits.size()<< endl;
                        hits.clear();
                        //cout << "hits.clear()" << endl;
                        //cout << "peptide=" << peptide << "\tsubkey=" << *s << endl;
                        pepIndex->retrieve(hits, peptide, *s);
                        insertOneEntry(hits[0], "SEMI-EMPIRICAL");
                        isFound = 1;
                        break;

                    } else if (distance == 1) {
                        //	  } else if ((distance <= m_params.allowablePredictionDistance) && (distance > 0)) {
                        pepIndex->retrieve(hits, peptide, *s);
                        isFound = 2;
                        break;
                    } else {
                        if (distance == -1) {
                            break;
                        }
                        // unexpected cases, do nothing
                    }
                }

            }
            curPepIndex++;
        }

        //cout << "exit while" << endl;
        // if (isFound) continue;
        if (isFound == 1) {
            continue;
        }

        if (hits.size() == 0) {
            // no candidates found
            g_log->log("SEMI-EMPIRICAL", "Cannot find template for predicting \"" + line + "\".");
            continue;
        }

        vector<SpectraSTLibEntry *> predicts;

        for (vector<SpectraSTLibEntry *>::iterator h = hits.begin(); h != hits.end(); h++) {
            Peptide *pep = new Peptide(targetPeptide, targetCharge, targetMods);

            // create semi-empirical spectrum
            SpectraSTLibEntry *closest = *h;

            Peptide *origPep = closest->getPeptidePtr();
            pep->prevAA = origPep->prevAA;
            pep->nextAA = origPep->nextAA;

            SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*closest);
            newEntry->makeSemiempiricalSpectrum(pep);

            stringstream dss;
            dss << "Perturb " << origPep->interactStyleWithCharge() << " to ";
            dss << pep->interactStyleWithCharge() << " .";
            g_log->log("SEMI-EMPIRICAL", dss.str());

            //insertOneEntry(newEntry, "SEMI-EMPIRICAL");
            predicts.push_back(newEntry);
        }

        for (vector<SpectraSTLibEntry *>::iterator h = hits.begin(); h != hits.end(); h++) {
            delete (*h);
        }

        if (predicts.size() == 1) {

            insertOneEntry(predicts[0], "SEMI-EMPIRICAL");

        } else if (predicts.size() > 1) {

            SpectraSTReplicates *replicates = new SpectraSTReplicates(predicts, m_params, m_denoisers);
            // make a consensus spectrum of replicates and insert that into the library
            SpectraSTLibEntry *consensus = replicates->makeConsensusSpectrum();

            insertOneEntry(consensus, "SEMI-EMPIRICAL");
            delete (replicates);
        }

        for (vector<SpectraSTLibEntry *>::iterator p = predicts.begin(); p != predicts.end(); p++) {
            delete (*p);
        }


    }

    pc.done();
}


void
SpectraSTSpLibImporter::parseLineOfPredictionTargetPeptide(string line, string &peptide, string &subkey, int &charge,
                                                           string &modstr) {

    vector<string> *li = split(line, " \t", " \t");
    peptide = (*li)[0];
    subkey = (*li)[1];

    vector<string> *si = split(subkey, "|", "|");
    charge = atoi((*si)[0].c_str());
    modstr = (*si)[1];

    li->clear();
    delete li;
    si->clear();
    delete si;

}

int SpectraSTSpLibImporter::calcPredictionDistance(string targetPeptide, string targetSubkey, string peptide,
                                                   string subkey) {

    if (targetPeptide.size() != peptide.size()) {
        return (-1);
    }

    //bool hasInvalidDistanceType = false;
    vector<string> *si = split(subkey, "|", "|");
    int charge = atoi((*si)[0].c_str());
    string modstr = (*si)[1];

    vector<string> *tsi = split(targetSubkey, "|", "|");
    int tc = atoi((*tsi)[0].c_str());
    string tmodstr = (*tsi)[1];

    //cout << "charge=" << charge << endl;
    //cout << "modstr=" << modstr << endl;
    //cout << "tc=" << tc << endl;
    //cout << "tmodstr=" << tmodstr << endl;

    if (charge != tc) {
        return (-2);
    }

    delete si;
    delete tsi;

    int distance = 0;

    int pepLen = peptide.size();

    map<int, int> aaMap;
    for (int i = 0; i < pepLen; i++) {
        if (targetPeptide[i] != peptide[i]) {
            string distanceType = "SNP:" + targetPeptide.substr(i, 1) + "/" + peptide.substr(i, 1);
            if (passValidDistanceType(distanceType)) {
                distance++;
                aaMap[i] = 1;
            } else {
                return (-3);
            }

        } else {
            aaMap[i] = 0;
        }
    }

    //for(map<int,int>::iterator a = aaMap.begin(); a != aaMap.end(); a ++){
    //  cout << a->first <<"\t" << a->second << endl;
    //}

    vector<string> *mods = split(modstr, "/", "/");
    mods->erase(mods->begin());
    map<int, string> modsMap;
    for (vector<string>::iterator j = mods->begin(); j != mods->end(); j++) {
        //cout << *j << endl;
        vector<string> *mi = split(*j, ",", ",");
        int pos = atoi((*mi)[0].c_str());
        string type = (*mi)[2];
        modsMap[pos] = type;
        delete mi;
    }

    //for (map<int,string>::iterator m = modsMap.begin(); m != modsMap.end(); m++){
    //  cout << "modsMap:" << m->first << "\t" << m->second << endl;
    //}

    vector<string> *tmods = split(tmodstr, "/", "/");
    tmods->erase(tmods->begin());
    map<int, string> tmodsMap;
    for (vector<string>::iterator j = tmods->begin(); j != tmods->end(); j++) {
        vector<string> *mi = split(*j, ",", ",");
        int pos = atoi((*mi)[0].c_str());
        string type = (*mi)[2];
        tmodsMap[pos] = type;
        delete mi;
    }
    //for (map<int,string>::iterator m = tmodsMap.begin(); m != tmodsMap.end(); m++){
    //  cout << "tmodsMap:" << m->first << "\t" << m->second << endl;
    //}

    delete mods;
    delete tmods;

    for (int i = -2; i < pepLen; i++) {
        if (i >= 0) {
            if (aaMap[i] == 1) {
                continue;
            }
        }

        string aa("");
        if (i == -1) {
            aa = "Nterm";
        } else if (i == -2) {
            aa = "Cterm";
        } else if (i >= 0) {
            aa = peptide.substr(i, 1);
        }

        //cout << "aa=" << aa << endl;
        /*debug below*/

        if (tmodsMap.find(i) != tmodsMap.end()) {
            if (modsMap.find(i) != modsMap.end()) {
                if (tmodsMap[i] != modsMap[i]) {
                    string distanceType = "PTM:" + aa + "," + tmodsMap[i] + "|" + aa + "," + modsMap[i];
                    if (passValidDistanceType(distanceType)) {
                        distance++;
                    } else {
                        return (-3);
                    }
                }
            } else {
                string distanceType = "PTM:+" + aa + "," + tmodsMap[i];
                if (passValidDistanceType(distanceType)) {
                    distance++;
                } else {
                    return (-3);
                }
            }
        } else {
            if (modsMap.find(i) != modsMap.end()) {
                string distanceType = "PTM:-" + aa + "," + modsMap[i];
                if (passValidDistanceType(distanceType)) {
                    distance++;
                } else {
                    return (-3);
                }
            } else {
                continue;
            }
        }
    }

    return distance;
}

bool SpectraSTSpLibImporter::passValidDistanceType(string type) {
    //cout << "distanceType=" << type << endl;

    if (type.substr(0, 4) == "SNP:") {
        if (m_params.allowableSNP == "ALL") {
            return (true);
        } else if (m_params.allowableSNP == "NONE") {
            return (false);
        } else {
            vector<string> *ti = split(type.substr(4), "/", "/");
            string target = (*ti)[0];
            string orig = (*ti)[1];
            delete ti;

            vector<string> *tsi = split(m_params.allowableSNP, "/", "/");
            for (vector<string>::iterator i = tsi->begin(); i != tsi->end(); i++) {
                vector<string> *ci = split(*i, ",", ",");
                string t = (*ci)[0];
                string o = (*ci)[1];
                delete ci;
                //cout << "t:" << t <<"target:" << target << endl;
                //cout << "o:" << o << "orig:" << orig << endl;
                if (((t == "X") || (t == target)) && ((o == "X") || (o == orig))) {
                    delete tsi;
                    return true;
                }
            }
            delete tsi;
            return false;
        }
    } else if (type.substr(0, 4) == "PTM:") {
        if (m_params.allowableMutatedModifications == "ALL") {
            return (true);
        } else if (m_params.allowableMutatedModifications == "NONE") {
            return (false);
        } else {
            type = type.substr(4);
            vector<string> *msi = split(m_params.allowableMutatedModifications, "/", "/");
            for (vector<string>::iterator i = msi->begin(); i != msi->end(); i++) {
                if ((*i) == type) {
                    return true;
                }
            }
            return false;
        }
    } else {
        return false;
    }
}
