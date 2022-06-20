#include "SpectraSTPepXMLLibImporter.hpp"
#include "SpectraSTReplicates.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "ProgressCount.hpp"
#include "XMLWalker.hpp"
#include <iostream>
#include <sstream>
#include <math.h>
#include <algorithm>


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

/* Class: SpectraSTPepXMLLibImporter
 * 
 * Implements a library importer for the .pepXML file format. Takes a pepXML file,
 * finds all queries with peptide identification over a certain confidence threshold, goes
 * back to the .mzXML file to retrieve the spectrum, and creates a library of these spectra.
 * 
 * Note that this importer will import all replicate spectra of the same peptide.
 * 
 * If multiple instances of the same spectra are present (e.g. as separate queries for different charge state
 * in the same file, or as queries in different searches), the importer will resolve the conflict by
 * only taking the one with the highest probability.
 * 
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTPepXMLLibImporter::SpectraSTPepXMLLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                       SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params),
        m_mzXMLFiles(),
        m_queries(),
        m_numMzXMLOpen(0),
        m_datasetName(""),
        m_numQuery(0),
        m_numQueryWithPeptideProphetProb(0),
        m_numQueryWithPercolatorProb(0),
        m_numQueryWithiProphetProb(0),
        m_numQueryPassedProbCutoff(0),
        m_rtLandmarks(NULL),
        m_numSkipped(0) {

    m_probCutoff = params.minimumProbabilityToInclude;

}

// destructor
SpectraSTPepXMLLibImporter::~SpectraSTPepXMLLibImporter() {

    for (map<string, cRamp *>::iterator i = m_mzXMLFiles.begin(); i != m_mzXMLFiles.end(); i++) {
        if (i->second) {
            delete (i->second);
        }
    }

}

// import - imports all the pepXML files
void SpectraSTPepXMLLibImporter::import() {

    // read the pepXML files
    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        readFromFile(*i);
    }

    // normalize RT
    if (!(m_params.normalizeRTWithLandmarks.empty())) {
        readRTLandmarksFromFile();
    }
    if (m_rtLandmarks) {
        populateRTLandmarks();
    }

    // by now we have all entries to be imported in memory

    m_lib->writePreamble(m_preamble);

    ProgressCount pc(!g_quiet, 1, (int) (m_queries.size()));

    stringstream importProbss;
    importProbss << "Importing all spectra with P>=" << m_probCutoff << " ";
    pc.start(importProbss.str());

    // for efficiency, the entries don't have the spectra yet. Now we load the spectra.
    for (map<string, pair<vector<string>, SpectraSTLibEntry *> >::iterator q = m_queries.begin();
         q != m_queries.end(); q++) {

        pc.increment();

        string query = q->first;
        string path = (q->second.first.size() > 0 ? q->second.first[0] : "");
        string altPath = (q->second.first.size() > 1 ? q->second.first[1] : "");

        vector<SpectraSTLibEntry *> entries;
        entries.push_back(q->second.second);

        if (!(loadSpectrum(query, path, entries, altPath))) {
            // can't load spectrum
            g_log->log("PEPXML IMPORT", "Problem loading spectrum. Skipped query \"" + query + "\".");
            m_numSkipped++;
            delete (entries[0]);
            continue;
        }

        for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {

            SpectraSTLibEntry *entry = (*en);
            SpectraSTPeakList *peakList = entry->getPeakList();

            if (entry->getMods().find("iTRAQ") != string::npos) {
                peakList->removeITRAQReporterPeaks();
            }

            if (entry->getMods().find("TMT") != string::npos) {
                peakList->removeTMTReporterPeaks();
            }

            // for phosphopeptides, get a handle on how good the site assignment is
            if (m_params.evaluatePhosphoSiteAssignment) {
                entry->evaluatePhosphoSiteAssignment();
            }

            // insert the entry into the library
            if (insertOneEntry(entry, "PEPXML")) {
                m_count++;
            } else {
                m_numSkipped++;
            }

            delete (entry);

        }

    }

    stringstream countss;
    countss << "Total of " << m_count << " spectra imported, ";
    countss << m_numSkipped << " spectra skipped.";
    g_log->log("PEPXML IMPORT", countss.str());

    pc.done();

    if ((m_numQueryWithPeptideProphetProb == 0 && m_numQueryWithiProphetProb == 0 &&
         m_numQueryWithPercolatorProb == 0) &&
        m_probCutoff > 0.0000000000001) {
        g_log->warning("PEPXML IMPORT",
                       "Importing a .pep.xml file with no probabilities. PeptideProphet probably needs to be run on .pep.xml first.");
    }

}

// readFromCurFile - reads one pepXML file
void SpectraSTPepXMLLibImporter::readFromFile(string &impFileName) {

    string fullImpFileName(impFileName);
    makeFullPath(fullImpFileName);

    FileName fn;
    parseFileName(fullImpFileName, fn);

    if (m_params.datasetName.empty()) {

        // not specified... use the full path and the base PepXML name as the dataset name
        // to be safe, replace special characters with _
        stringstream dnss;
        string forbidden("~`!@#$%^&*()+\\|:<>,.?/{}'\"");
        string pathPlusBase(fn.path + fn.name);
        for (string::size_type pos = 0; pos < pathPlusBase.length(); pos++) {
            if (forbidden.find(pathPlusBase[pos]) != string::npos) {
                dnss << '_';
            } else {
                dnss << pathPlusBase[pos];
            }
        }

        m_datasetName = dnss.str();

    } else {
        m_datasetName = m_params.datasetName;
    }

    string baseName;
    string searchEngine;
    string database;
    string databaseType;
    string instrumentType;
    string fragType;
    string enzyme;

    ProgressCount pc(!g_quiet, 500, 0);
    string message("Processing \"");
    message += impFileName + "\"";
    pc.start(message);

    map<string, string> r;

    XMLWalker w(impFileName, &r, true);
    if (!w.good()) {
        g_log->error("PEPXML IMPORT", "Cannot open PEPXML file \"" + impFileName + "\" for reading. File skipped.");
        return;
    }

    map<string, int> engines;
    vector<string> stops;

    if (m_params.maximumFDRToInclude < 1.0) {

        r.clear();
        r["min_prob@error_point"] = "";
        r["error@error_point"] = "";
        stops.push_back("<msms_run_summary");
        w.stopAt(stops);

        if (w.walk("<roc_error_data charge=\"all\"", "</roc_error_data>")) {
            if (r["min_prob@error_point"] != "_NOT_FOUND_" && r["error@error_point"] != "_NOT_FOUND_") {
                double actualFDR = 0.0;
                string errorsStr(r["error@error_point"]);
                string minProbsStr(r["min_prob@error_point"]);
                double probCutoffForFDR = getProbCutoffFromFDR(m_params.maximumFDRToInclude, errorsStr, minProbsStr,
                                                               actualFDR);

                stringstream fdrss;
                fdrss << "User-specified global FDR Cutoff is " << m_params.maximumFDRToInclude << ",";
                fdrss << " which corresponds most closely to probability cutoff of " << probCutoffForFDR;
                fdrss << " (Actual Prophet-estimated FDR = " << actualFDR << ")";
                g_log->log("PEPXML IMPORT", fdrss.str());

                if (probCutoffForFDR > m_probCutoff) {
                    m_probCutoff = probCutoffForFDR;
                } else {
                    stringstream fdr2ss;
                    fdr2ss << "Using more stringent probability cutoff of " << m_probCutoff
                           << " specified by -cP option.";
                    g_log->log("PEPXML IMPORT", fdr2ss.str());
                }
            }
        }

        stops.clear();
        w.stopAt(stops);
    }


    do {

        // first, parse out the msms_run_summary element to get the instrument
        r.clear();
        r["base_name"] = "";
        r["msMassAnalyzer"] = "";
        r["name@sample_enzyme"] = "";


        if (w.walk("<msms_run_summary", "<search_summary")) {
            if (r["msMassAnalyzer"] != "_NOT_FOUND_") {
                instrumentType = r["msMassAnalyzer"];
            }
            if (r["base_name"] != "_NOT_FOUND_") {
                baseName = r["base_name"];
            }
            if (r["name@sample_enzyme"] != "_NOT_FOUND_") {
                enzyme = r["name@sample_enzyme"];
            }
        }

        string altPath("");
        string::size_type lastSlashPos = baseName.rfind('/');
        if (lastSlashPos != string::npos) {
            altPath = baseName.substr(0, lastSlashPos + 1);
        }

        if (altPath.empty()) {
            lastSlashPos = baseName.rfind('\\');
            if (lastSlashPos != string::npos) {
                altPath = baseName.substr(0, lastSlashPos + 1);
            }
        }

        // next, parse out the search summary element
        r.clear();
        // r["base_name"] = "";
        r["search_engine"] = "";
        r["local_path@search_database"] = "";
        r["type@search_database"] = "";
        w.setFields(&r, false);
        stops.clear();
        w.stopAt(stops);

        if (w.walk("<search_summary ", "</search_summary>")) {

            // baseName = r["base_name"];
            searchEngine = r["search_engine"];
            database = r["local_path@search_database"];
            databaseType = r["type@search_database"];

        }

        stringstream sess;
        sess << searchEngine << " against \"" << database << "\" (" << databaseType << ")";

        if (engines.find(sess.str()) == engines.end()) {
            engines[sess.str()] = 1;
        }

        // now, parse out the spectrum queries
        r.clear();
        r["spectrum"] = "";
        r["experiment_label"] = "";
        r["ion_mobility"] = "";
        r["assumed_charge"] = "";
        r["retention_time_sec"] = "";
        r["activation_method"] = "";
        r["@search_hit"] = "";
        w.setFields(&r, false);

        stops.push_back("<msms_run_summary");
        w.stopAt(stops);

        while (w.walk("<spectrum_query ", "</spectrum_query>")) {

            if (r["spectrum"] == "_NOT_FOUND_" ||
                r["assumed_charge"] == "_NOT_FOUND_" ||
                r["@search_hit"] == "_NOT_FOUND_") {

                continue;
            }

            string rtstr("");
            if (r["retention_time_sec"] != "_NOT_FOUND_") {
                rtstr = r["retention_time_sec"];
            }

            string imbstr("");
            if (r["ion_mobility"] != "_NOT_FOUND_") {
                imbstr = r["ion_mobility"];
            }

            int charge = atoi((r["assumed_charge"]).c_str());

            string experimentLabel("");
            if (r["experiment_label"] != "_NOT_FOUND_") {
                experimentLabel = r["experiment_label"];
            }

            if (r["activation_method"] != "_NOT_FOUND_") {
                fragType = r["activation_method"];
            }

            m_numQuery++;

            if (processSearchHit(r["@search_hit"], r["spectrum"], charge, rtstr, searchEngine, instrumentType, fragType,
                                 fn.path, altPath, experimentLabel, imbstr)) {
                pc.increment();
            }

        } // END	while (w.walk("<spectrum_query ", "</spectrum_query>"))

    } while (w.current().find("<msms_run_summary") != string::npos);


    // describe what we are doing in the preamble
    stringstream filess;
    filess << "\"" << fullImpFileName << "\", (";
    for (map<string, int>::iterator se = engines.begin(); se != engines.end(); se++) {
        filess << se->first << "; ";
    }
    filess << ") ";
    string desc = m_params.constructDescrStr(filess.str(), ".pepXML");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    stringstream queryss;
    queryss << "Total of " << m_numQuery << " read (" << m_numQueryWithPeptideProphetProb
            << " with PeptideProphet probability; " << m_numQueryWithiProphetProb << " with iProphet probability";
    if (m_numQueryWithPercolatorProb > 0) {
        queryss << "; " << m_numQueryWithPercolatorProb << " with Percolator probability";
    }
    queryss << ").";
    g_log->log("PEPXML IMPORT", queryss.str());

    stringstream countss;
    countss << "Total of " << pc.count() << " spectra (" << m_numQueryPassedProbCutoff
            << " queries) have probability >= " << m_probCutoff << ".";
    g_log->log("PEPXML IMPORT", countss.str());

    pc.done();

}

bool SpectraSTPepXMLLibImporter::processSearchHit(string &searchHit, string &query, int charge, string &rtstr,
                                                  string &searchEngine, string &instrumentType, string &fragType,
                                                  string &path, string &altPath, string &experimentLabel,
                                                  string &imbstr) {

    istringstream *searchHitss = new istringstream(searchHit);
    map<string, string> r;

    r["peptide_prev_aa"] = "";
    r["peptide_next_aa"] = "";
    r["peptide"] = "";
    r["@modification_info"] = "";
    r["probability@peptideprophet_result"] = "";
    r["probability@interprophet_result"] = "";
    r["probability@percolator_result"] = "";
    r["fval"] = "";
    r["massdiff"] = "";

    if (searchEngine == "SEQUEST") {
        r["xcorr"] = "";
        r["deltacn"] = "";
        r["deltacnstar"] = "";
    }
    if (searchEngine == "X! Tandem (k-score)" || searchEngine == "X! Tandem (hrk-score)" ||
        searchEngine == "X! Tandem (native)") {
        r["hyperscore"] = "";
        r["nextscore"] = "";
        r["expect"] = "";
    }
    if (searchEngine == "MASCOT") {
        r["ionscore"] = "";
        r["identityscore"] = "";
        r["homologyscore"] = "";
        r["expect"] = "";
    }
    if (searchEngine == "MyriMatch") {
        r["mvh"] = "";
    }
    if (searchEngine == "InsPecT") {
        r["mqscore"] = "";
        r["expect"] = "";
        r["fscore"] = "";
        r["deltascore"] = "";
    }
    if (searchEngine == "OMSSA") {
        r["expect"] = "";
        r["pvalue"] = "";
    }
    if (searchEngine == "SpectraST") {
        r["dot"] = "";
        r["delta"] = "";
        r["dot_bias"] = "";
    }
    if (searchEngine == "MS-GFDB") {
        r["P-value"] = "";
        r["EFDR"] = "";
    }

    r["protein"] = "";

    XMLWalker w(searchHitss, &r, true, "$");

    w.walk("<search_hit", "<search_hit");

//  delete (searchHitss);

    double prob = 0.0;

    // if multiple probabilities are available, use iProphet over PeptideProphet over Percolator

    double peprob = 0.0;
    if (r["probability@percolator_result"] != "_NOT_FOUND_") {
        peprob = atof((r["probability@percolator_result"]).c_str());
        prob = peprob;
        m_numQueryWithPercolatorProb++;
    }

    double ppprob = 0.0;
    if (r["probability@peptideprophet_result"] != "_NOT_FOUND_") {
        ppprob = atof((r["probability@peptideprophet_result"]).c_str());
        prob = ppprob;
        m_numQueryWithPeptideProphetProb++;
    }

    double ipprob = 0.0;
    if (r["probability@interprophet_result"] != "_NOT_FOUND_") {
        ipprob = atof((r["probability@interprophet_result"]).c_str());
        prob = ipprob;
        m_numQueryWithiProphetProb++;
    }

    double fval = -99999.0;
    if (r["fval"] != "_NOT_FOUND_") {
        fval = atof((r["fval"]).c_str());
    }

    string peptide = r["peptide"];
    if (r["@modification_info"] != "_NOT_FOUND_") {
        addModificationsToPeptide(peptide, r["@modification_info"]);
    }

    // if probability is too low, don't include
    if (prob < m_probCutoff) {
        return (false);
    }

    m_numQueryPassedProbCutoff++;

    // if deltaCn is too small, don't include this
    if (m_params.minimumDeltaCnToInclude > 0.00001 &&
        r["deltacn"] != "_NOT_FOUND_" &&
        atof((r["deltacn"]).c_str()) < m_params.minimumDeltaCnToInclude) {
        g_log->log("PEPXML IMPORT", "DeltaCn too low. Skipped query \"" + query + "\".");
        m_numSkipped++;
        return (false);
    }

    // if massDiff is too big, don't include this
    if (r["massdiff"] != "_NOT_FOUND_" && fabs(atof(r["massdiff"].c_str())) > m_params.maximumMassDiffToInclude) {
        g_log->log("PEPXML IMPORT", "MassDiff too big. Skipped query \"" + query + "\".");
        m_numSkipped++;
        return (false);
    }

    // exclude non peptides
    if (peptide == "XNNPEPTIDEX") {
        g_log->log("PEPXML IMPORT", "Non-peptide ID. Skipped query \"" + query + "\".");
        m_numSkipped++;
        return (false);
    }

    vector<string> prevAAs;
    vector<string> nextAAs;
    if (r["peptide_prev_aa"] != "_NOT_FOUND_") {
        string nextPrevAA("");
        string::size_type divPos = 0;
        while (!((nextPrevAA = nextToken(r["peptide_prev_aa"], divPos, divPos, "$\r\n")).empty())) {
            prevAAs.push_back(nextPrevAA);
            divPos++;
        }
    }
    if (r["peptide_next_aa"] != "_NOT_FOUND_") {
        string nextNextAA("");
        string::size_type divPos = 0;
        while (!((nextNextAA = nextToken(r["peptide_next_aa"], divPos, divPos, "$\r\n")).empty())) {
            nextAAs.push_back(nextNextAA);
            divPos++;
        }
    }

    char prevAA = 'X';
    if (!(prevAAs.empty())) {
        prevAA = prevAAs[0][prevAAs[0].length() - 1];
    }
    char nextAA = 'X';
    if (!(nextAAs.empty())) {
        nextAA = nextAAs[0][0];
    }

    Peptide *pep = createPeptide(peptide, charge, "", query, "PEPXML");

    if (!pep) {
        m_numSkipped++;
        return (false);
    }

    pep->prevAA = prevAA;
    pep->nextAA = nextAA;

    unsigned int highestNTT = pep->NTT();
    if (highestNTT < 2) {

        if (prevAAs.size() > 1 && nextAAs.size() > 1 && prevAAs.size() == nextAAs.size()) {
            // see if there are alternatives with higher NTT
            Peptide *altPep = new Peptide(peptide, charge);
            vector<string>::iterator pr = prevAAs.begin();
            vector<string>::iterator ne = nextAAs.begin();
            for (; pr != prevAAs.end() && ne != nextAAs.end(); pr++, ne++) {
                altPep->prevAA = (*pr)[pr->length() - 1];
                altPep->nextAA = (*ne)[0];
                unsigned int newNTT = 0;
                if ((altPep->prevAA != pep->prevAA || altPep->nextAA != pep->nextAA) &&
                    (newNTT = altPep->NTT()) > highestNTT) {
                    pep->prevAA = altPep->prevAA;
                    pep->nextAA = altPep->nextAA;
                    highestNTT = newNTT;
                    if (highestNTT == 2) break;

                }
            }
            delete (altPep);
        }
    }


    if (m_params.setDeamidatedNXST) {
        setDeamidatedNXST(pep);
    }

    //     cerr << "Importing query " << query << ": " << pep->interactStyleFullWithCharge() << endl;
    SpectraSTLibEntry *entry = new SpectraSTLibEntry(pep, "", "Normal", NULL, fragType);
    entry->getPeakList()->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);

    entry->setOneComment("Spec", "Raw");
    if (pep->NTT() == 2) {
        entry->setOneComment("Pep", "Tryptic");
        entry->setOneComment("NTT", "2");
    } else if (pep->NTT() == 1) {
        entry->setOneComment("Pep", "Semi-tryp");
        entry->setOneComment("NTT", "1");
    } else {
        entry->setOneComment("Pep", "Non-tryp");
        entry->setOneComment("NTT", "0");
    }

    unsigned int nmc = pep->NMC();
    stringstream nmcss;
    nmcss << nmc;
    entry->setOneComment("NMC", nmcss.str());

    if (r["massdiff"] != "_NOT_FOUND_") {
        double massdiff = atof(r["massdiff"].c_str());
        stringstream massdiffss;
        massdiffss.precision(4);
        massdiffss << fixed << massdiff;
        entry->setOneComment("MassDiff", massdiffss.str());
    }

    stringstream naass;
    naass << pep->NAA();
    entry->setOneComment("NAA", naass.str());

    if (!(rtstr.empty())) {
        entry->setOneComment("RetentionTime", rtstr + "," + rtstr + "," + rtstr);
    }

    if (!(imbstr.empty())) {
        entry->setOneComment("IonMobility", imbstr + "," + imbstr + "," + imbstr);
    }

    //      entry->setOneComment("Fullname", pep->fullWithCharge());
    entry->setOneComment("Mods", pep->mspMods());

    stringstream parentss;
    parentss << entry->getPrecursorMz();
    entry->setOneComment("Parent", parentss.str());

    if (instrumentType == "Ion Trap") {
        entry->setOneComment("Inst", "1/it,1,1");
    } else if (instrumentType == "TOF") {
        entry->setOneComment("Inst", "1/tof,1,1");
    } else if (instrumentType == "Q-TOF") {
        entry->setOneComment("Inst", "1/qtof,1,1");
    } else if (instrumentType == "Quadruple Ion Trap") {
        entry->setOneComment("Inst", "1/qit,1,1");
    } else if (!instrumentType.empty()) {
        entry->setOneComment("Inst", "1/" + instrumentType + ",1,1");
    }

    if (!(fragType.empty())) {
        entry->setOneComment("Frag", fragType);
    }

    string proteinss("");
    unsigned int proteinCount = 0;
    if (r["protein"] != "_NOT_FOUND_") {
        string altProteins(r["protein"]);
        string nextProtein("");
        string::size_type divPos = 0;
        while (!((nextProtein = nextToken(altProteins, divPos, divPos, "$\r\n")).empty())) {

            proteinCount++;
            divPos++;

            if (proteinss.empty()) {
                proteinss = nextProtein;
                continue;
            }

            if (nextProtein.compare(0, 5, "DECOY") == 0 || nextProtein.compare(0, 3, "REV") == 0 ||
                nextProtein.compare(0, 3, "rev") == 0) {
                proteinss = proteinss + "/" + nextProtein;
            } else {
                proteinss = nextProtein + "/" + proteinss;
            }

        }
    }

    if (proteinCount > 0) {
        stringstream fullproteinss;
        fullproteinss << proteinCount << "/" << proteinss;
        entry->setOneComment("Protein", fullproteinss.str());
    } else {

        // TODO: FILTERING FROM SGD NAMES -- REMOVE LATER!

//    stringstream fullproteinss;
//    fullproteinss << "0/UNMAPPED_TO_SGD";
//    entry->setOneComment("Protein", fullproteinss.str());
    }

    stringstream searchss;

    char se = searchEngine[0];
    if (searchEngine == "SpectraST") {
        se = 'A';
    }
    if (searchEngine == "X! Tandem (k-score)" || searchEngine == "X! Tandem (hrk-score)") {
        se = 'K';
    }
    if (searchEngine == "MyriMatch") {
        se = 'Y';
    }
    if (searchEngine == "MS-GFDB" || searchEngine == "MS-GF+") {
        se = 'G';
    }


    searchss << "1^" << se << "1:pb=";
    searchss.precision(4);
    searchss << fixed << prob << "/0";

    if (fval > -99998.0) {
        searchss.precision(4);
        searchss << ",fv=" << fixed << fval << "/0";
    }


    if (se == 'S') {
        if (r["xcorr"] != "_NOT_FOUND_") {
            double xcorr = atof((r["xcorr"]).c_str());
            searchss.precision(4);
            searchss << ",xc=" << fixed << xcorr << "/0";

            stringstream xcorrss;
            xcorrss.precision(4);
            xcorrss << fixed << xcorr;
            entry->setOneComment("XCorr", xcorrss.str());
        }

        if (r["deltacn"] != "_NOT_FOUND_") {
            double deltacn = atof((r["deltacn"]).c_str());
            searchss.precision(4);
            searchss << ",dc=" << fixed << deltacn << "/0";

            if (r["deltacnstar"] != "_NOT_FOUND_") {
                double deltacnstar = atof((r["deltacnstar"]).c_str());
                if (deltacnstar > 0.000001 && deltacnstar < 0.999999) {
                    // special case, for Bernd's phospho data
                    deltacn = deltacnstar;
                }
            }

            stringstream deltacnss;
            deltacnss.precision(4);
            deltacnss << fixed << deltacn;
            entry->setOneComment("DeltaCn", deltacnss.str());
        }
    }

    if (se == 'K' || se == 'X') {
        if (r["hyperscore"] != "_NOT_FOUND_") {
            int hyperscore = atoi((r["hyperscore"]).c_str());
            searchss << ",hs=" << fixed << hyperscore << "/0";
        }
        if (r["nextscore"] != "_NOT_FOUND_") {
            int nextscore = atoi((r["nextscore"]).c_str());
            searchss << ",ns=" << fixed << nextscore << "/0";
        }
        if (r["expect"] != "_NOT_FOUND_") {
            double expect = atof((r["expect"]).c_str());
            searchss << ",ex=" << fixed << expect << "/0";
        }
    }

    if (se == 'M') {
        double ionscore = -1.0;
        if (r["ionscore"] != "_NOT_FOUND_") {
            ionscore = atof((r["ionscore"]).c_str());
            searchss << ",sc=" << fixed << ionscore << "/0";
        }
        if (r["identityscore"] != "_NOT_FOUND_") {
            double identityscore = atof((r["identityscore"]).c_str());
            searchss << ",is=" << fixed << identityscore << "/0";
        }
        if (r["homologyscore"] != "_NOT_FOUND_") {
            double homologyscore = atof((r["homologyscore"]).c_str());
            searchss << ",hs=" << fixed << homologyscore << "/0";
            if (ionscore >= 0.0) searchss << ",sr=" << fixed << ionscore - homologyscore << "/0";
        }
        if (r["expect"] != "_NOT_FOUND_") {
            double expect = atof((r["expect"]).c_str());
            searchss << ",ex=" << fixed << expect << "/0";
        }
    }

    if (se == 'Y') {
        double mvh = -1.0;
        if (r["mvh"] != "_NOT_FOUND_") {
            mvh = atof(r["mvh"].c_str());
            searchss << ",mv=" << fixed << mvh << "/0";
        }
    }

    if (se == 'I') {
        if (r["mqscore"] != "_NOT_FOUND_") {
            double mqscore = atoi((r["mqscore"]).c_str());
            searchss << ",mq=" << fixed << mqscore << "/0";
        }
        if (r["fscore"] != "_NOT_FOUND_") {
            double fscore = atoi((r["fscore"]).c_str());
            searchss << ",fs=" << fixed << fscore << "/0";
        }
        if (r["expect"] != "_NOT_FOUND_") {
            double expect = atof((r["expect"]).c_str());
            searchss << ",ex=" << fixed << expect << "/0";
        }
        if (r["deltascore"] != "_NOT_FOUND_") {
            double deltascore = atof((r["deltascore"]).c_str());
            searchss << ",ds=" << fixed << deltascore << "/0";
        }
    }

    if (se == 'O') {
        if (r["pvalue"] != "_NOT_FOUND_") {
            double pvalue = atoi((r["pvalue"]).c_str());
            searchss << ",pv=" << fixed << pvalue << "/0";
        }
        if (r["expect"] != "_NOT_FOUND_") {
            double expect = atof((r["expect"]).c_str());
            searchss << ",ex=" << fixed << expect << "/0";
        }
    }

    if (se == 'A') {
        if (r["dot"] != "_NOT_FOUND_") {
            double dot = atof((r["dot"]).c_str());
            searchss << ",do=" << fixed << dot << "/0";
        }
        if (r["delta"] != "_NOT_FOUND_") {
            double delta = atof((r["delta"]).c_str());
            searchss << ",dd=" << fixed << delta << "/0";
        }
        if (r["dot_bias"] != "_NOT_FOUND_") {
            double dotBias = atof((r["dot_bias"]).c_str());
            searchss << ",db=" << fixed << dotBias << "/0";
        }
    }

    if (se == 'G') {
        if (r["P-value"] != "_NOT_FOUND_") {
            double pvalue = atoi((r["pvalue"]).c_str());
            searchss << ",pv=" << fixed << pvalue << "/0";
        }
        if (r["EFDR"] != "_NOT_FOUND_") {
            double efdr = atof((r["expect"]).c_str());
            searchss << ",ef=" << fixed << efdr << "/0";
        }
    }

    // HAVE TO ADD SUPPORT FOR OTHER SEARCH ENGINES LATER

    entry->setOneComment("Se", searchss.str());


    string sample(m_datasetName);

    if (!experimentLabel.empty() && m_params.datasetName.empty()) {
        sample = experimentLabel;
    }

    if (m_params.addMzXMLFileToDatasetName) {
        string::size_type dummy = 0;
        string prefix = nextToken(query, 0, dummy, ". \t\r\n");
        // also tag on mzXML file name -- i.e. substring before first dot in query (e.g. to record which SCX/IEF fraction this is from)
        entry->setOneComment("Sample", "1/" + sample + "_" + prefix + ",1,1");
    } else {
        entry->setOneComment("Sample", "1/" + sample + ",1,1");
    }

    entry->setOneComment("Nreps", "1/1");

    string::size_type dotPos = query.rfind('.');
    query = query.substr(0, dotPos);

    entry->setOneComment("RawSpectrum", query); // RawSpectrum field excludes the .<charge> for consistency

    stringstream probss;
    probss.precision(4);
    probss << fixed << prob;
    entry->setOneComment("Prob", probss.str());

    // check to see if the query is seen before

    // take off the charge

//  query = path + query;

    // feature to extract spectra in a bracket (multiple MS2 of the same ion)
    if (m_params.bracketSpectra) {
        extractBracket(query, path, altPath);
    }

    map<string, pair<vector<string>, SpectraSTLibEntry *> >::iterator found = m_queries.find(query);
    if (found == m_queries.end()) {
        // query not seen before. add.
        // add to a structure in memory for now
        pair<vector<string>, SpectraSTLibEntry *> p;
        p.first.push_back(path);
        p.first.push_back(altPath);
        p.second = entry;
        m_queries[query] = p;

        return (true);

    } else {

        // query seen before. resolve.

        SpectraSTLibEntry *foundEntry = found->second.second;
        // note the probability of the existing one
        double oldProb = foundEntry->getProb();

        // check if the peptide ion is the same as the existing one
        Peptide *oldpep = foundEntry->getPeptidePtr();
        oldpep->considerIsobaricAASame = true;
        if (*oldpep == *pep) {
            // same peptide! add to Se field
            // cerr << "Same query! " << query << " ID: " << pep->interactStyleWithCharge() << " and " << oldpep->interactStyleWithCharge() << endl;

            map<char, map<string, pair<double, double> > *> seqs;
            foundEntry->getSeqInfo(seqs);
            map<char, map<string, pair<double, double> > *>::iterator foundse = seqs.find(se);
            if (foundse == seqs.end()) {
                // new search engine
                entry->getSeqInfo(seqs);
                /*
                map<string, pair<double, double> >* scs = new map<string, pair<double, double> >;
                pair<double, double> pbPair;
                pbPair.first = prob;
                pbPair.second = 0.0;
                (*scs)["pb"] = pbPair;
                pair<double, double> nrPair;
                nrPair.first = 1.0;
                nrPair.second = 0.0;
                (*scs)["nr"] = nrPair;
                seqs[se] = scs;
                */
            } else {
                // already have that search engine
                map<string, pair<double, double> >::iterator foundpb = foundse->second->find("pb");
                if (foundpb != foundse->second->end()) {
                    if (prob > foundpb->second.first) {
                        seqs.erase(se);
                        entry->getSeqInfo(seqs);
//            foundpb->second.first = prob;
                    }
                }
            }

            entry->setSeqInfo(seqs);

            // Prob should reflect the max probability of all sequence search results
            // since this new entry has lower probability than the existing one, set it back to the existing one
            if (oldProb > prob) {
                stringstream oldProbss;
                oldProbss.precision(4);
                oldProbss << fixed << showpoint << oldProb;
                entry->setOneComment("Prob", oldProbss.str());
            }

            // replace existing one with new entry (which already incorporated the existing one)
            delete (foundEntry);
            found->second.first.clear();
            found->second.first.push_back(path);
            found->second.first.push_back(altPath);
            found->second.second = entry;

        } else {

            /* TESTING -- Uncomment to keep both */

            /*
            cerr << "Same query! " << query << " ID: " << pep->interactStyleWithCharge() << " and " << oldpep->interactStyleWithCharge() << endl;

            found->second.second->setOneComment("Mixie", pep->interactStyleWithCharge());
            entry->setOneComment("Mixie", oldpep->interactStyleWithCharge());
            pair<vector<string>, SpectraSTLibEntry*> p;
            p.first.push_back(path);
            p.first.push_back(altPath);
            p.second = entry;

            // hack the query such that it won't collide with the old one
            string::size_type lastDotPos = query.rfind('.');   // .<lastScanNum>  (<charge> has been taken off already)
            if (lastDotPos == string::npos) {
              g_log->error("PEPXML IMPORT", "Illegal query name \"" + query + "\". Scan not imported.");
              return (false);
            }
            string::size_type secondLastDotPos = query.rfind('.', lastDotPos - 1);  // .<firstScanNum>.<lastScanNum>
            if (secondLastDotPos == string::npos) {
              g_log->error("PEPXML IMPORT", "Illegal query name \"" + query + "\". Scan not imported.");
              return (false);
            }

            unsigned int decr = 0;
            string newQuery("");
            do {
          decr++;
              stringstream newQueryss;
          newQueryss << query.substr(0, lastDotPos) << '.';
          newQueryss.width(5);
          newQueryss.fill(0);
          newQueryss << right << atoi(query.substr(lastDotPos + 1).c_str()) - decr;
          newQuery = newQueryss.str();
            } while (m_queries.find(newQuery) != m_queries.end());

            m_queries[newQuery] = p;

            return (true);

            */
            /* END TESTING -- Uncomment to keep both */

            /* TESTING -- Uncomment to delete both */
            /*
            g_log->log("PEPXML IMPORT", "Conflicting IDs (" + pep->interactStyleWithCharge() + " and " + oldpep->interactStyleWithCharge() + "). Skipped query \"" + query + "\".");
            m_numSkipped++; // this is for the original ID, also being deleted

            delete (found->second.second);
            m_queries.erase(found);
            delete (entry);
            return (false);
            */
            /* END TESTING -- Uncomment to delete both */

            // not the same - just pick the higher-probability one

            if (oldProb < prob) {
                // this one is better, replace.
                delete (found->second.second);
                found->second.first.clear();
                found->second.first.push_back(path);
                found->second.first.push_back(altPath);
                found->second.second = entry;

            } else {
                // the old one is better, ignore this one.
                delete entry;
            }

        } // END if (oldpep == pep)

        return (false);

    } // END if (found == m_queries.end())

}


// loadSpectrum - go to the mzXML file to load the spectrum
bool SpectraSTPepXMLLibImporter::loadSpectrum(string query, string path, vector<SpectraSTLibEntry *> &entries,
                                              string altPath) {

    if (entries.empty()) return (false);

    SpectraSTLibEntry *entry = entries[0];
    SpectraSTPeakList *peakList = entry->getPeakList();

    int firstScanNum = 0;
    int lastScanNum = 0;
    string baseName;

    if (!parseQuery(query, baseName, firstScanNum, lastScanNum)) {
        g_log->error("PEPXML IMPORT", "Illegal query name \"" + query + "\". Scan not imported.");
        return (false);
    }

    cRamp *cramp = openCRamp(baseName, path, altPath);

    if (!cramp) return (false);

    string fullFileName = path + baseName + ".mzXML";

    int scanNum = firstScanNum;

    // check if the searched scan is a merged one via the CombIon feature of Thermo
    if (firstScanNum != lastScanNum) {

        entry->setOneComment("OrigQuery", query);

        vector<SpectraSTLibEntry *> bracket;

        double firstPrecursorMz = 0.0;
        SpectraSTLibEntry *firstEntry = new SpectraSTLibEntry(*entry);
        if (readOneScan(cramp, scanNum, firstEntry, baseName, fullFileName)) {
            firstPrecursorMz = firstEntry->getPrecursorMz();
            stringstream rawss;
            rawss << baseName << "." << scanNum << "." << scanNum;
            firstEntry->setOneComment("RawSpectrum", rawss.str());
            bracket.push_back(firstEntry);
        } else {
            // even first one failed, just bail
            delete (firstEntry);
            return (false);
        }

        for (int scanNum = firstScanNum + 1; scanNum <= lastScanNum; scanNum++) {

            SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*entry);
            if (readOneScan(cramp, scanNum, newEntry, baseName, fullFileName, true, firstPrecursorMz - 1.0,
                            firstPrecursorMz + 1.0)) {
                stringstream rawss;
                rawss << baseName << "." << scanNum << "." << scanNum;
                newEntry->setOneComment("RawSpectrum", rawss.str());
                bracket.push_back(newEntry);
            } else {
                delete (newEntry);
            }
        }

        if (bracket.empty()) {
            // should not happen!
            return (false);

        } else if (bracket.size() == 1) {
            *entry = *(bracket[0]);
            delete (bracket[0]);

        } else {

            if (m_params.mergeBracket) {

                SpectraSTReplicates reps(bracket, m_params, NULL);
                reps.setRecordRawSpectra(
                        true); // overrides m_params.recordRawSpectra. In this case of CombIon, always record raw spectra
                SpectraSTLibEntry *consensus = reps.makeConsensusSpectrum();
                *entry = *consensus; // deep copy

                // delete the bracket
                for (vector<SpectraSTLibEntry *>::iterator den = bracket.begin(); den != bracket.end(); den++) {
                    if (*den) delete (*den);
                }

            } else {

                // keep all spectra in the bracket
                delete (entry);
                entries = bracket; // this replaces the whole entries vector

            }


        }

        return (true);

    } else {

        // this is the usual operation, when firstScanNum and lastScanNum are identical
        return (readOneScan(cramp, scanNum, entry, baseName, fullFileName));
    }

    return (true);

}

bool SpectraSTPepXMLLibImporter::readOneScan(cRamp *cramp, int scanNum, SpectraSTLibEntry *entry, string &baseName,
                                             string &fullFileName, bool silent, double minMz, double maxMz) {

    rampScanInfo *scanInfo = cramp->getScanHeaderInfo(scanNum);
    if (!scanInfo || scanInfo->m_data.acquisitionNum != scanNum || scanInfo->m_data.msLevel == 1) {
        // bad. not found!
        if (!silent) {
            stringstream errss;
            errss << "Cannot find MS2+ scan #" << scanNum << " in file \"" << fullFileName << "\". Scan not imported.";
            g_log->error("PEPXML IMPORT", errss.str());
        }
        if (scanInfo) delete (scanInfo);
        return (false);
    }

    double precursorMz = scanInfo->m_data.precursorMZ;

    if (precursorMz < minMz || precursorMz > maxMz) {
        if (!silent) {
            stringstream errss;
            errss << "Precursor M/Z: " << precursorMz << " outside range [" << minMz << ", " << maxMz
                  << "]. Scan not imported.";
            g_log->error("PEPXML IMPORT", errss.str());
        }
        delete (scanInfo);
        return (false);
    }

    // now we can read the peaks
    rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);
    if (!peaks) {
        if (!silent) {
            stringstream errss;
            errss << "Cannot read peaks for scan #" << scanNum << " in file \"" << fullFileName
                  << "\". Scan not imported.";
            g_log->error("PEPXML IMPORT", errss.str());
        }
        delete (scanInfo);
        return (false);
    }

    int peakCount = peaks->getPeakCount();
    int precursorCharge = scanInfo->m_data.precursorCharge;
    double precursorIntensity = scanInfo->m_data.precursorIntensity;
    double totIonCurrent = scanInfo->m_data.totIonCurrent;
    double retentionTime = scanInfo->m_data.retentionTime;
    string fragType(scanInfo->m_data.activationMethod);
    double collisionEnergy = scanInfo->m_data.collisionEnergy;

    if (!(fragType.empty())) {
        entry->setFragType(fragType);
    }

    // will overwrite the retention time from pepXML file with that from the mzXML file
    stringstream rtss;
    rtss.precision(1);
    rtss << fixed << retentionTime << ',' << retentionTime << ',' << retentionTime;
    entry->setOneComment("RetentionTime", rtss.str());

    if (m_rtLandmarks) {
        double iRT = normalizeRT(baseName, retentionTime);
        if (iRT >= -999999999.0) { // uninitialized value is -1000000000.0
            stringstream irtss;
            irtss.precision(1);
            irtss << fixed << iRT << ',' << iRT << ',' << iRT;
            entry->setOneComment("iRT", irtss.str());
        }

        // cerr << baseName << "\t" << entry->getPeptidePtr()->interactStyleWithCharge() << "\t" << retentionTime << "\t" << iRT << endl;

    }

    stringstream ticss;
    ticss.precision(2);
    ticss << totIonCurrent;
    entry->setOneComment("TotalIonCurrent", ticss.str());

    stringstream pintss;
    pintss.precision(2);
    pintss << precursorIntensity;
    entry->setOneComment("PrecursorIntensity", pintss.str());

    if (collisionEnergy >= 0.0001 && collisionEnergy <= 1000.0) {
        stringstream cess;
        cess.precision(1);
        cess << fixed << collisionEnergy;
        entry->setOneComment("CollisionEnergy", cess.str());
    }

    if (precursorCharge < 1) precursorCharge = 0;

//  cout << "inserting peaks " << peakCount << endl;
    // create the peak list and read the peaks one-by-one
    for (int j = 0; j < peakCount; j++) {
        double mz = peaks->getPeak(j)->mz;
        float intensity = (float) (peaks->getPeak(j)->intensity);

        if (intensity > 0.1) {
            entry->getPeakList()->insert(mz, intensity, "", "");
        }
    }

    if (precursorIntensity > 0.0) {
        entry->getPeakList()->setPrecursorIntensity(precursorIntensity);
    }

    stringstream maxss;
    maxss.precision(2);
    maxss << entry->getPeakList()->getOrigMaxIntensity();
    entry->setOneComment("OrigMaxIntensity", maxss.str());

    delete scanInfo;
    delete peaks;

    return (true);

}


// setStaticMods (Deprecated) - based on the parsed out static mods in the pepXML file, tell the Peptide class how
// to parse peptide strings. (i.e. if a static mod of CAM is specified on a cysteine, then all 'C' will be treated
// as CAM-cysteine.)
void SpectraSTPepXMLLibImporter::setStaticMods(string aas, string masses, string variables) {

    string::size_type c1 = 0;
    string::size_type c2 = 0;
    string::size_type c3 = 0;


    while (c1 < aas.length() &&
           c2 < masses.length() &&
           c3 < variables.length()) {

        string aa = nextToken(aas, c1, c1, "|/r/n");
        double mass = atof((nextToken(masses, c2, c2, "|/r/n")).c_str());
        string variable = nextToken(variables, c3, c3, "|/r/n");

        if (variable == "N") {
            if (aa == "C" && (int) (mass + 0.5) == 161) {
                Peptide::addModTokenToTables("C", "Carboxymethyl");
                g_log->log("PEPXML IMPORT", "Detected static modification C => C[161] (Carboxymethyl)");
            } else if (aa == "C" && (int) (mass + 0.5) == 160) {
                // CAM
                Peptide::addModTokenToTables("C", "Carbamidomethyl");
                g_log->log("PEPXML IMPORT", "Detected static modification C => C[160] (Carbamidomethyl)");

            } else if (aa == "C" && (int) (mass + 0.5) == 330) {
                // ICAT_cl
                Peptide::addModTokenToTables("C", "ICAT-C");
                g_log->log("PEPXML IMPORT", "Detected static modification C => C[330] (ICAT-C)");

            } else if (aa == "C" && ((int) (mass + 0.5) >= 545 && (int) (mass + 0.5) <= 546)) {
                // ICAT_uc
                Peptide::addModTokenToTables("C", "ICAT-D");
                g_log->log("PEPXML IMPORT", "Detected static modification C => C[546] (ICAT-D)");

            } else if (aa == "C" && ((int) (mass + 0.5) == 518)) {
                // PEO -- Hui's special alkylation agent that nobody else uses
                Peptide::addModTokenToTables("C", "PEO-Iodoacetyl-LC-Biotin");
                g_log->log("PEPXML IMPORT", "Detected static modification C => C[518] (PEO-Iodoacetyl-LC-Biotin)");

            } else {
                // unknown, don't bother at this point
            }
        }
        c1++;
        c2++;
        c3++;
    }
}

// addModificationsToPeptide - takes the <modification_info> element in the pepXML file, parses it
// and adds the modifications to the peptide string
void SpectraSTPepXMLLibImporter::addModificationsToPeptide(string &peptide, string &modInfo) {

    istringstream *modss = new istringstream(modInfo);
    map<string, string> modr;
    modr["mod_nterm_mass"] = "";
    modr["mod_cterm_mass"] = "";
    modr["position"] = "";
    modr["mass"] = "";
    XMLWalker modw(modss, &modr, true);

    modw.walk("<modification_info", "BLAH"); // won't find "BLAH" ever, so will walk until the end

    stringstream newpeptide;

    if (modr["mod_nterm_mass"] != "_NOT_FOUND_") {
        int ntermMass = (int) (atof(modr["mod_nterm_mass"].c_str()) + 0.5);
        newpeptide << "n[" << ntermMass << "]";
    }

    if (modr["position"] == "_NOT_FOUND_") {
        newpeptide << peptide;
    } else {
        string positions = modr["position"];
        string masses = modr["mass"];

        map<int, string> modaa;
        string::size_type c1 = 0;
        string::size_type c2 = 0;
        while (c1 < positions.length() &&
               c2 < masses.length()) {

            int position = atoi(nextToken(positions, c1, c1, "|\r\n").c_str()) - 1;
            string massStr = nextToken(masses, c2, c2, "|\r\n");
            modaa[position] = massStr;
            c1++;
            c2++;
        }

        for (string::size_type i = 0; i < peptide.length(); i++) {
            map<int, string>::iterator found = modaa.find((int) i);
            if (found != modaa.end()) {
                string massStr = found->second;
                if (massStr.length() > 0 && (!isdigit(massStr[0]))) {
                    // Not a number. Assume it is a SpectraST-recognized mod tag
                    newpeptide << peptide[i] << '[' << massStr << ']';
                } else {
                    int nominalMass = (int) (atof(massStr.c_str()) + 0.5);
                    newpeptide << peptide[i] << '[' << nominalMass << ']';
                }
            } else {
                newpeptide << peptide[i];
            }
        }
    }

    if (modr["mod_cterm_mass"] != "_NOT_FOUND_") {
        int ctermMass = (int) (atof(modr["mod_cterm_mass"].c_str()) + 0.5);
        newpeptide << "c[" << ctermMass << "]";
    }

    peptide = newpeptide.str();
}

bool SpectraSTPepXMLLibImporter::parseQuery(string &query, string &baseName, int &firstScanNum, int &lastScanNum) {


    string::size_type lastDotPos = query.rfind('.');   // .<lastScanNum>  (<charge> has been taken off already)
    if (lastDotPos == string::npos) {
        return (false);
    }
    string::size_type secondLastDotPos = query.rfind('.', lastDotPos - 1);  // .<firstScanNum>.<lastScanNum>
    if (secondLastDotPos == string::npos) {
        return (false);
    }

    firstScanNum = atoi((query.substr(secondLastDotPos + 1, lastDotPos - secondLastDotPos - 1)).c_str());
    lastScanNum = atoi((query.substr(lastDotPos + 1)).c_str());

    baseName = query.substr(0, secondLastDotPos);

    return (true);
}


cRamp *SpectraSTPepXMLLibImporter::openCRamp(string &baseName, string &path, string &altPath) {

    // even though the scans may not be in this fullFileName, still use this as the key in the map
    string fullFileName = path + baseName + ".mzXML";

    // try to see if the file is already opened
    map<string, cRamp *>::iterator found = m_mzXMLFiles.find(fullFileName);
    cRamp *cramp = NULL;
    if (found == m_mzXMLFiles.end()) {
        // can't find this mzXML file
        if (m_numMzXMLOpen >= 3) {
            // already at max number of open files... need to close one
            // pick the first file in m_mzXMLFiles to close -- note that this will not cause
            // incessant opening and closing of that first file, because m_queries is sorted
            // by query name (i.e. <mzXMLFileName>.<scanNum>.<scanNum>), and this loadSpectrum
            // is called in the order of query name. Therefore, it should be done with the first
            // file before it starts on the second, and so forth.
            map<string, cRamp *>::iterator dead = m_mzXMLFiles.begin();
            delete (dead->second);
            m_mzXMLFiles.erase(dead);
            m_numMzXMLOpen--;
        }

        string triedStr(baseName + ".mzXML");
        string tryFileName(fullFileName);
        cramp = new cRamp(tryFileName.c_str());

        const char **rampSupportedTypes = rampListSupportedFileTypes();
        unsigned int rampSupportedTypesIndex = 0;

        // if the mzXML file is not there, look for other RAMP-supported formats (.mzData, .mzML, .mzXML.gz, etc)
        while (!cramp->OK()) {
            const char *nextTypeCPtr = rampSupportedTypes[rampSupportedTypesIndex++];
            if (!nextTypeCPtr) break;
            string nextType(nextTypeCPtr);
            if (nextType == ".mzXML") continue; // already tried this one
            tryFileName = path + baseName + nextType;
            triedStr += "|" + nextType;
            delete (cramp);
            cramp = new cRamp(tryFileName.c_str());
        }
        triedStr += " in " + path;

        // if still cannot find the spectrum file, try find it in altPath
        rampSupportedTypesIndex = 0;
        if (!(altPath.empty()) && altPath != path) {
            while (!cramp->OK()) {
                const char *nextTypeCPtr = rampSupportedTypes[rampSupportedTypesIndex++];
                if (!nextTypeCPtr) break;
                string nextType(nextTypeCPtr);
                tryFileName = altPath + baseName + nextType;
                delete (cramp);
                cramp = new cRamp(tryFileName.c_str());
                if (cramp->OK()) {
                    cerr << "INFO:  Reading file: " << tryFileName << endl;
                }
            }
            triedStr += " or in " + altPath;
        }


        if (!cramp->OK()) {
            g_log->error("PEPXML IMPORT",
                         "Cannot open file \"" + triedStr + "\". No scan from this file will be imported.");
            delete (cramp);
            m_mzXMLFiles[fullFileName] = NULL; // note this so that future attempts to read scan from this file will die silently
            return (NULL);
        } else {
            m_numMzXMLOpen++;
            m_mzXMLFiles[fullFileName] = cramp;
            return (cramp);
        }

    } else {
        if (found->second) {
            cramp = found->second;
            return (cramp);
        } else {
            // tried to open this file previously and failed
            return (NULL);
        }
    }

}


// extractBracket - given a query, goes into the corresponding mzXML file, find the scan, and look in neighboring
// MS2 scans to see if the precursor m/z is the same. If so, they are considered part of the bracket. The query
// string will be changed to <mzXML>.<first scan in the bracket>.<last scan in the bracket>
void SpectraSTPepXMLLibImporter::extractBracket(string &query, string &path, string &altPath) {

    int firstScanNum = 0;
    int lastScanNum = 0;
    string baseName;

    if (!parseQuery(query, baseName, firstScanNum, lastScanNum)) {
        return;
    }

    if (firstScanNum != lastScanNum) return; // CombIon mode on. This gets complicated, don't extract bracket

    cRamp *cramp = openCRamp(baseName, path, altPath);

    if (!cramp) return;

    // get precursor m/z of this scan
    rampScanInfo *scanInfo = cramp->getScanHeaderInfo(firstScanNum);
    if (!scanInfo || scanInfo->m_data.acquisitionNum != firstScanNum || scanInfo->m_data.msLevel == 1) {
        // bad. not found!
        if (scanInfo) delete (scanInfo);
        return;
    }
    double precursorMz = scanInfo->m_data.precursorMZ;
    delete (scanInfo);

    // look before
    int scanNum = firstScanNum;
    int bracketStart = firstScanNum;

    while (scanNum > 1) {
        scanNum--;
        scanInfo = cramp->getScanHeaderInfo(scanNum);
        if (!scanInfo || scanInfo->m_data.acquisitionNum != scanNum) {
            if (scanInfo) delete (scanInfo);
            continue;
        }
        if (scanInfo->m_data.msLevel == 1) break; // MS1 scan, stop searching

        double prevPrecursorMz = scanInfo->m_data.precursorMZ;
        delete (scanInfo);

        if (fabs(prevPrecursorMz - precursorMz) < 0.01) {
            // same precursor m/z! it belongs to the bracket
            bracketStart = scanNum;
        } else {
            // different precursor m/z. stop searching
            break;
        }
    }

    scanNum = firstScanNum;
    int bracketEnd = firstScanNum;
    // look after
    while (scanNum < cramp->getLastScan()) {
        scanNum++;
        scanInfo = cramp->getScanHeaderInfo(scanNum);
        if (!scanInfo || scanInfo->m_data.acquisitionNum != scanNum) {
            if (scanInfo) delete (scanInfo);
            continue;
        }
        if (scanInfo->m_data.msLevel == 1) break; // MS1 scan, stop searching

        double nextPrecursorMz = scanInfo->m_data.precursorMZ;
        delete (scanInfo);

        if (fabs(nextPrecursorMz - precursorMz) < 0.01) {
            // same precursor m/z! it belongs to the bracket
            bracketEnd = scanNum;
        } else {
            // different precursor m/z. stop searching
            break;
        }
    }

    if (bracketStart == bracketEnd) {
        // no bracket, just return without changing the query
        return;
    }

    // rewrite query
    stringstream newQueryss;
    newQueryss << baseName << '.';
    newQueryss.width(5);
    newQueryss.fill('0');
    newQueryss << right << bracketStart;
    newQueryss << '.';
    newQueryss.width(5);
    newQueryss.fill('0');
    newQueryss << right << bracketEnd;

    query = newQueryss.str();

}

double SpectraSTPepXMLLibImporter::getProbCutoffFromFDR(double desiredFDR, string &errorsStr, string &minProbsStr,
                                                        double &actualFDR) {

    map<double, double, std::greater<double> > cutoffs;
    string::size_type c1 = 0;
    string::size_type c2 = 0;

    while (c1 < errorsStr.length() &&
           c2 < minProbsStr.length()) {

        double error = atof(nextToken(errorsStr, c1, c1, "|/r/n").c_str());
        double minProb = atof(nextToken(minProbsStr, c2, c2, "|/r/n").c_str());
        cutoffs[error] = minProb;
        c1++;
        c2++;
    }

    map<double, double>::iterator found = cutoffs.upper_bound(desiredFDR + 0.00000001);

    if (found != cutoffs.end()) {
        actualFDR = found->first;
        return (found->second);
    } else {
        actualFDR = 0.0;
        return (1.0);
    }


}


double SpectraSTPepXMLLibImporter::normalizeRT(string &baseName, double rRT) {

    double iRT = -1000000000.0; // uninitialized value -- hopefully this won't happen in reality

    if (m_params.normalizeRTLinearRegression) {

        map<string, pair<double, double> >::iterator foundCoeff = m_normalizeRTCoeff.find(baseName);

        if (foundCoeff != m_normalizeRTCoeff.end()) {
            iRT = calculateIRTLinearRegression(rRT, foundCoeff->second);
        }

    } else {

        map<string, map<double, double> >::iterator foundPivots = m_normalizeRTPivots.find(baseName);

        if (foundPivots != m_normalizeRTPivots.end()) {
            iRT = calculateIRTLinearInterpolation(rRT, foundPivots->second);
        }

    }

    return (iRT);

}


double SpectraSTPepXMLLibImporter::calculateIRTLinearRegression(double rRT, pair<double, double> &coeff) {

    return ((rRT - coeff.second) / coeff.first);

}

double SpectraSTPepXMLLibImporter::calculateIRTLinearInterpolation(double rRT, map<double, double> &pivots) {

    map<double, double>::iterator upper = pivots.upper_bound(rRT);

    double rRT1 = 0.0;
    double iRT1 = 0.0;
    double rRT2 = 0.0;
    double iRT2 = 0.0;

    if (upper == pivots.begin()) {
        // rRT is before the first pivot: extrapolate line between first and second pivot

        rRT1 = upper->first;
        iRT1 = upper->second;

        upper++;
        rRT2 = upper->first;
        iRT2 = upper->second;

    } else if (upper == pivots.end()) {
        // rRT is after the last pivot: extrapolate line between last and next-to-last

        upper--;
        rRT2 = upper->first;
        iRT2 = upper->second;

        upper--;
        rRT1 = upper->first;
        iRT1 = upper->second;

    } else {

        rRT2 = upper->first;
        iRT2 = upper->second;

        upper--;
        rRT1 = upper->first;
        iRT1 = upper->second;

    }

    if (fabs(iRT2 - iRT1) < 0.00001) {
        return (iRT1);
    }

    double slope = (rRT2 - rRT1) / (iRT2 - iRT1);

    if (slope < 0.00001) {
        return (iRT1);
    }

    return (iRT1 - (rRT1 - rRT) / slope);

}

void SpectraSTPepXMLLibImporter::readRTLandmarksFromFile() {

    // read standard peptides
    ifstream fin;
    if (!myFileOpen(fin, m_params.normalizeRTWithLandmarks)) {
        g_log->error("PEPXML IMPORT", "Cannot read landmark table. No RT normalization will be performed.");
        return;
    }

    if (m_rtLandmarks) delete m_rtLandmarks;
    m_rtLandmarks = new map<string, pair<double, vector<double> > >;

    string line("");
    stringstream pepss;
    pepss << "Read RT landmarks from file \"" << m_params.normalizeRTWithLandmarks << "\" : ";

    while (nextLine(fin, line, "_EOF_", "")) {

        string::size_type pos = 0;
        string pepStr = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
        if (pepStr.empty() || pepStr[0] == '#') {
            continue;
        }

        string iRTStr = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");

        if (!iRTStr.empty() && (iRTStr[0] == '-' || iRTStr[0] == '.' || (iRTStr[0] >= '0' && iRTStr[0] <= '9'))) {
            double iRT = atof(iRTStr.c_str());
            pepss << pepStr << " (";
            pepss.precision(1);
            pepss << fixed << iRT << ") | ";

            (*m_rtLandmarks)[pepStr].first = iRT;
        }


    }

    double sumIRT = 0.0;
    double sumIRTSq = 0.0;
    for (map<string, pair<double, vector<double> > >::iterator i = m_rtLandmarks->begin();
         i != m_rtLandmarks->end(); i++) {
        sumIRT += i->second.first;
        sumIRTSq += i->second.first * i->second.first;
    }
    int count = m_rtLandmarks->size();
    sumIRT /= (double) (count);
    sumIRTSq /= (double) (count);
    double variance = sumIRTSq - sumIRT * sumIRT;

    if (count < 2 || variance < 0.0001) {
        pepss << "\b\b. Too few landmarks with distinct specified iRT values to perform RT normalization.";
        delete (m_rtLandmarks);
        m_rtLandmarks = NULL;
        g_log->error("PEPXML IMPORT", pepss.str());
    } else {
        pepss << "\b\b.";
        g_log->log("PEPXML IMPORT", pepss.str());
    }

}

void SpectraSTPepXMLLibImporter::populateRTLandmarks() {

    if (!m_rtLandmarks) return;

    string curBaseName("");

    for (map<string, pair<vector<string>, SpectraSTLibEntry *> >::iterator q = m_queries.begin();
         q != m_queries.end(); q++) {

        string query = q->first;
        SpectraSTLibEntry *entry = q->second.second;
        Peptide *pep = entry->getPeptidePtr();

        // cerr << query << ": " << pep->interactStyle() << endl;

        if (!pep) continue;

        string pepStr = pep->interactStyle();
        map<string, pair<double, vector<double> > >::iterator found = m_rtLandmarks->find(pepStr);

        if (found == m_rtLandmarks->end()) continue;

        // cerr << "  IS LANDMARK!" << endl;

        string baseName("");
        int firstScanNum = 0;
        int lastScanNum = 0;

        if (!parseQuery(query, baseName, firstScanNum, lastScanNum)) {
            continue;
        }

        if (!(curBaseName.empty()) && curBaseName != baseName) {
            // a new run, produce normalization function, then clear m_rtLandmarks

            if (m_params.normalizeRTLinearRegression) {
                pair<double, double> coeff;
                if (doNormalizeRTLinearRegression(coeff, curBaseName)) {
                    m_normalizeRTCoeff[curBaseName] = coeff;
                    // cerr << "PUSH " << curBaseName << " : " << coeff.first << ", " << coeff.second << endl;
                }
            } else {
                map<double, double> pivots;
                if (doNormalizeRTLinearInterpolation(pivots, curBaseName)) {
                    m_normalizeRTPivots[curBaseName] = pivots;
                }
            }

            for (map<string, pair<double, vector<double> > >::iterator i = m_rtLandmarks->begin();
                 i != m_rtLandmarks->end(); i++) {
                i->second.second.clear();
            }

        }

        string rRTStr("");
        double rRT = 0.0;
        if (entry->getOneComment("RetentionTime", rRTStr)) {

            string::size_type commaPos = 0;
            rRT = atof(nextToken(rRTStr, commaPos, commaPos, ", \t\r\n").c_str());
            found->second.second.push_back(rRT);

        } else {
            // TODO: go to mzXML to find it,
        }

        curBaseName = baseName;

        entry->setOneComment("Remark", "Landmark");

    }

    // last one

    if (m_params.normalizeRTLinearRegression) {
        pair<double, double> coeff;
        if (doNormalizeRTLinearRegression(coeff, curBaseName)) {
            m_normalizeRTCoeff[curBaseName] = coeff;
        }
    } else {
        map<double, double> pivots;
        if (doNormalizeRTLinearInterpolation(pivots, curBaseName)) {
            m_normalizeRTPivots[curBaseName] = pivots;
        }
    }
}

bool SpectraSTPepXMLLibImporter::doNormalizeRTLinearRegression(pair<double, double> &coeff, string &baseName) {

    if (!m_rtLandmarks) return (false);

    double minRsq = 0.95;
    double minCoverage = 0.6;

//  cerr << "LINEAR REGRESSION" << endl;

    vector<double> vX;
    vector<double> vY;
    int count = 0;
    stringstream log1ss;
    
    for (map<string, pair<double, vector<double> > >::iterator i = m_rtLandmarks->begin();
         i != m_rtLandmarks->end(); i++) {

        double x = i->second.first;

        if (i->second.second.empty()) continue;

        double y = SpectraSTReplicates::getMedian(i->second.second);

          cerr << "LANDMARKS: " << i->first << "\t" << x << "\t" << y << endl;
          log1ss  << "LANDMARKS: " << i->first << "\t" << x << "\t" << y << endl;

        vX.push_back(x);
        vY.push_back(y);
        count++;

    }

    double rsq = linearRegress(vX, vY, coeff);

    
    log1ss << "RT normalization by linear regression. Found " << count << " landmarks in MS run \"" << baseName
           << "\". ";
    g_log->log("PEPXML IMPORT", log1ss.str());

    if (rsq < 0.000001) {

        g_log->error("PEPXML IMPORT", "Too few landmarks with distinct iRTs to perform RT normalization.");
        return (false);

    }

    //   cerr << "iRT = (rRT - " << coeff.second << ") / " << coeff.first << endl;

    while (rsq < minRsq && ((double) (vX.size()) > (double) (count) * minCoverage)) {

        vector<double> residuals;
        int pos = calcLinearRegressionResiduals(vX, vY, coeff, residuals);

        // remove outlier -- defined here as the one with largest absolute residual
        vX.erase(vX.begin() + pos);
        vY.erase(vY.begin() + pos);

        // redo regression
        rsq = linearRegress(vX, vY, coeff);
    }

    stringstream log2ss;
    log2ss << "Final fitted equation: iRT = (rRT ";

    log2ss.precision(4);
    if (coeff.second >= 0.0) {
        log2ss << "- " << coeff.second;
    } else {
        log2ss << "+ " << -coeff.second;
    }
    log2ss.precision(4);
    log2ss << ") / (" << coeff.first << "); ";
    log2ss.precision(4);
    log2ss << "R^2 = " << fixed << rsq << "; ";
    log2ss << count - vX.size() << " outliers removed.";

    g_log->log("PEPXML_IMPORT", log2ss.str());

    if (rsq >= minRsq) {

        return (true);

    } else {

        // not enough coverage

        g_log->error("PEPXML_IMPORT",
                     "R^2 still too low at required coverage. No RT normalization performed. Consider interpolation instead.");
        return (false);
    }

    return (false);

}

double SpectraSTPepXMLLibImporter::linearRegress(vector<double> &vX, vector<double> &vY, pair<double, double> &coeff) {

    double sumX = 0.0;
    double sumY = 0.0;
    double sumXX = 0.0;
    double sumXY = 0.0;
    double sumYY = 0.0;
    int count = 0;

    coeff.first = 0.0;
    coeff.second = 0.0;

    for (int pos = 0; pos < vX.size() && pos < vY.size(); pos++) {

        double x = vX[pos];
        double y = vY[pos];

        count++;
        sumX += x;
        sumY += y;
        sumXX += x * x;
        sumXY += x * y;
        sumYY += y * y;
    }

    if (count < 2) {

        return (0.0);

    }

    sumX /= (double) count;
    sumY /= (double) count;
    sumXX /= (double) count;
    sumXY /= (double) count;
    sumYY /= (double) count;

    double varX = sumXX - sumX * sumX;
    double varY = sumYY - sumY * sumY;
    double covar = sumXY - sumX * sumY;

    if (varX < 0.0001) {

        return (0.0);

    }

    coeff.first = covar / varX;
    coeff.second = sumY - coeff.first * sumX;
    double rsq = covar * covar / (varX * varY);

    return (rsq);
}

int SpectraSTPepXMLLibImporter::calcLinearRegressionResiduals(vector<double> &vX, vector<double> &vY,
                                                              pair<double, double> &coeff, vector<double> &residuals) {

    double maxAbsResidual = 0.0;
    int maxAbsResidualPos = -1;

    residuals.clear();

    for (int pos = 0; pos < vX.size() && pos < vY.size(); pos++) {

        double x = vX[pos];
        double y = vY[pos];

        double res = y - (coeff.first * x - coeff.second);
        residuals.push_back(res);

        if (fabs(res) > maxAbsResidual) {
            maxAbsResidual = fabs(res);
            maxAbsResidualPos = pos;
        }

    }

    return (maxAbsResidualPos);

}


bool SpectraSTPepXMLLibImporter::doNormalizeRTLinearInterpolation(map<double, double> &pivots, string &baseName) {

    if (!m_rtLandmarks) return (false);

    pivots.clear();

    // cerr << "LINEAR INTERPOLATION" << endl;

    for (map<string, pair<double, vector<double> > >::iterator i = m_rtLandmarks->begin();
         i != m_rtLandmarks->end(); i++) {
        double x = i->second.first;

        if (i->second.second.empty()) continue;

        double y = SpectraSTReplicates::getMedian(i->second.second);

        // cerr << "LANDMARKS: " << i->first << "\t" << x << "\t" << y << endl;

        pivots[y] = x;
    }

    stringstream logss;
    logss << "RT normalization by linear interpolation. Found " << pivots.size() << " landmarks in MS run \""
          << baseName << "\". ";


    bool success = false;
    if (pivots.size() < 2) {
        logss << "Too few to perform RT normalization.";
    } else {

        success = true;

        // check for inversion
        double lastIRT = -1000000000.0;
        for (map<double, double>::iterator j = pivots.begin(); j != pivots.end(); j++) {
            // rRT (first) is ascending, expect iRT (second) to be ascending too
            if (j->second < lastIRT) {
                // inversion detected
                success = false;
                logss << "RT inversion detected at iRT = " << j->second << ". No RT normalization performed.";
                break;
            }
            lastIRT = j->second;
        }
    }

    if (success) {
        g_log->log("PEPXML IMPORT", logss.str());
    } else {
        g_log->error("PEPXML IMPORT", logss.str());
    }

    return (success);

}


