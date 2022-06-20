#include "SpectraSTTsvLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "ProgressCount.hpp"
#include "XMLWalker.hpp"
#include <iostream>
#include <sstream>


/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
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
 * If multiple instances of the same spectra are present, the importer will end up importing the last instance.
 * 
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTTsvLibImporter::SpectraSTTsvLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                 SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params),
        m_mzXMLFiles(),
        m_queries(),
        m_numMzXMLOpen(0),
        m_datasetName(""),
        m_numSkipped(0),
        m_numPassedProbCutoff(0) {
//    cout <<__FILE__ << "\t"  << __FUNCTION__  << ": " << __LINE__ << endl;

}

// destructor
SpectraSTTsvLibImporter::~SpectraSTTsvLibImporter() {

    for (map<string, cRamp *>::iterator i = m_mzXMLFiles.begin(); i != m_mzXMLFiles.end(); i++) {
        if (i->second) {
            delete (i->second);
        }
    }

}

// import - imports all the pepXML files
void SpectraSTTsvLibImporter::import() {
//    cout <<__FILE__ << "\t"  << __FUNCTION__  << ": " << __LINE__ << endl;

    // read the pepXML files
    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        readFromFile(*i);
    }

    // by now we have all entries to be imported in memory
    m_lib->writePreamble(m_preamble);

    ProgressCount pc(!g_quiet, 1, (int) (m_queries.size()));
    pc.start("Importing spectra");

    // for efficiency, the entries don't have the spectra yet. Now we load the spectra.
    for (map<string, pair<vector<string>, SpectraSTLibEntry *> >::iterator q = m_queries.begin();
         q != m_queries.end(); q++) {

        pc.increment();

        string query = q->first;
        string path = (q->second.first.size() > 0 ? q->second.first[0] : "");
        string altPath = (q->second.first.size() > 1 ? q->second.first[1] : "");
//        cout <<  query << " " << path << " " << altPath << endl;
//        cout <<__FILE__ << "\t"  << __FUNCTION__  << ": " << __LINE__ << endl;
        SpectraSTLibEntry *entry = q->second.second;
        SpectraSTPeakList *peakList = entry->getPeakList();

        if (!(loadSpectrum(query, path, entry, altPath))) {
            // can't load spectrum
            m_numSkipped++;
            delete (entry);
            continue;
        }

        if (entry->getMods().find("iTRAQ") != string::npos) {
            peakList->removeITRAQReporterPeaks();
        }

        if (entry->getMods().find("TMT") != string::npos) {
            peakList->removeTMTReporterPeaks();
        }

        stringstream maxss;
        maxss.precision(2);
        maxss << peakList->getOrigMaxIntensity();
        entry->setOneComment("OrigMaxIntensity", maxss.str());

        if (insertOneEntry(entry, "TSV ")) {
            m_count++;
        } else {
            m_numSkipped++;
        }

        delete (entry);

    }

    stringstream countss;
    countss << "Total of " << m_count << " spectra imported, ";
    countss << m_numSkipped << " spectra skipped.";
    g_log->log("TSV IMPORT", countss.str());

    pc.done();

}

// readFromCurFile - reads one pepXML file
void SpectraSTTsvLibImporter::readFromFile(string &impFileName) {

    FileName fn;
    parseFileName(impFileName, fn);

    string fullImpFileName(impFileName);
    makeFullPath(fullImpFileName);

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

    string searchEngine("UNKNOWN");
    string database("UNKNOWN");
    string databaseType("AA");
    string instrumentType("UNKNOWN");

    ifstream fin;
    if (!myFileOpen(fin, impFileName)) {
        g_log->error("CREATE", "Cannot open TSV file \"" + impFileName + " for reading IDs. File Skipped.");
        return;
    }

    string line("");

    ProgressCount pc(!g_quiet, 500, 0);
    string message("Processing \"");
    message += impFileName + "\"";
    pc.start(message);

    while (nextLine(fin, line)) {
        string::size_type pos = 0;
        pos = line.find('#');
        if (pos != string::npos) {
            line = line.substr(0, pos);
        }
        string mzXMLFile = nextToken(line, 0, pos, " \t\r\n", " \t");
        if (mzXMLFile.empty()) {
            continue;
        }
        string scanNumStr = nextToken(line, pos, pos, " \t\r\n", " \t");
        int scanNum = atoi(scanNumStr.c_str());
        if (scanNum <= 0) {
            m_numSkipped++;
            continue;
        }
        string id = nextToken(line, pos, pos, " \t\r\n", " \t");

        string probStr = nextToken(line, pos, pos, " \t\r\n", " \t");
        double prob = atof(probStr.c_str());

        string scores = nextToken(line, pos, pos, " \t\r\n", " \t");

        string protein = nextToken(line, pos, pos, " \t\r\n", " \t");

        if (processSearchHit(mzXMLFile, scanNum, id, searchEngine, instrumentType, fn.path, prob, scores, protein)) {
            pc.increment();
        }

    }

    // describe what we are doing in the preamble
    stringstream filess;
    filess << "\"" << fullImpFileName << "\"";
    string desc = m_params.constructDescrStr(filess.str(), ".tsv");
    m_preamble.push_back(desc);
    g_log->log("CREATE", desc);

    pc.done();

}

bool SpectraSTTsvLibImporter::processSearchHit(string &mzXMLFile, int scanNum, string &peptide, string &searchEngine,
                                               string &instrumentType, string &path, double &prob, string &scores,
                                               string &protein) {

    if (prob < m_params.minimumProbabilityToInclude) {
        return (false);
    }
    m_numPassedProbCutoff++;

    FileName fn;
    parseFileName(mzXMLFile, fn);

    stringstream queryss;
    queryss << fn.name << '.' << scanNum << '.' << scanNum;
    string query = queryss.str();

    if (peptide[0] == '_') {
        // non peptide -- unknown precursor m/z, put in 0.0 for now
        SpectraSTLibEntry *entry = new SpectraSTLibEntry(peptide, 0.0, "", "Normal", NULL);
        entry->getPeakList()->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);
        entry->setOneComment("Spec", "Raw");

        if (instrumentType == "Ion Trap") {
            entry->setOneComment("Inst", "1/it,1");
        } else if (instrumentType == "TOF") {
            entry->setOneComment("Inst", "1/tof,1");
        } else if (instrumentType == "Q-TOF") {
            entry->setOneComment("Inst", "1/qtof,1");
        } else if (instrumentType == "Quadruple Ion Trap") {
            entry->setOneComment("Inst", "1/it,1");
        } else if (!instrumentType.empty()) {
            entry->setOneComment("Inst", "1/" + instrumentType + ",1");
        }

        string sample(m_datasetName);
        if (m_params.addMzXMLFileToDatasetName) {
            string::size_type dummy = 0;
            entry->setOneComment("Sample", "1/" + sample + "_" + fn.name + ",1,1");
        } else {
            entry->setOneComment("Sample", "1/" + sample + ",1,1");
        }

        entry->setOneComment("Nreps", "1/1");

        entry->setOneComment("RawSpectrum", query);

        stringstream probss;
        probss.precision(4);
        probss << fixed << prob;
        entry->setOneComment("Prob", probss.str());

        string altPath = fn.path;
        makeFullPath(altPath);

        // add to a structure in memory for now
        pair<vector<string>, SpectraSTLibEntry *> p;
        p.first.push_back(path);
        p.first.push_back(altPath);
        p.second = entry;
        m_queries[query] = p;

        return (true);

    }

    Peptide *pep = createPeptide(peptide, 0, "", query, "TSV");

    if (!pep) {
        m_numSkipped++;
        return (false);
    }

    if (m_params.setDeamidatedNXST) {
        setDeamidatedNXST(pep);
    }

    SpectraSTLibEntry *entry = new SpectraSTLibEntry(pep, "", "Normal", NULL, "");
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


    stringstream naass;
    naass << pep->NAA();
    entry->setOneComment("NAA", naass.str());

    entry->setOneComment("Mods", pep->mspMods());

    if (!(protein.empty())) {
        entry->setOneComment("Protein", protein);
    }

    stringstream parentss;
    parentss << entry->getPrecursorMz();
    entry->setOneComment("Parent", parentss.str());

    if (instrumentType == "Ion Trap") {
        entry->setOneComment("Inst", "1/it,1");
    } else if (instrumentType == "TOF") {
        entry->setOneComment("Inst", "1/tof,1");
    } else if (instrumentType == "Q-TOF") {
        entry->setOneComment("Inst", "1/qtof,1");
    } else if (instrumentType == "Quadruple Ion Trap") {
        entry->setOneComment("Inst", "1/it,1");
    } else if (!instrumentType.empty()) {
        entry->setOneComment("Inst", "1/" + instrumentType + ",1");
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

    // parse scores
    map<string, double> scs;
    string sc("");
    string::size_type commaPos = 0;
    while (!(sc = nextToken(scores, commaPos, commaPos, ",\t\r\n", ",")).empty()) {
        string::size_type equalPos = sc.find('=');
        if (equalPos != string::npos) {
            string scName = sc.substr(0, equalPos);
            double scValue = atof(sc.substr(equalPos + 1).c_str());
            scs[scName] = scValue;
            searchss.precision(4);
            searchss << "," << scName << "=" << fixed << scValue << "/0";
        }
    }

    if (se == 'S') {
        if (scs.find("xc") != scs.end()) {
            double xcorr = scs["xc"];
            stringstream xcorrss;
            xcorrss.precision(4);
            xcorrss << fixed << xcorr;
            entry->setOneComment("XCorr", xcorrss.str());
        }

        if (scs.find("dc") != scs.end()) {
            double deltacn = scs["dc"];
            stringstream deltacnss;
            deltacnss.precision(4);
            deltacnss << fixed << deltacn;
            entry->setOneComment("DeltaCn", deltacnss.str());
        }
    }

    entry->setOneComment("Se", searchss.str());

    string sample(m_datasetName);


    if (m_params.addMzXMLFileToDatasetName) {
        string::size_type dummy = 0;
        entry->setOneComment("Sample", "1/" + sample + "_" + fn.name + ",1,1");
    } else {
        entry->setOneComment("Sample", "1/" + sample + ",1,1");
    }

    entry->setOneComment("Nreps", "1/1");

    entry->setOneComment("RawSpectrum", query);

    stringstream probss;
    probss.precision(4);
    probss << fixed << prob;
    entry->setOneComment("Prob", probss.str());

    string altPath = fn.path;
    makeFullPath(altPath);

    // add to a structure in memory for now
    pair<vector<string>, SpectraSTLibEntry *> p;
    p.first.push_back(path);
    p.first.push_back(altPath);
    p.second = entry;
    m_queries[query] = p;

    return (true);

}


// loadSpectrum - go to the mzXML file to load the spectrum
bool SpectraSTTsvLibImporter::loadSpectrum(string query, string path, SpectraSTLibEntry *entry, string altPath) {
//    cout <<__FILE__ << "\t"  << __FUNCTION__  << ": " << __LINE__ << endl;
    SpectraSTPeakList *peakList = entry->getPeakList();

    string::size_type lastDotPos = query.rfind('.');   // .<lastScanNum>  (<charge> has been taken off already)
    if (lastDotPos == string::npos) {
        g_log->error("TSV IMPORT", "Illegal query name \"" + query + "\". Scan not imported.");
        return (false);
    }
    string::size_type secondLastDotPos = query.rfind('.', lastDotPos - 1);  // .<firstScanNum>.<lastScanNum>
    if (secondLastDotPos == string::npos) {
        g_log->error("TSV IMPORT", "Illegal query name \"" + query + "\". Scan not imported.");
        return (false);
    }

    int scanNum = atoi((query.substr(secondLastDotPos + 1, lastDotPos - secondLastDotPos - 1)).c_str());

    string baseName = query.substr(0, secondLastDotPos);

    // even though the scans may not be in fullFileName, still use this as the key in the map
    string fullFileName = path + baseName + ".mzXML";

    // try to see if the file is already opened
    map<string, cRamp *>::iterator found = m_mzXMLFiles.find(fullFileName);
    cRamp *cramp = NULL;
    if (found == m_mzXMLFiles.end()) {
        // can't find this mzXML file
        if (m_numMzXMLOpen >= MAX_NUM_OPEN_FILES) {
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

        string triedFileNames(fullFileName);
        string tryFileName(fullFileName);
        cramp = new cRamp(tryFileName.c_str());


        if (!cramp->OK()) {
            // try mzML instead
            delete (cramp);
            tryFileName = path + baseName + ".mzML";
            triedFileNames += "|" + tryFileName;
            cramp = new cRamp(tryFileName.c_str());
        }

        if (!(altPath.empty())) {
            if (!cramp->OK()) {
                delete (cramp);
                tryFileName = altPath + baseName + ".mzXML";
                triedFileNames += "|" + tryFileName;
                cramp = new cRamp(tryFileName.c_str());
            }

            if (!cramp->OK()) {
                delete (cramp);
                tryFileName = altPath + baseName + ".mzML";
                triedFileNames += "|" + tryFileName;
                cramp = new cRamp(tryFileName.c_str());
            }
        }

        if (!cramp->OK()) {
            g_log->error("TSV IMPORT",
                         "Cannot open file \"" + triedFileNames + "\". No scan from this file will be imported.");
            delete (cramp);
            m_mzXMLFiles[fullFileName] = NULL; // note this so that future attempts to read scan from this file will die silently
            return (false);
        } else {
            m_numMzXMLOpen++;
            m_mzXMLFiles[fullFileName] = cramp;
        }

    } else {
        if (found->second) {
            cramp = found->second;
        } else {
            // tried to open this file previously and failed
            return (false);
        }
    }

    rampScanInfo *scanInfo = cramp->getScanHeaderInfo(scanNum);
    if (!scanInfo || scanInfo->m_data.acquisitionNum != scanNum || scanInfo->m_data.msLevel == 1) {
        // bad. not found!
        stringstream errss;
        errss << "Cannot find MS2+ scan #" << scanNum << " in file \"" << fullFileName << "\". Scan not imported.";
        g_log->error("PEPXML IMPORT", errss.str());

        if (scanInfo) delete (scanInfo);
        return (false);
    }

//  cout << "get peaks" << endl;
    // now we can read the peaks
    rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);
    if (!peaks) {
        stringstream errss;
        errss << "Cannot read peaks for scan #" << scanNum << " in file \"" << fullFileName << "\". Scan not imported.";
        g_log->error("PEPXML IMPORT", errss.str());
        return (false);
    }

    int peakCount = peaks->getPeakCount();
    double precursorMz = scanInfo->m_data.precursorMZ;
    int precursorCharge = scanInfo->m_data.precursorCharge;
    if (precursorCharge < 1) precursorCharge = 0;
    string fragType(scanInfo->m_data.activationMethod);
    double rt = scanInfo->m_data.retentionTime;
    entry->setOneComment("RetentionTime",to_string(rt)+","+to_string(rt)+","+to_string(rt));
    // set mass diff also
    entry->setOneComment("MassDiff","0.01");

    if (!(fragType.empty())) {
        entry->setFragType(fragType);
    }

    if (entry->getPrecursorMz() < 0.0001) entry->setPrecursor(precursorMz, precursorCharge);

//  cout << "inserting peaks " << peakCount << endl;
    // create the peak list and read the peaks one-by-one
    for (int j = 0; j < peakCount; j++) {
        if (peaks->getPeak(j)->intensity > 0.1) {
            peakList->insert(peaks->getPeak(j)->mz, (float) (peaks->getPeak(j)->intensity), "", "");
        }
    }

    delete peaks;

    return (true);

}


