#include "SpectraSTMzXMLLibImporter.hpp"
#include "SpectraSTReplicates.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>


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

/* Class: SpectraSTMzXMLLibImporter
 * 
 * Implements a library importer for the .mzXML file format.
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor - will open the first file
SpectraSTMzXMLLibImporter::SpectraSTMzXMLLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                     SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params),
        m_numScansInFile(0),
        m_numMissingInFile(0),
        m_numMS1InFile(0),
        m_numFailedFilterInFile(0),
        m_numImportedInFile(0),
        m_numConsensusInFile(0),
        m_numBadConsensusInFile(0),
        m_datasetName("") {


}

// destructor 
SpectraSTMzXMLLibImporter::~SpectraSTMzXMLLibImporter() {

}

// import - prints the preamble, then loops over all files and import them one by one
void SpectraSTMzXMLLibImporter::import() {

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        string fullName(*i);
        makeFullPath(fullName);


        if (m_params.datasetName.empty()) {
            // not specified... use the full path and the directory name as the dataset name
            // to be safe, replace special characters with _
            stringstream dnss;
            string forbidden("~`!@#$%^&*()+\\|:<>,.?/{}'\"");
            string path("");
            getPath(fullName, path);
            for (string::size_type pos = 0; pos < path.length() - 1; pos++) {  // -1 to take off the last /
                if (forbidden.find(path[pos]) != string::npos) {
                    dnss << '_';
                } else {
                    dnss << path[pos];
                }
            }
            m_datasetName = dnss.str();
        } else {
            m_datasetName = m_params.datasetName;
        }

        string quoted("\"" + fullName + "\"");
        string desc = m_params.constructDescrStr(quoted, ".mzXML");
        m_preamble.push_back(desc);
    }

    m_lib->writePreamble(m_preamble);

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        readFromFile(*i);
    }

}

// readFromFile - reads one .mzXML file
void SpectraSTMzXMLLibImporter::readFromFile(string &impFileName) {

    g_log->log("MZXML IMPORT", "Importing .mzXML file \"" + impFileName + "\".");


    // open the file using cRamp
    cRamp *cramp = new cRamp(impFileName.c_str());

    if (!cramp->OK()) {
        g_log->error("MZXML IMPORT", "Cannot open file \"" + impFileName + "\". File skipped.");
        delete (cramp);
        return;
    }

    // Read the run info to extract the number of scans
    rampRunInfo *runInfo = cramp->getRunInfo();

    if (!runInfo) {
        // probably an empty file...
        g_log->error("MZXML IMPORT", "Cannot open file \"" + impFileName + "\". File skipped.");
        delete (cramp);
        return;
    }

    //  rampInstrumentInfo* instr = cramp->getInstrumentInfo();

    int numScans = cramp->getLastScan();
    delete (runInfo);


    if (g_verbose) {
        cout << "\nImporting spectra from \"" << impFileName << "\"" << endl;
    }

    // start the progress count
    ProgressCount pc(!g_quiet && !g_verbose, 1, numScans);
    pc.start("\nImporting spectra from " + impFileName);

    // parse out the file name to determine the query prefix. Note that the query string
    // has the form <mzXML file name>.<scan num>.<scan num>.0
    FileName fn;
    parseFileName(impFileName, fn);

    m_numScansInFile = numScans;
    m_numMissingInFile = 0;
    m_numMS1InFile = 0;
    m_numFailedFilterInFile = 0;
    m_numImportedInFile = 0;
    m_numConsensusInFile = 0;
    m_numBadConsensusInFile = 0;

    for (int k = 1; k <= numScans; k++) {

        pc.increment();

        // get the scan header (no peak list) first to check whether it's MS2.
        // it'd be a waste of time if we read all scans, including MS1
        rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);

        // check to make sure the scan is good, and is not MS1
        if (!scanInfo || (scanInfo->m_data.acquisitionNum != k)) {
            m_numMissingInFile++;

            if (scanInfo) delete (scanInfo);
            continue;
        }

        if (scanInfo->m_data.msLevel == 1) {
            m_numMS1InFile++;
            delete (scanInfo);
            continue;
        }

        // now we can import
        importOne(cramp, scanInfo, fn.name);
        // done, can delete scanInfo
        delete scanInfo;

    }

    // flush out last remaining clusters
    for (map<double, vector<SpectraSTLibEntry *> *>::iterator i = m_clusters.begin(); i != m_clusters.end(); i++) {
        formConsensusEntry(i->second);
        delete (i->second);
    }

    pc.done();

    m_clusters.clear();

    // log the import of this file
    stringstream importLogss;
    importLogss << "Imported \"" << impFileName + "\" ";
    importLogss << "(Max " << m_numScansInFile << " scans; " << m_numImportedInFile << " imported, ";
    importLogss << m_numFailedFilterInFile << " failed filter; " << m_numMissingInFile << " missing; " << m_numMS1InFile
                << " MS1)";
    importLogss << " (" << m_numConsensusInFile << " clusters, " << m_numBadConsensusInFile << " bad)";
    g_log->log("MZXML IMPORT", importLogss.str());

    // we can delete the cRamp object now that we're done with this file.
    delete (cramp);
}

// import - reads one MS2 spectrum  
void SpectraSTMzXMLLibImporter::importOne(cRamp *cramp, rampScanInfo *scanInfo, string &prefix) {

    unsigned int scanNum = scanInfo->m_data.acquisitionNum;

    stringstream namess;
    namess << '_' << prefix << '_';
    namess.width(5);
    namess.fill('0');
    namess << right << scanNum;

    // Go back to the mzXML file and get the peaks using Ramp
    rampPeakList *peaks = cramp->getPeakList(scanNum);
    if (!peaks) {
        m_numFailedFilterInFile++;
        return;
    }

    int peakCount = peaks->getPeakCount();
    double precursorMz = scanInfo->m_data.precursorMZ;
    int precursorCharge = scanInfo->m_data.precursorCharge;
    if (precursorCharge < 1) precursorCharge = 0;
    double retentionTime = scanInfo->getRetentionTimeSeconds();
    double precursorIntensity = scanInfo->m_data.precursorIntensity;
    double totIonCurrent = scanInfo->m_data.totIonCurrent;

    string fragType("");
    if (scanInfo->m_data.activationMethod) fragType = scanInfo->m_data.activationMethod;

    if (!(m_params.setFragmentation.empty())) {
        // allow override of frag type if the user explicitly specifies it
        fragType = m_params.setFragmentation;
    }
    // create the peak list and read the peaks one-by-one
    SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, precursorCharge, peakCount, false, fragType);
    peakList->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);

    for (int j = 0; j < peakCount; j++) {
        double mz = peaks->getPeak(j)->mz;
        float inten = (float) (peaks->getPeak(j)->intensity);
        peakList->insert(mz, inten, "", "");
    }

    delete peaks;

    if (m_params.centroidPeaks) {
        // check instrument? doesn't seem to matter much, so just assume it's TOF for now
        peakList->centroid("TOF");
    }

    // normalize
    if (!(m_params.keepRawIntensities)) {
        peakList->normalizeTo(10000.0, m_params.rawSpectraMaxDynamicRange);
    }

    if (!(peakList->passFilterUnidentified(m_params))) {
        m_numFailedFilterInFile++;
        delete peakList;
        return;
    }

    stringstream commentss;

    commentss << " Nreps=1/1";
    commentss.precision(2);
    commentss << " PrecursorIntensity=" << precursorIntensity;
    commentss << " Prob=1.000";
    commentss << " RawSpectrum=" << prefix << '.';
    commentss.width(5);
    commentss.fill('0');
    commentss << right << scanNum;
    commentss << '.';
    commentss.width(5);
    commentss.fill('0');
    commentss << right << scanNum;
    commentss.precision(1);
    commentss << " RetentionTime=" << fixed << retentionTime << ',' << retentionTime << ',' << retentionTime;

    string sample(m_datasetName);
    if (m_params.addMzXMLFileToDatasetName) {
        commentss << " Sample=1/" << sample << "_" << prefix << ",1,1";
    } else {
        commentss << " Sample=1/" << sample << ",1,1";
    }

    double sn = peakList->calcSignalToNoise();
    commentss.precision(1);
    commentss << " SN=" << fixed << sn;

    commentss.precision(2);
    commentss << " TotalIonCurrent=" << totIonCurrent;

    commentss.precision(2);
    commentss << " OrigMaxIntensity=" << peakList->getOrigMaxIntensity();


    double xrea = peakList->calcXrea(true);
    commentss.precision(3);
    commentss << " Xrea=" << fixed << xrea;

    namess << '/' << precursorCharge;

    SpectraSTLibEntry *entry = new SpectraSTLibEntry(namess.str(), precursorMz, commentss.str(), "Normal", peakList,
                                                     fragType);

    m_numImportedInFile++;


    // clustering
    //
    // The following code tries to detect replicates, merge them by consensus creation, then insert the
    // merged spectrum into the library instead of the raw spectra. It is based on a "lazy" or "smart" (depending
    // on your perspective!) clustering approach, whereby only scans acquired at the same precursor m/z AND
    // without a dissimilar scan acquired between them (on the time axis) are potential candidates to be merged.
    // This implies that clusters have to be kept only as long as there is no dissimilar scan has been acquired
    // at the same precursor m/z.
    //

    if (m_params.unidentifiedClusterIndividualRun) {

        map<double, vector<SpectraSTLibEntry *> *>::iterator lower = m_clusters.upper_bound(precursorMz - 1.0);
        map<double, vector<SpectraSTLibEntry *> *>::iterator upper = m_clusters.lower_bound(precursorMz + 1.0);

        if (lower == upper) {
            // form own cluster
            //    cerr << "Form cluster: " << namess.str() << " (" << precursorMz << ")" << endl;
            vector<SpectraSTLibEntry *> *cluster = new vector<SpectraSTLibEntry *>;
            cluster->push_back(entry);
            m_clusters[precursorMz] = cluster;

        } else {

            bool clustered = false;
            for (map<double, vector<SpectraSTLibEntry *> *>::iterator i = lower; i != upper; i++) {
                if (!(i->second) || i->second->empty()) continue; // should not happen
                vector<SpectraSTLibEntry *> *cluster = i->second;
                if (!clustered && (*cluster)[cluster->size() - 1]->getPeakList()->compare(peakList) >
                                  m_params.unidentifiedClusterMinimumDot) {
                    // cluster!
                    // cerr << "Join cluster: " << namess.str() << " (" << precursorMz << ") => " << i->first << endl;
                    cluster->push_back(entry);
                    clustered = true;

                } else {
                    // a new scan at same m/z but cannot cluster, this cluster should not carry on any more
                    // cerr << "Close cluster: " << cluster->size()<< " entries (" << i->first << ")" << endl;
                    formConsensusEntry(cluster);
                    delete (cluster);
                    m_clusters.erase(i);
                }
            }
            if (!clustered) {
                // form own cluster
                // cerr << "Form cluster: " << namess.str() << " (" << precursorMz << ")" << endl;
                vector<SpectraSTLibEntry *> *cluster = new vector<SpectraSTLibEntry *>;
                cluster->push_back(entry);
                m_clusters[precursorMz] = cluster;
            }
        }

    } else {
        // no clustering, just put in library

        if (passAllFilters(entry)) {
            m_lib->insertEntry(entry);
            m_count++;
        }

        delete (entry);
    }

}

// formConsensusEntry - uses the consensus creation algorithm to merge all spectra in a cluster, then inserts
// the merged spectrum into the library. It will also delete all SpectraSTLibEntry's passed in.
void SpectraSTMzXMLLibImporter::formConsensusEntry(vector<SpectraSTLibEntry *> *cluster) {

    SpectraSTReplicates reps(*cluster, m_params);
    reps.setRecordRawSpectra(true);
    SpectraSTLibEntry *consensus = reps.makeConsensusSpectrum();
    if (consensus && passAllFilters(consensus)) {
        if (!(consensus->getPeakList()->passFilterUnidentified(m_params))) {
            // bad one after consensus, likely noise?
            m_numBadConsensusInFile++;
        } else {
            m_lib->insertEntry(consensus);
            m_count++;
        }

        m_numConsensusInFile++;
    }

    for (vector<SpectraSTLibEntry *>::iterator en = cluster->begin(); en != cluster->end(); en++) {
        delete (*en);
    }


}
