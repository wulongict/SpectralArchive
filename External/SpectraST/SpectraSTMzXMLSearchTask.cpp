#include "SpectraSTMzXMLSearchTask.hpp"
#include "SpectraSTPeakList.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "ProgressCount.hpp"

#include <sstream>
#include <algorithm>
#include <cmath>

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

/* Class: SpectraSTMzXMLSearchTask
 * 
 * Subclass of SpectraSTSearchTask that handles.mzXML file types. Note that because of
 * the availability of a random access index at the end of the mzXML, conveniently one can
 * sort the spectra by ascending precursor m/z and search them in that order to take advantage
 * of caching.
 * 
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;


// constructor - opens all mzXML files with cRAMP. For efficiency, the queries from multiple files are merged and
// sorted by precursor m/z before searching. 
SpectraSTMzXMLSearchTask::SpectraSTMzXMLSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                                   SpectraSTLib *lib) :
        SpectraSTSearchTask(searchFileNames, params, lib),
        m_files(searchFileNames.size()),
        m_scans(),
        m_isMzData(false),
        m_numScansInFile(0),
        m_numNotSelectedInFile(0),
        m_numMissingInFile(0),
        m_numMS1InFile(0),
        m_numFailedFilterInFile(0),
        m_numSearchedInFile(0),
        m_numLikelyGoodInFile(0) {

    char *rampExt = rampValidFileType(m_searchFileNames[0].c_str());
    if (strstri(rampExt, ".mzdata")) { // accommodates .mzdata.gz
        m_isMzData = true;
    }


}


// destructor - deletes the cRamp objects -- essentially closing the mzXML files too
SpectraSTMzXMLSearchTask::~SpectraSTMzXMLSearchTask() {

}


// search - runs the searches
void SpectraSTMzXMLSearchTask::search() {

    if (!m_params.indexCacheAll) {

        // Not caching all entries. In this case, the queries have to be sorted by precursor m/z first, such that
        // the cached window slides from low to high precursor m/z only once. (Otherwise, the cached entries will need to be swapped
        // in and out repeatedly, defeating the purpose of caching.)

        // There is a catch however. To be able to search out of order, one has to keep many mzXML files open, and most systems
        // have a max file opened limit. In such case, we will need to divide the mzXML files into smaller batches. The library will
        // will need to be read numBatches times, but the tradeoff is we won't need to keep all entries cached in memory by selecting
        // the indexCacheAll option.

        // Divide the mzXML files into equal batches of at most MAX_NUM_OPEN_FILES files.
        unsigned int numBatches = ((unsigned int) m_searchFileNames.size() - 1) / MAX_NUM_OPEN_FILES + 1;
        unsigned int batchStart = 0;
        for (unsigned int b = 0; b < numBatches; b++) {
            m_batchBoundaries.push_back(batchStart);
            batchStart += (unsigned int) m_searchFileNames.size() / numBatches;
        }
        m_batchBoundaries.push_back((unsigned int) m_searchFileNames.size());

        // For each batch, sort all spectra by precursor m/z, open the files, and set the search in motion
        for (unsigned int batch = 0; batch < (unsigned int) m_batchBoundaries.size() - 1; batch++) {

            // this will do the sorting. the vector m_scans will be populated with the sorted scans.
            prepareSortedSearch((unsigned int) batch);

            // open the output files and print the headers (e.g. the xml definitions, ms run info, etc)
            for (unsigned int n = m_batchBoundaries[batch]; n < m_batchBoundaries[batch + 1]; n++) {
                m_outputs[n]->openFile();
                m_outputs[n]->printHeader();
            }

            // tracking search progress
            ProgressCount pc(!g_quiet && !g_verbose, 1, (int) (m_scans.size()));
            string msg("Searching");
            pc.start(msg);

            // create searches from the m_scans one-by-one, and search them
            for (vector<pair<unsigned int, rampScanInfo *> >::iterator i = m_scans.begin(); i != m_scans.end(); i++) {

                searchOne((*i).first, (*i).second);
                pc.increment();

                // done. we can delete the rampScanInfo object now.
                delete (*i).second;
            }
            pc.done();

            // log the search of the batch
            stringstream searchLogss;
            searchLogss << "Searched sorted scans ";
            searchLogss << "(Max " << m_numScansInFile << " scans; " << m_numSearchedInFile << " searched, ";
            searchLogss << m_numLikelyGoodInFile << " likely good; ";
            if (m_numNotSelectedInFile > 0) {
                searchLogss << m_numNotSelectedInFile << " not selected; ";
            }
            searchLogss << m_numFailedFilterInFile << " failed filter; " << m_numMissingInFile << " missing; "
                        << m_numMS1InFile << " MS1)";
            g_log->log("MZXML SEARCH", searchLogss.str());

            // done with this batch. close the files so that we can open more files in the next batch
            for (unsigned int n = m_batchBoundaries[batch]; n < m_batchBoundaries[batch + 1]; n++) {
                delete (m_files[n].second); // the cramp objects, which will close the mzXML files
                m_files[n].second = NULL;
                m_outputs[n]->printFooter(); // the output files
                m_outputs[n]->closeFile();
            }
        }

    } else {

        // This is the case where we're caching everything anyway. In this case, it is not necessary to sort
        // by precursor m/z before searching. We simply open the files one by one and search the queries
        // in the order they are read.

        // Fingerprinting
        m_searchFingerprintIndexedByLibID.assign(m_lib->getMzLibIndexPtr()->getEntryCount() + 1, 0);
        m_searchFingerprintIndexedBySearchID.assign(99999999, 0);
        // END Fingerprinting

        for (unsigned int n = 0; n < (unsigned int) m_searchFileNames.size(); n++) {


            // open the file using cRamp
            cRamp *cramp = new cRamp(m_searchFileNames[n].c_str());

            if (!cramp->OK()) {
                g_log->error("MZXML SEARCH", "Cannot open file \"" + m_searchFileNames[n] + "\". File skipped.");
                delete (cramp);
                continue;
            }

            // Read the run info to extract the number of scans
            rampRunInfo *runInfo = cramp->getRunInfo();

            if (!runInfo) {
                // probably an empty file...
                g_log->error("MZXML SEARCH", "Cannot open file \"" + m_searchFileNames[n] + "\". File skipped.");
                delete (cramp);
                continue;
            }

            rampInstrumentInfo *instr = cramp->getInstrumentInfo();
            if (instr) {
                m_outputs[n]->setInstrInfo(instr);
//        delete (instr);
            }

            // open the output file
            m_outputs[n]->openFile();
            m_outputs[n]->printHeader();

            int numScans = cramp->getLastScan();

            delete (runInfo);

            // parse out the file name to determine the query prefix. Note that the query string
            // has the form <mzXML file name>.<scan num>.<scan num>.0
            FileName fn;
            parseFileName(m_searchFileNames[n], fn);

            // m_files is a vector of (FileName, cRamp*)
            m_files[n].first = fn;
            m_files[n].second = cramp;

            ProgressCount pc(!g_quiet && !g_verbose, 1, numScans);
            stringstream msg;
            msg << "Searching \"" << m_searchFileNames[n] << "\" " << "(" << n + 1 << " of " << m_searchFileNames.size()
                << ")";
            pc.start(msg.str());

            m_numScansInFile = numScans;
            m_numNotSelectedInFile = 0;
            m_numMissingInFile = 0;
            m_numMS1InFile = 0;
            m_numFailedFilterInFile = 0;
            m_numSearchedInFile = 0;
            m_numLikelyGoodInFile = 0;


            for (int k = 1; k <= numScans; k++) {

                pc.increment();

                // Filter out all scans not in selected list (in this case,
                // the selected list contains a list of scan numbers as strings
                if (!m_searchAll && !isInSelectedList(SpectraSTQuery::constructQueryTPPStyle(fn.name, k, k, 0))) {
                    m_numNotSelectedInFile++;
                    continue;
                }
                // get the scan header (no peak list) first to check whether it's MS2.
                // it'd be a waste of time if we read all scans, including MS1
                rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);

                // check to make sure the scan is good, and is not MS1
                if (!scanInfo || (!m_isMzData && scanInfo->m_data.acquisitionNum != k)) {
                    m_numMissingInFile++;

                    if (scanInfo) delete (scanInfo);
                    continue;
                }

                if (scanInfo->m_data.msLevel == 1) {
                    m_numMS1InFile++;
                    delete (scanInfo);
                    continue;
                }

                // now we can search
                searchOne(n, scanInfo);
                // done, can delete scanInfo
                delete scanInfo;

            }
            pc.done();



            // log the search of this file
            stringstream searchLogss;
            searchLogss << "Searched \"" << m_searchFileNames[n] + "\" ";
            searchLogss << "(Max " << m_numScansInFile << " scans; " << m_numSearchedInFile << " searched, ";
            searchLogss << m_numLikelyGoodInFile << " likely good; ";
            if (m_numNotSelectedInFile > 0) {
                searchLogss << m_numNotSelectedInFile << " not selected; ";
            }
            searchLogss << m_numFailedFilterInFile << " failed filter; " << m_numMissingInFile << " missing; "
                        << m_numMS1InFile << " MS1)";
            g_log->log("MZXML SEARCH", searchLogss.str());

            // we can delete the cRamp object now that we're done with this file.
            // this is in contrast to the case where we're opening all the files at once for
            // sorting -- in that case the cRamp objects will be deleted at the end of all
            // searches
            delete (m_files[n].second);
            m_files[n].second = NULL;

            m_outputs[n]->printFooter();
            m_outputs[n]->closeFile(); // just so we won't hit the File Open limit if there are too many files
        }

        // Fingerprinting
        if (!(m_params.printFingerprintingSummary.empty())) {

            string fingerprintFileName(m_params.printFingerprintingSummary);
            ofstream fingerprintFout;

            if (!myFileOpen(fingerprintFout, fingerprintFileName)) {
                g_log->error("CREATE", "Cannot open file \"" + fingerprintFileName +
                                       "\" for writing fingerprinting summary. Writing skipped.");
                return;
            }

            SpectraSTLibEntry *tempEntry = NULL;
            map<string, pair<unsigned int, unsigned int> > samples;

            m_lib->getMzLibIndexPtr()->reset();

            while (tempEntry = m_lib->getMzLibIndexPtr()->nextEntry()) {
                tempEntry->getSampleInfo(samples, "USED_ONLY");
                delete (tempEntry);
            }

            // record all the sample source information into sampleSource with its index
            // (e.g._media_data_project_clustering_spc_zoo_blood_African_Lion  0
            // _media_data_project_clustering_spc_zoo_blood_Amur_Tiger  1
            // _media_data_project_clustering_spc_zoo_blood_Barred_Owl  2) They are sorted in an alphabetical order.

            map<string, int> sampleSource;
            map<string, pair<unsigned int, unsigned int> >::iterator iter;
            int tempCount = 0;

            for (iter = samples.begin(); iter != samples.end(); iter++) {
                sampleSource.insert(pair<string, int>(iter->first, tempCount));
                tempCount++;
            }

            m_bootstrapSupport.assign(sampleSource.size(), 0.0);

            fingerprintFout << "Columns: ";
            for (map<string, int>::iterator test = sampleSource.begin(); test != sampleSource.end(); test++) {
                fingerprintFout << "[" << test->second << "]" << test->first << "\t";
            }

            fingerprintFout << endl;

//      calcFingerprint(true, fingerprintFout);

//      fingerprintFout << "***bootstrapping starts***" << endl;

            vector<float> originalsearchFingerprintIndexedByLibID = m_searchFingerprintIndexedByLibID;

            for (int loopTime = 0; loopTime != 4999; loopTime++) {

                doBootstrapping();
                calcFingerprint(false, fingerprintFout);

            }

//      fingerprintFout << "***bootstrapping ends***" << endl;

            fingerprintFout << "bootstrap_support:  ";


            for (int i = 0; i < sampleSource.size(); i++) {
                fingerprintFout.precision(3);
                fingerprintFout << fixed << m_bootstrapSupport[i] / 5000.0 << " ";
            }

            fingerprintFout << endl;

            m_searchFingerprintIndexedByLibID = originalsearchFingerprintIndexedByLibID;

            calcFingerprint(true, fingerprintFout);

        }

        // END Fingerprinting

    }

    m_searchTaskStats.logStats();


}


// Fingerprinting

void SpectraSTMzXMLSearchTask::calcFingerprint(bool printFingerprint, ofstream &fingerprintFout) {

    vector<vector<float> > libFP = m_lib->getFingerprint();

    vector<float> dot, dot_sqrt, dot_binary, dot_unique, dot_unique_sqrt, spectralCount, spectralCount_unique, sq, sq_sqrt, sq_unique, sq_unique_sqrt, libSize, libSum, sqItself_unique, sqItself_unique_sqrt;
    dot.assign(libFP.size(), 0);
    dot_sqrt.assign(libFP.size(), 0);
    dot_binary.assign(libFP.size(), 0);
    dot_unique.assign(libFP.size(), 0);
    dot_unique_sqrt.assign(libFP.size(), 0);
    sq.assign(libFP.size(), 0);
    sq_sqrt.assign(libFP.size(), 0.001);
    sq_unique.assign(libFP.size(), 0.001);
    sq_unique_sqrt.assign(libFP.size(), 0001);
    spectralCount.assign(libFP.size(), 0);
    spectralCount_unique.assign(libFP.size(), 0);
    libSize.assign(libFP.size(), 0);
    libSum.assign(libFP.size(), 0);
    sqItself_unique.assign(libFP.size(), 0.001);
    sqItself_unique_sqrt.assign(libFP.size(), 0.001);
    float sqItself = 0;
    float sqItself_sqrt = 0;
    float sqItself_binary = 0;

    float total = 0;

    for (int i = 0; i != libFP[0].size(); i++) {
        for (int j = 0; j != libFP.size(); j++) {
            dot[j] += libFP[j][i] * m_searchFingerprintIndexedByLibID[i];
            dot_sqrt[j] += sqrt(libFP[j][i] * m_searchFingerprintIndexedByLibID[i]);
            if (libFP[j][i] > 0 && m_searchFingerprintIndexedByLibID[i] > 0) dot_binary[j]++;
            sq[j] += libFP[j][i] * libFP[j][i];
            sq_sqrt[j] += sqrt(libFP[j][i] * libFP[j][i]);
            if (libFP[j][i] > 0) spectralCount[j] += m_searchFingerprintIndexedByLibID[i];

            double numSpecies = 0;

            for (int k = 0; k != libFP.size(); k++) {
                if (libFP[k][i] != 0) {
                    numSpecies++;
                }
            }

            if (numSpecies == 1) {

                dot_unique[j] += m_searchFingerprintIndexedByLibID[i] * libFP[j][i];
                dot_unique_sqrt[j] += sqrt(m_searchFingerprintIndexedByLibID[i] * libFP[j][i]);
                sq_unique[j] += libFP[j][i] * libFP[j][i];
                sq_unique_sqrt[j] += sqrt(libFP[j][i] * libFP[j][i]);
                sqItself_unique[j] += m_searchFingerprintIndexedByLibID[i] * m_searchFingerprintIndexedByLibID[i];
                sqItself_unique_sqrt[j] += sqrt(
                        m_searchFingerprintIndexedByLibID[i] * m_searchFingerprintIndexedByLibID[i]);

                if (libFP[j][i] > 0) {

                    spectralCount_unique[j] += m_searchFingerprintIndexedByLibID[i];
                    libSum[j] += libFP[j][i];

                }

            }

            if (libFP[j][i] > 0) {
                libSize[j]++;
            }

        }


        sqItself += m_searchFingerprintIndexedByLibID[i] * m_searchFingerprintIndexedByLibID[i];
        sqItself_sqrt += sqrt(m_searchFingerprintIndexedByLibID[i] * m_searchFingerprintIndexedByLibID[i]);
        if (m_searchFingerprintIndexedByLibID[i] > 0) sqItself_binary++;

        total += m_searchFingerprintIndexedByLibID[i];
    }

    float sum = 0;
    for (int j = 0; j != libFP.size(); j++) {
        sum += dot[j] / libSize[j];
    }

    int maxIndex = 0;
    float max = 0;


    for (int x = 0; x != dot.size(); x++) {

        if (max < dot_unique_sqrt[x] / sqrt(sq_unique_sqrt[x] * sqItself_unique_sqrt[x])) {
            max = dot_unique_sqrt[x] / sqrt(sq_unique_sqrt[x] * sqItself_unique_sqrt[x]);
            maxIndex = x;
        }

    }
    m_bootstrapSupport[maxIndex]++;

    if (printFingerprint) {

        fingerprintFout.precision(3);
        fingerprintFout << "spectrum-spectrum_match_ratio: " << fixed << total / (double) (m_searchCount) << " ";
        fingerprintFout.precision(0);
        fingerprintFout << total << "/" << m_searchCount << endl;

        if (total / (double) (m_searchCount) < 0.025) {
            fingerprintFout << "WARNING: Spectrum-spectrum match ratio is too low. Results may be unreliable." << endl;
        }

        fingerprintFout << "spectral_counting(unique/total): ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout << spectralCount_unique[x] << "/" << spectralCount[x] << " ";
        }
        fingerprintFout << endl;

        fingerprintFout << "dot_product:             ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout.precision(3);
            fingerprintFout << fixed << dot[x] / sqrt(sq[x] * sqItself) << " ";
        }
        fingerprintFout << endl;

        fingerprintFout << "dot_product_sqrt:        ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout.precision(3);
            fingerprintFout << fixed << dot_sqrt[x] / sqrt(sq_sqrt[x] * sqItself_sqrt) << " ";
        }
        fingerprintFout << endl;

        fingerprintFout << "dot_product_binary:      ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout.precision(3);
            fingerprintFout << fixed << dot_binary[x] / sqrt(libSize[x] * sqItself_binary) << " ";
        }
        fingerprintFout << endl;

        fingerprintFout << "dot_product_unique:      ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout.precision(3);
            fingerprintFout << fixed << dot_unique[x] / sqrt(sq_unique[x] * sqItself_unique[x]) << " ";
        }
        fingerprintFout << endl;

        int max_index = 0;
        float max = 0;

        fingerprintFout << "dot_product_unique_sqrt: ";
        for (int x = 0; x != dot.size(); x++) {
            fingerprintFout.precision(3);
            fingerprintFout << fixed << dot_unique_sqrt[x] / sqrt(sq_unique_sqrt[x] * sqItself_unique_sqrt[x]) << " ";

            if (max < dot_unique_sqrt[x] / sqrt(sq_unique_sqrt[x] * sqItself_unique_sqrt[x])) {
                max = dot_unique_sqrt[x] / sqrt(sq_unique_sqrt[x] * sqItself_unique_sqrt[x]);
                max_index = x;
            }

        }

        fingerprintFout << endl;

        fingerprintFout.precision(0);
        for (int count = 0; count != libFP[0].size(); count++) {
            fingerprintFout << "* " << count << " * " << m_searchFingerprintIndexedByLibID[count] << " * ";
            for (int count1 = 0; count1 != libFP.size(); count1++) {
                fingerprintFout << libFP[count1][count] << " ";
            }
            fingerprintFout << endl;
        }

    }


}

void SpectraSTMzXMLSearchTask::doBootstrapping() {

    m_searchFingerprintIndexedByLibID.assign(m_lib->getMzLibIndexPtr()->getEntryCount(), 0);

    for (int count = 0; count < m_searchCount; count++) {

        unsigned int randomSearchID = (unsigned int) ((double) rand() / (double) RAND_MAX * (double) (m_searchCount));

        m_searchFingerprintIndexedByLibID[(unsigned int) (m_searchFingerprintIndexedBySearchID[randomSearchID])]++;

    }

}


void SpectraSTMzXMLSearchTask::prepareSortedSearch(unsigned int batch) {

    // display and log messages
    stringstream searchLogss;
    if (m_batchBoundaries.size() > 2) {
        if (!g_quiet) {
            cout << "BATCH " << batch + 1 << " of " << m_batchBoundaries.size() - 1 << ": Sorting query spectra in \"";
            cout << m_searchFileNames[m_batchBoundaries[batch]] << "\"..\"";
            cout << m_searchFileNames[m_batchBoundaries[batch + 1] - 1];
            cout << "\" by precursor m/z before searching...";
            cout.flush();
        }

        searchLogss << "BATCH #" << batch + 1 << " - Sorted query spectra in \"";
        searchLogss << m_searchFileNames[m_batchBoundaries[batch]] << "\"..\"";
        searchLogss << m_searchFileNames[m_batchBoundaries[batch + 1] - 1];
        searchLogss << "\" by precursor m/z";
        g_log->log("MZXML SEARCH", searchLogss.str());

    } else {
        // only one batch
        if (!g_quiet) {
            cout << "Sorting query spectra in all mzXML files by precursor m/z before searching...";
            cout.flush();
        }

        searchLogss << "Sorted query spectra in \"";
        searchLogss << m_searchFileNames[0] << "\"";
        if (m_searchFileNames.size() > 1) {
            searchLogss << "..\"" << m_searchFileNames[m_searchFileNames.size() - 1];
            searchLogss << "\"";
        }
        searchLogss << " by precursor m/z";
        g_log->log("MZXML SEARCH", searchLogss.str());
    }

    m_numScansInFile = 0;
    m_numNotSelectedInFile = 0;
    m_numMissingInFile = 0;
    m_numMS1InFile = 0;
    m_numFailedFilterInFile = 0;
    m_numSearchedInFile = 0;
    m_numLikelyGoodInFile = 0;

    m_scans.clear();

    // open all mzXML files with CRAMP
    for (unsigned int n = m_batchBoundaries[batch]; n < m_batchBoundaries[batch + 1]; n++) {

        cRamp *cramp = new cRamp(m_searchFileNames[n].c_str());
        if (!cramp->OK()) {
            g_log->error("MZXML SEARCH", "Cannot open file \"" + m_searchFileNames[n] + "\". File skipped.");
            delete (cramp);
            continue;
        }

        rampInstrumentInfo *instr = cramp->getInstrumentInfo();
        if (instr) {
            m_outputs[n]->setInstrInfo(instr);
        }

        FileName fn;
        parseFileName(m_searchFileNames[n], fn);

        // m_files is a vector of (FileName, cRamp*)
        m_files[n].first = fn;
        m_files[n].second = cramp;

        // Read the number of scans
        rampRunInfo *runInfo = cramp->getRunInfo();

        if (!runInfo) {
            // probably an empty file...
            g_log->error("MZXML SEARCH", "Cannot read run info from \"" + m_searchFileNames[n] + "\". File skipped.");
            continue;
        }

        int numScans = cramp->getLastScan();
        delete (runInfo);

        m_numScansInFile += numScans;

        // Read all scan headers into memory (excluding the peak lists to save memory)
        for (int k = 1; k <= numScans; k++) {

            // Filter out all scans not in selected list (in this case,
            // the selected list contains a list of scan numbers as strings
            if (!m_searchAll && !isInSelectedList(SpectraSTQuery::constructQueryTPPStyle(fn.name, k, k, 0))) {
                m_numNotSelectedInFile++;
                continue;
            }
            rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);

            if (!scanInfo || scanInfo->m_data.acquisitionNum != k) {
                m_numMissingInFile++;
                // the middle predicate is to deal with the case where RAMP returns a bogus scan when
                // given a nonexistent scan number -- this should become unnecessary eventually if cramp becomes smart enough
                if (scanInfo) delete scanInfo;
                continue;
            }

            if (scanInfo->m_data.msLevel == 1) {
                m_numMS1InFile++;
                delete (scanInfo);
                continue;
            }


            pair<unsigned int, rampScanInfo *> ms;
            ms.first = n;
            ms.second = scanInfo;
            m_scans.push_back(ms);

        }
    }

    // sort all the MS2 scans by precursor m/z
    sort(m_scans.begin(), m_scans.end(), SpectraSTMzXMLSearchTask::sortRampScanInfoPtrsByPrecursorMzAsc);

    // display DONE sorting message
    if (!g_quiet) {
        cout << "DONE!" << endl;
    }


}


// searchOne - search one spectrum, specified by the cRamp object that points to that mzXML file,
// and a rampScanInfo object that points to that scan.
void SpectraSTMzXMLSearchTask::searchOne(unsigned int fileIndex, rampScanInfo *scanInfo) {

    cRamp *cramp = m_files[fileIndex].second;

    // Go back to the mzxml file and get the peaks using Ramp
    rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);
    if (!peaks) {
        m_numFailedFilterInFile++;
        return;
    }

    int peakCount = peaks->getPeakCount();
    double precursorMz = scanInfo->m_data.precursorMZ;
    int precursorCharge = scanInfo->m_data.precursorCharge;
    if (precursorCharge < 1) precursorCharge = 0;

    string fragType("");
    if (scanInfo->m_data.activationMethod) fragType = scanInfo->m_data.activationMethod;

    // create the peak list and read the peaks one-by-one
    SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, precursorCharge, peakCount, false, fragType);
    peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);

    // construct the query
    string prefix = m_files[fileIndex].first.name;
    int scanNum = scanInfo->m_data.acquisitionNum;
    SpectraSTQuery *query = new SpectraSTQuery(prefix, scanNum, scanNum, precursorMz, precursorCharge, "", peakList);

    query->setRetentionTime(scanInfo->getRetentionTimeSeconds());

    for (int j = 0; j < peakCount; j++) {
        double mz = peaks->getPeak(j)->mz;
        float intensity = (float) (peaks->getPeak(j)->intensity);
        peakList->insertForSearch(mz, intensity, "");
    }

    delete peaks;


    // see if peak list passes filter; if not, ignore
    if (!(peakList->passFilter(m_params))) {
        m_numFailedFilterInFile++;
        delete query;
        return;
    }

/*
  /*
  if (m_params.filterMaxPeaksUsed < 50) {
    double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange, true, true);
  } else {
    double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange);
  }

  */

    int possibleChargeCount = 0;
    int ch = 1;
    while (possibleChargeCount < scanInfo->m_data.numPossibleCharges) {
        while (ch < 8 && !(scanInfo->m_data.possibleChargesArray[ch])) ch++;
        if (ch >= 8) break;
        query->addPossibleCharge(ch);
        possibleChargeCount++;
        ch++;
    }

    // create the Search object and search!
    SpectraSTSearch *s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
    s->search(m_lib);
    m_numSearchedInFile++;
    if (s->isLikelyGood()) {
        m_numLikelyGoodInFile++;
    }

    m_searchCount++; // counter in parent class SpectraSTSearchTask; counting all searches in search task

    if (!m_params.printFingerprintingSummary.empty()) {
        if (s->isLikelyGood()) {
            //1. spectral counts
            m_searchFingerprintIndexedByLibID[s->getTopHit()->getLibId()]++;
            m_searchFingerprintIndexedBySearchID[m_searchCount - 1] = s->getTopHit()->getLibId();


            //    cerr << "libid: " << candidates[0]->getEntry()->getLibId() << " query: " << queryName << " dot: " << candidates[0]->getSimScoresRef().dot
            //    << " delta: " << candidates[0]->getSimScoresRef().delta << " dotBias: " << candidates[0]->getSimScoresRef().dotBias
            //    << " mzDiff: " << candidates[0]->getSimScoresRef().precursorMzDiff << endl;

            //2. precursorIntensity
            //    m_searchFingerprintIndexedByLibID[candidates[0]->getEntry()->getLibId()] += query->getPrecursorIntensity();
            // log_transform
            //if (query->getPrecursorIntensity()!=0){m_searchFingerprintIndexedByLibID[candidates[0]->getEntry()->getLibId()] += log2(query->getPrecursorIntensity());}
        } else {
            m_searchFingerprintIndexedBySearchID[m_searchCount - 1] = m_lib->getMzLibIndexPtr()->getEntryCount() + 1;
        }
    }

    m_searchTaskStats.processSearch(s);

    s->print();
    delete s;

}

// sortRampScanInfoPtrsByPrecursorMzAsc - comparison function used by sort() to sort rampScanInfo pointers
// by precursor m/z
bool SpectraSTMzXMLSearchTask::sortRampScanInfoPtrsByPrecursorMzAsc(pair<unsigned int, rampScanInfo *> a,
                                                                    pair<unsigned int, rampScanInfo *> b) {

    return (a.second->m_data.precursorMZ < b.second->m_data.precursorMZ);
}

