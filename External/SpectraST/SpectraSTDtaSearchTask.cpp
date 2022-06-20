#include "SpectraSTDtaSearchTask.hpp"
#include "SpectraSTConstants.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTPeakList.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
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

/* Class: SpectraSTDtaSearchTask
 * 
 * Subclass of SpectraSTSearchTask that handles .msp file types
 * 
 */

extern bool g_quiet;
extern bool g_verbose;
extern SpectraSTLog *g_log;

// constructor
SpectraSTDtaSearchTask::SpectraSTDtaSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                               SpectraSTLib *lib) :
        SpectraSTSearchTask(searchFileNames, params, lib) {
}

// destructor
SpectraSTDtaSearchTask::~SpectraSTDtaSearchTask() {
}

// search - run the searches
void SpectraSTDtaSearchTask::search() {

    for (unsigned int n = 0; n < (unsigned int) m_searchFileNames.size(); n++) {
        searchOneFile(n);
    }
    stringstream logss;
    logss << "Searched " << m_searchCount << " out of " << m_searchFileNames.size() << " DTA files. " << m_numLikelyGood
          << " likely good.";
    g_log->log("SEARCH", logss.str());

    m_searchTaskStats.logStats();
}

// searchOneFile - search one dta file
void SpectraSTDtaSearchTask::searchOneFile(unsigned int fileIndex) {

    string searchFileName(m_searchFileNames[fileIndex]);

    ifstream fin;
    if (!myFileOpen(fin, searchFileName)) {
        g_log->error("SEARCH",
                     "Cannot open DTA file \"" + m_searchFileNames[fileIndex] + "\" for reading. File skipped.");
        return;
    }

    m_outputs[fileIndex]->openFile();
    m_outputs[fileIndex]->printHeader();

    FileName fn;
    parseFileName(searchFileName, fn);
    string name = fn.name;

    string line("");

    if (g_verbose) {
        cout << "Searching query spectrum in file \"" << searchFileName << "\" ..." << endl;
    }

    if (!nextLine(fin, line, "", "")) {
        // nothing in the file!
        return;
    }

    // parse first line, supposed to be <precursorMH> <charge>
    string::size_type pos = 0;
    double precursorMH = atof(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());
    int charge = 0;
    string chargeStr = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
    if (!chargeStr.empty()) {
        charge = atoi(chargeStr.c_str());
    }

    // special case, if no charge or zero charge is given, then the first number is assumed to be the precursorMZ
    double precursorMz = 0.0;
    if (charge == 0) {
        precursorMz = precursorMH;
    } else {
        precursorMz = (precursorMH + (double) (charge - 1)) / (double) charge;
    }

    SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, charge);
    peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);

    SpectraSTQuery *query = new SpectraSTQuery(name, precursorMz, charge, "", peakList);

    while (nextLine(fin, line, "", "")) {
        pos = 0;
        double mz = atof((nextToken(line, 0, pos)).c_str());
        float intensity = atof((nextToken(line, pos, pos)).c_str());
        peakList->insertForSearch(mz, intensity, "");

    }

    if (!peakList->passFilter(m_params)) {
        // Bad spectra, ignore
        delete query;
        m_outputs[fileIndex]->printAbortedQuery(name, "BAD_SPECTRUM_NOT_SEARCHED");

    } else {

        // double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange);

        // create the search based on what is read, then search
        SpectraSTSearch *s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
        s->search(m_lib);

        m_searchCount++; // counting searches in all dta files

        if (s->isLikelyGood()) {
            m_numLikelyGood++;
        }

        m_searchTaskStats.processSearch(s);

        // print search result
        s->print();
        delete s;
    }

    m_outputs[fileIndex]->printFooter();
    m_outputs[fileIndex]->closeFile();
}

