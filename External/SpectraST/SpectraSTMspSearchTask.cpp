#include "SpectraSTMspSearchTask.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"


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

/* Class: SpectraSTMspSearchTask
 * 
 * Subclass of SpectraSTSearchTask that handles .msp file types
 * 
 */

extern bool g_quiet;
extern bool g_verbose;
extern SpectraSTLog *g_log;

// constructor
SpectraSTMspSearchTask::SpectraSTMspSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                               SpectraSTLib *lib) :
        SpectraSTSearchTask(searchFileNames, params, lib) {
}

// destructor
SpectraSTMspSearchTask::~SpectraSTMspSearchTask() {
}

// search - run the searches
void SpectraSTMspSearchTask::search() {

    for (unsigned int n = 0; n < (unsigned int) m_searchFileNames.size(); n++) {

        searchOneFile(n);
    }

    m_searchTaskStats.logStats();
}

// searchOneFile - search one msp file
void SpectraSTMspSearchTask::searchOneFile(unsigned int fileIndex) {

    string searchFileName(m_searchFileNames[fileIndex]);

    ifstream fin;
    if (!myFileOpen(fin, searchFileName)) {
        g_log->error("SEARCH",
                     "Cannot open MSP file \"" + searchFileName + " for reading query spectra. File Skipped.");
        return;
    }

    m_outputs[fileIndex]->openFile();
    m_outputs[fileIndex]->printHeader();

    string line("");
    string::size_type pos = 0;

    if (g_verbose) {
        cout << "Searching query spectra in file \"" << searchFileName << "\" ..." << endl;
    }

    ProgressCount pc(!g_quiet && !g_verbose, 200);
    string msg("Searching query spectra in file \"");
    msg += searchFileName + "\"";
    pc.start(msg);

    while (true) {

        if (line == "_EOF_") {
            // no more record
            pc.done();
            m_outputs[fileIndex]->printFooter();
            m_outputs[fileIndex]->closeFile();
            return;
        }

        // skips over all lines until the line with Name: (will stop when it reaches either "Name:" or the end-of-file)
        if (line.compare(0, 6, "Name: ") != 0) {
            while (nextLine(fin, line, "Name: ", ""));
            if (line == "_EOF_") {
                // no more record
                pc.done();
                m_outputs[fileIndex]->printFooter();
                m_outputs[fileIndex]->closeFile();
                return;
            }
        }


        // line should now start with "Name: "
        string name = nextToken(line, 5, pos, "\r\n");

        bool selected;
        if (!m_searchAll) {
            selected = isInSelectedList(name);
            // not to mess up the parsing, the following fields will still be
            // parsed even though the record is not selected to be searched.
            // the search simply will not be initiated.
        } else {
            selected = true;
        }


        double mw = 0.0;
        double precursorMz = 0.0;
        string comments("");
        int charge = 0;

        // read the rest of the header fields until Num peaks:
        while (nextLine(fin, line, "Num peaks: ", "")) {
            if (line.compare(0, 3, "MW:") == 0) {
                mw = atof((nextToken(line, 3, pos, "\r\n")).c_str());
            } else if (line.compare(0, 8, "Comment:") == 0) {
                comments = nextToken(line, 8, pos, "\r\n");
            } else if (line.compare(0, 12, "PrecursorMZ:") == 0) {
                precursorMz = atof((nextToken(line, 12, pos, "\r\n")).c_str());
            } else if (line.compare(0, 7, "Charge:") == 0) {
                charge = atoi((nextToken(line, 7, pos, "\r\n")).c_str());
            } else if (line.compare(0, 9, "NumPeaks:") == 0) {
                // reach here because this .msp probably is converted from a .sptxt,
                // hack the line (pad one char) so that the number of peaks will be correctly parsed later
                line = " " + line;
                break;
            } else if (line.compare(0, 7, "Status:") == 0 ||
                       line.compare(0, 9, "FullName:") == 0 ||
                       line.compare(0, 6, "LibID:") == 0) {
                // reach here because this .msp probably is converted from a .sptxt,
                // just ignore this line
            } else if (line == "_EOF_" || line.compare(0, 5, "Name:") == 0) {
                // reach the end unexpectedly, or see another name field before the Num peaks field.
                // ignore this incomplete record, and return
                cerr << "\nBadly formatted .msp file! Library creation truncated." << endl;
                return;
            } else {
                cerr << "Unrecognized header field. Ignored." << endl;
            }
        }

        if (line == "_EOF_") {
            // no "Num peaks:" field. ignore this incomplete record, and return
            cerr << "\nBadly formatted .msp file! Library creation truncated." << endl;
            m_outputs[fileIndex]->printFooter();
            m_outputs[fileIndex]->closeFile();

            return;
        }


        if (precursorMz < 0.0001) {
            // Precursor m/z not specified as a header field, try to find it in the comments
            string::size_type parentPos = comments.find("Parent", 0);
            if (parentPos != string::npos) {
                precursorMz = atof((nextToken(comments, parentPos + 7, parentPos, " \t\r\n")).c_str());
            } else {
                if (charge != 0) {
                    // Charge is specified. Can calculate precursor m/z
                    precursorMz = (mw + charge * 1.00728) / (double) charge;
                } else {
                    // nothing we can do, maybe the mw is really the precursor m/z
                    precursorMz = mw;
                }
            }
        }

        // line now starts with "Num peaks: "
        int numPeaks = atoi((nextToken(line, 10, pos, "\r\n")).c_str());

        SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, charge, numPeaks);
        peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);

        // hack -- if the query file is an .sptxt, then put the binaryFileOffset as
        // part of the query name. This allows the HTML output file to link to the spectrum
        string::size_type binaryFileOffsetPos = comments.find("BinaryFileOffset", 0);
        if (binaryFileOffsetPos != string::npos) {
            string binaryFileOffset = nextToken(comments, binaryFileOffsetPos + 17, binaryFileOffsetPos, " \t\r\n");
            name += "." + binaryFileOffset + "." + binaryFileOffset + ".0";
        }

        SpectraSTQuery *query = new SpectraSTQuery(name, precursorMz, charge, comments, peakList);

        while (nextLine(fin, line, "Name: ", "")) { // will stop when it reaches the next "Name:" or end-of-file

            // here are the peaks
            double mz = atof((nextToken(line, 0, pos)).c_str());
            // pos1 now stores the position of the first space after the mz
            float intensity = (float) (atof((nextToken(line, pos, pos)).c_str()));
            // pos2 now stores the position of the first space after the intensity

            // annotation has quotes around it, remove them by adding the quote char to the skipover and delimiter strings passed into nextToken
            string annotation = nextToken(line, pos, pos, "\"\r\n", "\"\t");

            string info("");

            string::size_type spacePos = annotation.find_first_of(" \t", 0);
            string::size_type dummyPos = 0;
            if (spacePos != string::npos) {
                info = nextToken(annotation, spacePos + 1, dummyPos, "\r\n", " \t");
                annotation = annotation.substr(0, spacePos);
            }
            // annotation will get an empty string if there's no annotation
            peakList->insertForSearch(mz, intensity, annotation);


        }

        if (!selected || !peakList->passFilter(m_params)) {
            // not selected to be searched, or bad spectra, ignore
            delete query;
        } else {

//      if (m_params.filterMaxPeaksUsed <= 50) {
//        double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange, true, true);
//      } else {
//        double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange);
//      }

            FileName fn;
            parseFileName(searchFileName, fn);

            // create the search based on what is read, then search
            SpectraSTSearch *s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
            s->search(m_lib);

            m_searchCount++; // counting searches in all msp files
            m_searchTaskStats.processSearch(s);

            pc.increment();

            // print search result
            s->print();
            delete s;
        }

    }
    m_outputs[fileIndex]->printFooter();
    m_outputs[fileIndex]->closeFile();

}

