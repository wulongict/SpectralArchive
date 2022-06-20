#include "SpectraSTMs2LibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>
#include <stdlib.h>

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

/* Class: SpectraSTMs2LibImporter
 * 
 * Implements a library importer for the .ms2 file format (used by BiblioSpec).
 * Note that there is no guarantee that BiblioSpec libraries will work well with SpectraST!
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTMs2LibImporter::SpectraSTMs2LibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                 SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params) {}

// destructor
SpectraSTMs2LibImporter::~SpectraSTMs2LibImporter() {

}

// import - prints the preamble, then loops over all files and import them one by one
void SpectraSTMs2LibImporter::import() {

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        string fullName(*i);
        makeFullPath(fullName);
        string quoted("\"" + fullName + "\"");
        string desc = m_params.constructDescrStr(quoted, ".ms2");
        m_preamble.push_back(desc);
    }

    m_lib->writePreamble(m_preamble);

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        readFromFile(*i);
    }
}

// readFromCurFile - reads one .ms2 file
void SpectraSTMs2LibImporter::readFromFile(string &impFileName) {


    ifstream fin;
    if (!myFileOpen(fin, impFileName)) {
        g_log->error("CREATE", "Cannot open .ms2 file \"" + impFileName + "\" for reading. File skipped.");
        return;
    }

    g_log->log("MS2 IMPORT", "Importing .ms2 file \"" + impFileName + "\".");

    if (g_verbose) {
        cout << "\nImporting spectra from .ms2 library file..." << endl;
    }

    // start the progress count
    ProgressCount pc(!g_quiet && !g_verbose, 500, 0);
    pc.start("\nImporting spectra from .ms2 library file");

    string line("");
    unsigned int scan1 = 0;
    unsigned int scan2 = 0;
    double precursorMz = 0.0;
    int charge = 1;
    double mw = 0.0;

    string dummy("");
    string seq("");
    string modifiedSeq("");

    SpectraSTPeakList *peakList = NULL;

    // looks like their cysteines are actually C[160]'s
    // Peptide::addModTokenToTables("C", "Carbamidomethyl");

    while (nextLine(fin, line)) {

        if (line.empty()) continue;

        string::size_type pos = 0;

        if (line[0] == 'H') {
            continue;

        } else if (line[0] == 'S') {
            if (peakList) {
                if (peakList->getNumPeaks() == 0 || seq.empty()) {
                    delete peakList;
                } else {

                    pc.increment();
                    if (!(modifiedSeq.empty())) seq = modifiedSeq;

                    Peptide *pep = createPeptide(seq, charge, "", "", "MS2");

                    if (pep) {

                        if (m_params.setDeamidatedNXST) {
                            setDeamidatedNXST(pep);
                        }

                        stringstream commentss;
                        commentss << "Fullname=" << pep->interactStyleFullWithCharge();
                        commentss << " Prob=1.0000";
                        commentss << " ScanNum=" << scan1 << '.' << scan2;
                        commentss << " Spec=Raw";

                        SpectraSTLibEntry *entry = new SpectraSTLibEntry(pep, commentss.str(), "Normal", peakList);


                        if (insertOneEntry(entry, "MS2")) {
                            if (g_verbose) {
                                cout << "Importing record " << m_count << ": " << pep->interactStyleWithCharge()
                                     << endl;
                            }
                            m_count++;
                        }
                        delete (entry);

                    }
                }

            }

            seq = "";
            modifiedSeq = "";
            peakList = NULL;
            scan1 = atoi(nextToken(line, 1, pos, " \t\r\n", " \t\r\n").c_str());
            scan2 = atoi(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());
            precursorMz = atof(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());

        } else if (line[0] == 'Z') {

            charge = atoi(nextToken(line, 1, pos, " \t\r\n", " \t\r\n").c_str());
            mw = atof(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());

        } else if (line[0] == 'D') {
            dummy = nextToken(line, 1, pos, "\t\r\n", " \t\r\n");
            if (dummy == "seq") {
                seq = nextToken(line, pos, pos, "\r\n", " \t\r\n");
            } else if (dummy == "modified seq") {
                modifiedSeq = nextToken(line, pos, pos, "\r\n", " \t\r\n");
            }

        } else {
            // should be a peak
            if (!peakList) {
                peakList = new SpectraSTPeakList(precursorMz, 0);
                peakList->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);
            }
            double mz = atof(nextToken(line, 0, pos, " \t\r\n", " \t\r\n").c_str());
            float intensity = atof(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());
            peakList->insert(mz, intensity, "", "");

        }


    }

    // finish last record
    if (peakList) {
        if (peakList->getNumPeaks() == 0 || seq.empty()) {
            delete peakList;
        } else {

            pc.increment();
            if (!(modifiedSeq.empty())) seq = modifiedSeq;

            Peptide *pep = createPeptide(seq, charge, "", "", "MS2");

            if (pep) {

                if (m_params.setDeamidatedNXST) {
                    setDeamidatedNXST(pep);
                }

                stringstream commentss;
                commentss << "Fullname=" << pep->interactStyleFullWithCharge();
                commentss << " Prob=1.0000";
                commentss << " ScanNum=" << scan1 << '.' << scan2;
                commentss << " Spec=Raw";

                SpectraSTLibEntry *entry = new SpectraSTLibEntry(pep, commentss.str(), "Normal", peakList);

                if (insertOneEntry(entry, "MS2")) {
                    if (g_verbose) {
                        cout << "Importing record " << m_count << ": " << pep->interactStyleWithCharge() << endl;
                    }
                    m_count++;
                }
                delete (entry);

            }
        }

    }


    pc.done();
}





