#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include "SpectraSTPeptideLibIndex.hpp"
#include "SpectraSTMzLibIndex.hpp"
#include "SpectraSTLibEntry.hpp"
#include "SpectraSTLog.hpp"
#include "Peptide.hpp"
#include "FileUtils.hpp"

#define NORMALCELLCOLOR   "#FFDDDD"
#define HEADERCELLCOLOR   "#42D4FD"

/*

Program       : Lib2HTML
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

/* Program: Lib2HTML
 * 
 * A standalone program to parse a SpectraST library and display it as a web page. Provides hyperlinks for the user
 * to view spectra interactively.
 */

using namespace std;

static void printUsage();

static void
writeHTMLFile(SpectraSTPeptideLibIndex *pi, string splibFileName, unsigned int maxNumEntries, string &plotspectrastCGI,
              bool moreInfo = false, SpectraSTMzLibIndex *mi = NULL);

static void writeHTMLEntries(ofstream &fout, string &libBaseName, string &plotspectrastCGI,
                             vector<SpectraSTLibEntry *> &entries, unsigned int maxNumEntries, bool isUnique,
                             bool moreInfo);

static void
writeHTMLHeader(ofstream &fout, string splibFileName, unsigned int maxNumEntries, bool isUnique, bool moreInfo);

static void writeHTMLFooter(ofstream &fout);

static bool sortProbEntriesDesc(pair<double, SpectraSTLibEntry *> a, pair<double, SpectraSTLibEntry *> b);


SpectraSTLog *g_log; // since we are using SpectraST code, we need to have a log file.

// main
int main(int argc, char **argv) {

    if (argc < 2) {
        printUsage();
        return (1);
    }

    unsigned int maxNumEntries = 10;
    bool moreInfo = false;
    bool sortByMz = false;
    string plotspectrastCGI(getCGIPath() + "plotspectrast.cgi");

    string splibFileName = argv[argc - 1];

    // convert relative paths to absolute paths
    makeFullPath(splibFileName);

    FileName fn;
    parseFileName(splibFileName, fn);

    if (fn.ext != ".splib") {
        printUsage();
        return (1);
    }

    for (int ar = 1; ar < argc - 1; ar++) {
        if (argv[ar][0] != '-') {
            printUsage();
        }
        if (argv[ar][1] == 'N') {
            maxNumEntries = atoi(argv[ar] + 2);
            if (maxNumEntries < 1) maxNumEntries = 1;
        } else if (argv[ar][1] == 'V') {
            moreInfo = true;
        } else if (argv[ar][1] == 'm') {
            sortByMz = true;
        } else if (argv[ar][1] == 'P') {
            plotspectrastCGI = argv[ar] + 2;
        }
    }

    g_log = new SpectraSTLog("Lib2HTML.log");

    ifstream splibFin;
    if (!myFileOpen(splibFin, splibFileName, true)) {
        cerr << "Cannot open SPLIB file \"" + splibFileName + "\" for reading. Exiting..." << endl;
        return (1);
    }

    // peek at the library to see if it's in binary format or not.
    bool binary = true;
    char firstChar = splibFin.peek();
    if (firstChar == '#' || firstChar == 'N') {
        binary = false;
    }

    // read the peptide index
    SpectraSTPeptideLibIndex *pi = new SpectraSTPeptideLibIndex(fn.path + fn.name + ".pepidx", &splibFin, binary);

    SpectraSTMzLibIndex *mi = NULL;

    if (sortByMz) {
        mi = new SpectraSTMzLibIndex(fn.path + fn.name + ".spidx", &splibFin, 1.0, binary);
    }

    // write HTML file
    writeHTMLFile(pi, splibFileName, maxNumEntries, plotspectrastCGI, moreInfo, mi);

    if (mi) {
        delete (mi);
    }

    delete (pi);
}

// printUsage - prints out usage information
void printUsage() {
    cout << "Lib2HTML v.1.0 - Henry Lam, March 2007" << endl;
    cout << "Creates an HTML file for viewing a SpectraST library over HTTP." << endl;
    cout << "Usage: Lib2HTML [options] <library file name>.splib" << endl;
    cout << endl;
    cout << "Options: -P<URL to plotspectrast.cgi> - in case it is not where I expect it to be." << endl;
    cout << "         -N<Maximum number of spectra linked per peptide ion> - default is 10." << endl;
    cout << "         -V - verbose; provides more information for each entry." << endl;
    cout << "         -m - sort by precursor m/z. default is off, i.e. sort by peptide sequence." << endl;
    cout << endl;
    cout
            << "Note: If the environment variable WEBSERVER_ROOT is not set, include the path from the webserver root to library file."
            << endl;
    cout
            << "      e.g. If library is at /wwwroot/lib/foo.splib and webserver root is /wwwroot/, supply /lib/foo.splib as argument. Don't forget the preceding '/'."
            << endl;
    cout << endl;

}

// writeHTMLFile - the meat of the program. walks down the peptide index, parses out entries, and serves them up in HTML
void
writeHTMLFile(SpectraSTPeptideLibIndex *pi, string splibFileName, unsigned int maxNumEntries, string &plotspectrastCGI,
              bool moreInfo, SpectraSTMzLibIndex *mi) {

    FileName fn;
    parseFileName(splibFileName, fn);

    // the output is always of the same name and path as the .splib, but with .html extension.
    string htmlFileName(fn.path + fn.name + ".html");
    string libBaseName = fn.path + fn.name; // should be absolute paths

    std::string lbn(libBaseName);
    translate_absolute_filesystem_path_to_relative_webserver_root_path(lbn);
    string libBaseNameFromWebserverRoot(lbn);

    ofstream fout;
    if (!myFileOpen(fout, htmlFileName)) {
        cerr << "Cannot open file \"" << htmlFileName << "\" for writing HTML-style peptide index. Exiting..." << endl;
        exit(1);
    }

    // check if it's a unique library (i.e. each entry having only one replicate)
    bool isUnique = pi->isUniqueLibrary();

    writeHTMLHeader(fout, splibFileName, maxNumEntries, isUnique, moreInfo);

    if (!mi) {

        string spep;
        vector<string> subkeys;


        while (pi->nextPeptide(spep, subkeys)) {
            // loop over all peptide sequences

            for (vector<string>::iterator sk = subkeys.begin(); sk != subkeys.end(); sk++) {
                // loop over all peptide ions derived from that sequence (different charge states, different modifications)

                // retrieves all entries of this peptide ion
                vector<SpectraSTLibEntry *> entries;
                pi->retrieve(entries, spep, *sk);

                // write out the entries
                writeHTMLEntries(fout, libBaseName, plotspectrastCGI, entries, maxNumEntries, isUnique, moreInfo);

                // done, deletes the entries.
                for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {
                    delete (*en);
                }
            }
            subkeys.clear();
        }

    } else {

        // use m/z index instead

        SpectraSTLibEntry *entry = NULL;
        vector<SpectraSTLibEntry *> entries;
        string lastEntry("");

        while ((entry = mi->nextEntry())) {

            if (lastEntry.empty()) {

                entries.push_back(entry);
                lastEntry = entry->getFullName();

            } else if (lastEntry == entry->getFullName()) {

                entries.push_back(entry);

            } else {

                // write out the entries
                writeHTMLEntries(fout, libBaseName, plotspectrastCGI, entries, maxNumEntries, isUnique, moreInfo);

                for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {
                    delete (*en);
                }

                entries.clear();
                entries.push_back(entry);
                lastEntry = entry->getFullName();

            }
        }

        // write out the last entries
        writeHTMLEntries(fout, libBaseName, plotspectrastCGI, entries, maxNumEntries, isUnique, moreInfo);

        for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {
            delete (*en);
        }


    }

    writeHTMLFooter(fout);

}

// writeHTMLEntries - prints out an HTML table row for each entry
void writeHTMLEntries(ofstream &fout, string &libBaseName, string &plotspectrastCGI,
                      vector<SpectraSTLibEntry *> &entries, unsigned int maxNumEntries, bool isUnique, bool moreInfo) {

    fout << "<TR BGCOLOR=\"" << NORMALCELLCOLOR << "\" ALIGN=\"LEFT\">" << endl;

    if (isUnique) {

        // parse out number of replicates
        string nrepsStr("");
        int numUsedReps = 1;
        int numTotalReps = 1;
        if (entries[0]->getOneComment("Nreps", nrepsStr)) {
            string::size_type slashPos = 0;
            numUsedReps = atoi(nextToken(nrepsStr, 0, slashPos, "/").c_str());
            numTotalReps = atoi(nextToken(nrepsStr, slashPos + 1, slashPos).c_str());
        }

        //parse out probability
        double prob = entries[0]->getProb();

        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries[0]->getFullName() << "</TT></TD>" << endl;
        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries[0]->getPrecursorMz() << "</TT></TD>"
             << endl;
        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries[0]->getStatus() << "</TT></TD>" << endl;

        int proteinCount = 0;
        string protein = entries[0]->getFirstProtein(proteinCount);
        if (protein.length() > 20) protein = protein.substr(0, 20);
        if (proteinCount > 1) {
            stringstream proteinss;
            proteinss << " (+" << proteinCount - 1 << ")";
            protein += proteinss.str();
        }


        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << protein << "</TT></TD>" << endl;

        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << prob << "</TT></TD>" << endl;
        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << numUsedReps << "</TT></TD>" << endl;

        if (moreInfo) {

            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << numTotalReps << "</TT></TD>" << endl;

            string spec("");
            bool good = entries[0]->getOneComment("Spec", spec);

            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << spec << "</TT></TD>" << endl;


            string fracUnassigned("");
            if (entries[0]->getOneComment("FracUnassigned", fracUnassigned)) {
                string::size_type semicolonPos = fracUnassigned.find(';');
                fracUnassigned = nextToken(fracUnassigned, semicolonPos + 1, semicolonPos, "; \t\r\n");
            }

            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << fracUnassigned << "</TT></TD>" << endl;

            string xrea("");
            good = entries[0]->getOneComment("Xrea", xrea);

            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << xrea << "</TT></TD>" << endl;
        }

        // link to plotspectrast.cgi to display the spectrum
        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>";
        fout << "<A HREF=\"" << plotspectrastCGI << "?LibFile=" << libBaseName << ".splib";
        fout << "&LibFileOffset=" << entries[0]->getLibFileOffset();
        fout << "&QueryFile=" << libBaseName << "_" << entries[0]->getLibId() << ".none\">";
        fout << entries[0]->getLibId() << "</A>";

        fout << "</TT></TD>" << endl;

    } else {
        // not unique, spectrum links will be listed as columns following the peptide ion

        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries[0]->getFullName() << "</TT></TD>" << endl;
        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries[0]->getPrecursorMz() << "</TT></TD>"
             << endl;

        int proteinCount = 0;
        string protein = entries[0]->getFirstProtein(proteinCount);
        if (protein.length() > 20) protein = protein.substr(0, 20);
        if (proteinCount > 1) {
            stringstream proteinss;
            proteinss << " (+" << proteinCount - 1 << ")";
            protein += proteinss.str();
        }

        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << protein << "</TT></TD>" << endl;

        fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entries.size() << "</TT></TD>" << endl;

        // add all the entries into a vector of (prob, SpectraSTLibEntry*) pairs, this is for sorting
        vector<pair<double, SpectraSTLibEntry *> > probEntries;
        double maxXrea = 0.0;
        double minXrea = 1.0;
        for (vector<SpectraSTLibEntry *>::iterator en = entries.begin(); en != entries.end(); en++) {

            string xreaStr("");
            if ((*en)->getOneComment("Xrea", xreaStr)) {
                double xrea = atof(xreaStr.c_str());
                if (xrea < minXrea) minXrea = xrea;
                if (xrea > maxXrea) maxXrea = xrea;
            }

            double prob = (*en)->getProb();
            pair<double, SpectraSTLibEntry *> p;
            p.first = prob;
            p.second = (*en);
            probEntries.push_back(p);
        }

        if (moreInfo) {
            stringstream xreass;
            if (maxXrea < minXrea) {
                xreass << "--";
            } else {
                xreass.precision(3);
                xreass << fixed << minXrea << "-" << fixed << maxXrea;
            }
            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << xreass.str() << "</TT></TD>" << endl;
        }

        // sort the replicates by probabilities, higher-probability replicates will be linked first if there're not enough spaces.
        std::sort(probEntries.begin(), probEntries.end(), sortProbEntriesDesc);

        unsigned int numPrinted = 0;
        for (vector<pair<double, SpectraSTLibEntry *> >::iterator pe = probEntries.begin();
             pe != probEntries.end() && numPrinted <= maxNumEntries; pe++, numPrinted++) {

            // link to plotspectrast.cgi to display the spectrum
            fout << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>";
            fout << "<A HREF=\"" << plotspectrastCGI << "?LibFile=" << libBaseName << ".splib";
            fout << "&LibFileOffset=" << pe->second->getLibFileOffset();
            fout << "&QueryFile=" << libBaseName << "_" << pe->second->getLibId() << ".none\">";
            fout << pe->second->getLibId() << "(" << pe->first << ")</A>";
            fout << "</TT></TD>" << endl;
        }
    }
    fout << "</TR>" << endl;
}

// writeHTMLHeader - writes the HTML header, plus opens the table and prints the table headers too.
void writeHTMLHeader(ofstream &fout, string splibFileName, unsigned int maxNumEntries, bool isUnique, bool moreInfo) {

    fout << "Content-type: text/html\n\n";
    fout << "<HTML>" << endl;
    fout << "  <HEAD>" << endl;
    fout << "    <TITLE>" << "SpectraST Library Viewer (v.1.0, by Henry Lam, ISB, 2007): " << splibFileName
         << "</TITLE>" << endl;
    fout << "  </HEAD>" << endl;
    fout << endl;
    fout << "<BODY BGCOLOR=\"#EEEEEE\" OnLoad=\"self.focus();\">" << endl;


    fout << "<TABLE BORDER=0 CELLPADDING=\"2\">" << endl;
    fout << "<TBODY>";
    fout << "<TR BGCOLOR=\"" << HEADERCELLCOLOR << "\" ALIGN=\"CENTER\">";

    if (isUnique) {
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "PeptideIon" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "PrecMz" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Status" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Protein" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Prob" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Nreps" << "</TT></TH>" << endl;
        if (moreInfo) {
            fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "NrepsTot" << "</TT></TH>" << endl;

            fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Spec" << "</TT></TH>" << endl;
            fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "FrUnas20" << "</TT></TH>" << endl;
            fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Xrea" << "</TT></TH>" << endl;
        }
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "LibID" << "</TT></TH>" << endl;
    }
    else {
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "PeptideIon" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "PrecMz" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Protein" << "</TT></TH>" << endl;
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Nreps" << "</TT></TH>" << endl;
        if (moreInfo) {
            fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "XreaRange" << "</TT></TH>" << endl;
        }
        fout << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\" COLSPAN=" << maxNumEntries << "><TT>" << "LibID(s) (Prob)"
             << "</TT></TH>" << endl;
    }


    fout << "</TR>" << endl;

}

// writeHTMLFooter - closes the table and the HTML file.
void writeHTMLFooter(ofstream &fout) {

    fout << "</TABLE>" << endl;
    fout << "</BODY>" << endl;
    fout << "</HTML>" << endl;
}

// sortProbEntriesDesc - comparator for sorting
bool sortProbEntriesDesc(pair<double, SpectraSTLibEntry *> a, pair<double, SpectraSTLibEntry *> b) {
    return (a.first > b.first);
}
