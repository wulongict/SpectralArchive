#include "SpectraSTHtmlSearchOutput.hpp"
#include "Peptide.hpp"
#include "FileUtils.hpp"
#include <iostream>
#include <sstream>

#define NORMALCELLCOLOR   "#FFDDDD"
#define HEADERCELLCOLOR   "#42D4FD"

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

/* Class: SpectraSTHtmlSearchOutput
 * 
 * Outputter to .html format
 *    
 */

string SpectraSTHtmlSearchOutput::m_plotspectrastCGI = "";

// constructor
SpectraSTHtmlSearchOutput::SpectraSTHtmlSearchOutput(string outputFileName, string inputFileExt,
                                                     SpectraSTSearchParams &searchParams) :
        SpectraSTSearchOutput(outputFileName, inputFileExt, searchParams) {

    FileName fn;
    parseFileName(m_outputFileName, fn);

    string searchFile(fn.path + fn.name + inputFileExt);
    makeFullPath(searchFile);
    m_searchFile = searchFile;

    string libFile(m_searchParams.libraryFile);
    makeFullPath(libFile);
    m_libFile = libFile;


}

// destructor
SpectraSTHtmlSearchOutput::~SpectraSTHtmlSearchOutput() {

}

void SpectraSTHtmlSearchOutput::printHeader() {
    if (!m_fout) openFile();

    // (*m_fout) << "Content-type: text/html\n\n";
    (*m_fout) << "<HTML>" << endl;
    (*m_fout) << "  <HEAD>" << endl;
    (*m_fout) << "    <TITLE>" << "SpectraST Search Results: " << m_searchFile << " against "
              << m_searchParams.libraryFile << "</TITLE>" << endl;
    (*m_fout) << "  </HEAD>" << endl;
    (*m_fout) << endl;
    (*m_fout) << "<BODY BGCOLOR=\"#EEEEEE\" OnLoad=\"self.focus();\">" << endl;

    // TODO: Put in some more info

    (*m_fout) << "<TABLE BORDER=0 CELLPADDING=\"2\">" << endl;
    (*m_fout) << "<TBODY>";
    (*m_fout) << "<TR BGCOLOR=\"" << HEADERCELLCOLOR << "\" ALIGN=\"CENTER\">";

    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Query" << "</TT></TH>" << endl;
    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "PrecMZ" << "</TT></TH>" << endl;
    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Rank" << "</TT></TH>" << endl;
    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "ID" << "</TT></TH>" << endl;

    SpectraSTSimScores::printHeaderHtml(*m_fout);

    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Status" << "</TT></TH>" << endl;
    (*m_fout) << "  <TH BGCOLOR=\"" << HEADERCELLCOLOR << "\"><TT>" << "Protein(s)" << "</TT></TH>" << endl;

    (*m_fout) << "</TR>" << endl;

}

// printStartQuery - prints the query information
void
SpectraSTHtmlSearchOutput::printStartQuery(string query, double precursorMz, int assumedCharge, double retentionTime) {

    string::size_type dotpos;
    dotpos = query.rfind('.');
    dotpos = query.rfind('.', dotpos - 1);
    dotpos = query.rfind('.', dotpos - 1);
    string startScan = nextToken(query, dotpos + 1, dotpos, ".");

    (*m_fout) << "<TR BGCOLOR=\"" << NORMALCELLCOLOR << "\" ALIGN=\"LEFT\">" << endl;

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << query << "</TT></TD>" << endl;
    (*m_fout).precision(4);
    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << fixed << precursorMz << "</TT></TD>" << endl;

    m_query = query;
    m_queryScanNum = startScan;


}

void SpectraSTHtmlSearchOutput::printFooter() {
    (*m_fout) << "</TABLE>" << endl;
    (*m_fout) << "</BODY>" << endl;
    (*m_fout) << "</HTML>" << endl;
}

// printHit - prints the hit
void SpectraSTHtmlSearchOutput::printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry,
                                         SpectraSTSimScores &simScores) {

    if (hitRank != 1) {
        // lower hits, print two empty cells to take the place of the query name and precursor m/z
        (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT></TT></TD>" << endl;
        (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT></TT></TD>" << endl;
    }

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << hitRank << "</TT></TD>" << endl;

    if (!entry) {
        (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << "NO_MATCH" << "</TT></TD>" << endl;
        (*m_fout) << "</TR>" << endl;
        return;
    }

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>";
    (*m_fout) << "<A HREF=\"" << getCGIPath() << "plotspectrast.cgi";
    (*m_fout) << "?LibFile=" << m_libFile;
    (*m_fout) << "&LibFileOffset=" << entry->getLibFileOffset();

    if (m_searchFile.substr(m_searchFile.length() - 6) == ".sptxt") {
        (*m_fout) << "&QueryFile=" << m_searchFile.substr(0, m_searchFile.length() - 6) << ".splib";
    } else {
        (*m_fout) << "&QueryFile=" << m_searchFile;
    }
    (*m_fout) << "&QueryScanNum=" << m_queryScanNum << "\">";
    (*m_fout) << entry->getFullName() << "</A></TT></TD>" << endl;

    simScores.printHtml((*m_fout));

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << entry->getStatus() << "</TT></TD>" << endl;

    int proteinCount = 0;
    string protein = entry->getFirstProtein(proteinCount);
    if (protein.length() > 40) protein = protein.substr(0, 40);
    if (proteinCount > 1) {
        stringstream proteinss;
        proteinss << " (+" << proteinCount - 1 << ")";
        protein += proteinss.str();
    }

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << protein << "</TT></TD>" << endl;

    (*m_fout) << "</TR>" << endl;

}

// printAbortedQuery - displays a message for a query that is not searched (used for dta input, for example)
void SpectraSTHtmlSearchOutput::printAbortedQuery(string query, string message) {

    (*m_fout) << "<TR BGCOLOR=\"" << NORMALCELLCOLOR << "\" ALIGN=\"LEFT\">" << endl;

    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << query << "</TT></TD>" << endl;
    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << "0.0000" << "</TT></TD>" << endl;
    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << "1" << "</TT></TD>" << endl;
    (*m_fout) << "  <TD BGCOLOR=\"" << NORMALCELLCOLOR << "\"><TT>" << "NOT_SEARCHED" << "</TT></TD>" << endl;

    (*m_fout) << "</TR>" << endl;

}

