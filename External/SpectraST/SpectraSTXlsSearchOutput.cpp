#include "SpectraSTXlsSearchOutput.hpp"
#include "Peptide.hpp"
#include "FileUtils.hpp"
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

/* Class: SpectraSTXlsSearchOutput
 * 
 * Outputter to .xls format
 *    
 */


// constructor
SpectraSTXlsSearchOutput::SpectraSTXlsSearchOutput(string outFullFileName, string inputFileExt,
                                                   SpectraSTSearchParams &searchParams, bool printHeader) :
        SpectraSTSearchOutput(outFullFileName, inputFileExt, searchParams),
        m_printHeader(printHeader) {

}

// destructor
SpectraSTXlsSearchOutput::~SpectraSTXlsSearchOutput() {

}

void SpectraSTXlsSearchOutput::printHeader() {
    if (!m_fout) openFile();

    if (m_printHeader) {
        (*m_fout) << "### Query" << '\t';
        (*m_fout) << "Rk" << '\t';
        (*m_fout) << "ID" << '\t';
        SpectraSTSimScores::printHeaderTabDelimited(*m_fout);
        (*m_fout) << "Status" << '\t';
        (*m_fout) << "Inst" << '\t';
        (*m_fout) << "Spec" << '\t';
        (*m_fout) << "#Pr" << '\t';
        (*m_fout) << "Proteins" << '\t';
        (*m_fout) << "LibFileOffset" << endl;
    }
}

// printStartQuery - prints the query information
void
SpectraSTXlsSearchOutput::printStartQuery(string query, double precursorMz, int assumedCharge, double retentionTime) {

    (*m_fout) << query << '\t';

    m_query = query;

}

// printHit - prints the hit
void SpectraSTXlsSearchOutput::printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry,
                                        SpectraSTSimScores &simScores) {

    if (hitRank != 1) {
        // lower hits, print -- to take the place of the query name for nicer output
        (*m_fout) << "-----------" << '\t';
    }

    if (!entry) {
        (*m_fout) << "0" << '\t' << "NO_MATCH" << '\t';
    } else {
        (*m_fout) << hitRank << '\t' << entry->getName() << '\t';

        simScores.printTabDelimited((*m_fout));

        (*m_fout) << entry->getStatus() << '\t';

        string inst("");
        if (entry->getOneComment("Inst", inst)) {
            (*m_fout) << inst << '\t';
        } else {
            (*m_fout) << "Unk" << '\t';
        }

        string spectrumType("");
        if (entry->getOneComment("Spec", spectrumType)) {
            (*m_fout) << spectrumType.substr(0, spectrumType.length() < 3 ? spectrumType.length() : 3) << '\t';
        } else {
            (*m_fout) << "Unk" << '\t';
        }

        string protein("");
        stringstream proteinss;
        string sample("");
        string::size_type pos = 0;
        if (entry->getOneComment("Protein", protein)) {
            if (protein.find('/') != string::npos) {
                int proteinCount = atoi(nextToken(protein, 0, pos, "/\t\r\n").c_str());
                proteinss << proteinCount << '\t';
                bool isFirstProtein = true;
                while (pos < protein.length()) {
                    if (!isFirstProtein) proteinss << ';';
                    proteinss << nextToken(protein, pos, pos, "/\t\r\n", "/");
                    isFirstProtein = false;
                }
            } else {
                proteinss << '1' << '\t' << protein;
            }

        } else if (entry->getOneComment("Sample", sample)) {
            // if no protein field is given, use the sample field. (this is useful for unidentified spectra)
            if (sample.find('/') != string::npos) {
                int sampleCount = atoi(nextToken(sample, 0, pos, "/\t\r\n").c_str());
                proteinss << sampleCount << '\t';
                bool isFirstSample = true;
                while (pos < sample.length()) {
                    if (!isFirstSample) proteinss << ';';
                    proteinss << nextToken(sample, pos, pos, "/\t\r\n", "/");
                    isFirstSample = false;
                }
            }

        } else {
            proteinss << '0';
        }

        (*m_fout) << proteinss.str() << '\t';

        (*m_fout) << entry->getLibFileOffset();
    }
    (*m_fout) << endl;

}

// printAbortedQuery - displays a message for a query that is not searched (used for dta input, for example)
void SpectraSTXlsSearchOutput::printAbortedQuery(string query, string message) {

    (*m_fout) << query << '\t' << message << endl;
}

