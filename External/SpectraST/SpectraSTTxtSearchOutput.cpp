#include "SpectraSTTxtSearchOutput.hpp"
#include "SpectraSTConstants.hpp"
#include "Peptide.hpp"
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

/* Class: SpectraSTTxtSearchOutput
 * 
 * Outputter to .txt format
 *    
 */


// constructor
SpectraSTTxtSearchOutput::SpectraSTTxtSearchOutput(string outFullFileName, string inputFileExt,
                                                   SpectraSTSearchParams &searchParams, bool printHeader) :
        SpectraSTSearchOutput(outFullFileName, inputFileExt, searchParams),
        m_printHeader(printHeader) {

}

// destructor
SpectraSTTxtSearchOutput::~SpectraSTTxtSearchOutput() {
}

void SpectraSTTxtSearchOutput::printHeader() {
    if (!m_fout) openFile();

    if (m_printHeader) {
        (*m_fout).width(MAX_NAME_LEN);
        (*m_fout) << left << "### Query";
        (*m_fout).width(4);
        (*m_fout) << left << "Rk";
        (*m_fout).width(MAX_NAME_LEN);
        (*m_fout) << left << "ID";

        SpectraSTSimScores::printHeaderFixedWidth(*m_fout);

        (*m_fout).width(10);
        (*m_fout) << left << "Status";
        (*m_fout).width(10);
        (*m_fout) << left << "Inst";
        (*m_fout).width(6);
        (*m_fout) << left << "Spec";
        (*m_fout).width(20);
        (*m_fout) << left << "LibFileOffset";
        (*m_fout).width(MAX_NAME_LEN);
        (*m_fout) << left << "Proteins";

        (*m_fout) << endl;
    }

}


// printStartQuery - prints the query information
void
SpectraSTTxtSearchOutput::printStartQuery(string query, double precursorMz, int assumedCharge, double retentionTime) {

    (*m_fout).width(MAX_NAME_LEN);
    (*m_fout) << left << query;

    m_query = query;
}

// printStartQuery - prints the hit
void SpectraSTTxtSearchOutput::printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry,
                                        SpectraSTSimScores &simScores) {

    if (!entry) {

        (*m_fout).width(4);
        (*m_fout) << "0";
        (*m_fout).width(MAX_NAME_LEN);
        (*m_fout) << left << "NO_MATCH";
        (*m_fout) << endl;

        return;
    }

    if (hitRank > 1) {
        (*m_fout).width(MAX_NAME_LEN);
        (*m_fout) << left << "--";
    }

    (*m_fout).width(4);
    (*m_fout) << left << hitRank;

    stringstream peptidess;
    peptidess << entry->getName();
    (*m_fout).width(MAX_NAME_LEN);
    (*m_fout) << left << peptidess.str();
    simScores.printFixedWidth((*m_fout));

    (*m_fout).width(10);
    (*m_fout) << left << entry->getStatus();

    (*m_fout).width(10);
    string inst("");
    if (entry->getOneComment("Inst", inst)) {
        (*m_fout) << left << inst.substr(0, (inst.length() < 8 ? inst.length() : 8));
    } else {
        (*m_fout) << left << "Unk";
    }

    (*m_fout).width(6);
    string spectrumType("");
    if (entry->getOneComment("Spec", spectrumType)) {
        (*m_fout) << left << spectrumType.substr(0, (spectrumType.length() < 3 ? spectrumType.length() : 3));
    } else {
        (*m_fout) << left << "Unk";
    }

    (*m_fout).width(20);
    (*m_fout) << left << entry->getLibFileOffset();

    string protein("");
    stringstream proteinss;
    string sample("");
    string::size_type pos = 0;
    if (entry->getOneComment("Protein", protein)) {
        if (protein.find('/') != string::npos) {
            int proteinCount = atoi(nextToken(protein, 0, pos, "/\t\r\n").c_str());
            proteinss << proteinCount << ':';
            bool isFirstProtein = true;
            while (pos < protein.length()) {
                if (!isFirstProtein) proteinss << ';';
                proteinss << nextToken(protein, pos, pos, "/\t\r\n", "/");
                isFirstProtein = false;
            }
        } else {
            proteinss << '1' << ':' << protein;
        }

    } else if (entry->getOneComment("Sample", sample)) {
        // if no protein field is given, use the sample field. (this is useful for unidentified spectra)
        if (sample.find('/') != string::npos) {
            int sampleCount = atoi(nextToken(sample, 0, pos, "/\t\r\n").c_str());
            proteinss << sampleCount << ':';
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

    (*m_fout) << proteinss.str();

    (*m_fout) << endl;

}

// printAbortedQuery - displays a message for a query that is not searched (used for dta input, for example)
void SpectraSTTxtSearchOutput::printAbortedQuery(string query, string message) {

    (*m_fout).width(MAX_NAME_LEN);
    (*m_fout) << left << query;
    (*m_fout).width(MAX_NAME_LEN);
    (*m_fout) << left << message;
    (*m_fout) << endl;
}
