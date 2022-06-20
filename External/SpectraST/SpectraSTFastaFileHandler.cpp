#pragma warning (disable: 4503)

#include "SpectraSTFastaFileHandler.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"

#include <iostream>
#include <sstream>
#include <stdlib.h>

#ifndef __LGPL__

#ifdef STANDALONE_LINUX

#include "SpectraST_kwset.h"

#else
#include "refresh_interact/kwset.h"
#endif

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.08 


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

#else // #ifndef __LGPL__

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.08 


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

#endif // #ifndef __LGPL__


/* Class: SpectraSTFastaFileHandler
 * 
 * Manages a .fasta protein sequence file.
 * 
 * MAJOR WORK IN PROGRESS
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTFastaFileHandler::SpectraSTFastaFileHandler(string fastaFileName) :
        m_fin(),
        m_offsets(NULL),
        m_fastaFileName(fastaFileName),
        m_nextProteinOffset(0) {

    if (!myFileOpen(m_fin, m_fastaFileName, true)) {
        g_log->error("CREATE", "Cannot open .fasta file \"" + m_fastaFileName + "\". Exiting.");
        g_log->crash();
    }

}

// destructor
SpectraSTFastaFileHandler::~SpectraSTFastaFileHandler() {
    m_fin.close();
    if (m_offsets) {
        delete (m_offsets);
    }

}

// enableRandomAccess - goes through the entire .fasta file, and remembers the file offsets of all
// protein headers. Stores the offsets with the protein IDs, so that the protein can be retrieved by
// random access if the ID is specified.
void SpectraSTFastaFileHandler::enableRandomAccess() {

    if (!m_offsets) {
        m_offsets = new map<string, fstream::off_type>;
    } else {
        return;
    }

    // reset stream to beginning of file
    if (m_fin.eof()) {
        m_fin.close();
        if (!myFileOpen(m_fin, m_fastaFileName, true)) {
            g_log->error("CREATE",
                         "Error opening the .fasta file \"" + m_fastaFileName + "\" for enabling random access.");
            g_log->crash();
            return;
        }
    } else {
        m_fin.seekg(0);
    }

    string line("");

    do {

        fstream::off_type offset = 0;
        while (nextLineCropped(m_fin, line, ">", "")) {
            offset = m_fin.tellg();
        }

        if (line != "_EOF_") {
            // reach a line beginning with >
            string::size_type spacePos = 0;
            string id = nextToken(line, 1, spacePos, " \t\r\n");

            // if a protein ID appears multiple times in a .fasta file, only
            // remembers the first occurrence

            if (m_offsets->find(id) == m_offsets->end()) {
                // remember offset -- this is before the last nextLineCropped is called, so
                // should point to the start of the line with the >
                (*m_offsets)[id] = offset;
            }
        }

    } while (line != "_EOF_");

}

// findProtein - use random-access to find the protein with a specified ID. Returns true if found (and
// can successfully parse out desired protein), false if not found.
// The "ID" is the string immediately following the '>' in the header, until the first whitespace.
bool SpectraSTFastaFileHandler::findProtein(string &id, string &fullDescr, string &seq) {

    if (!m_offsets) {
        enableRandomAccess();
    }

    if (m_fin.eof()) {
        m_fin.close();
        if (!myFileOpen(m_fin, m_fastaFileName, true)) {
            g_log->error("CREATE",
                         "Error opening the .fasta file \"" + m_fastaFileName + "\" for enabling random access.");
            g_log->crash();
            return (false);
        }
    }

    map<string, fstream::off_type>::iterator found = m_offsets->find(id);
    if (found != m_offsets->end()) {
        m_fin.seekg(found->second);
        return (parseNextProtein(id, fullDescr, seq));
    } else {
        fullDescr = "";
        seq = "";
        return (false);
    }

}

bool SpectraSTFastaFileHandler::nextProtein(string &id, string &fullDescr, string &seq) {

    if (m_nextProteinOffset != -1 && m_fin.eof()) {
        m_fin.close();
        if (!myFileOpen(m_fin, m_fastaFileName, true)) {
            g_log->error("CREATE", "Error opening the .fasta file \"" + m_fastaFileName + "\".");
            g_log->crash();
            return (false);
        }
    }

    m_fin.seekg(m_nextProteinOffset);
    bool success = parseNextProtein(id, fullDescr, seq);
    m_nextProteinOffset = m_fin.tellg();
    return (success);

}

// parseNextProtein - helper method to parse the next protein (starting from the current location of the input stream). Will
// first find the nearest header, then concatenate all sequence lines, then stop at the next '>' or EOF.
// Returns true if successful, false if no protein can be found downstream from the current location.
bool SpectraSTFastaFileHandler::parseNextProtein(string &id, string &fullDescr, string &seq) {

    string line("");

    // go to the first '>' line from this point on
    while (m_fin.peek() != '>') {
        if (!nextLine(m_fin, line)) {
            // end of file before finding any '>' -- no more protein!
            id = "";
            fullDescr = "";
            seq = "";
            return (false);
        }
    }

    // parse the '>' line
    bool dummy = nextLine(m_fin, line);  // should never be false
    line = crop(line);

    string::size_type pos = 0;
    id = nextToken(line, 1, pos, " \t\r\n");  // remove the '>'
    fullDescr = line.substr(1);  // remove the '>'
    seq = "";

    // parse the sequence until we hit another '>'
    while (m_fin.peek() != '>') {
        if (!nextLine(m_fin, line)) {
            return (true);
        }
        seq += crop(line);
    }
    return (true);
}

void
SpectraSTFastaFileHandler::refresh(map<string, vector<pair<string, string> > *> &pepMappings, bool refreshTrypticOnly) {

#ifdef __LGPL__

    g_log->error("CREATE", "LGPL version of SpectraST has no REFRESH function. REFRESH not performed.");
    return;

#else

    kwset_t kws = kwsalloc(NULL);
    if (!kws) {
        g_log->error("CREATE",
                     "Error allocating memory for keyword search trie for peptide-protein remapping. Not performed.");
        return;
    }

    vector<map<string, vector<pair<string, string> > *>::iterator> peps;

    for (map<string, vector<pair<string, string> > *>::iterator i = pepMappings.begin(); i != pepMappings.end(); i++) {
        if (kwsincr(kws, i->first.c_str(), i->first.length())) {
            g_log->error("CREATE", "Error adding peptide sequence \"" + i->first +
                                   "\" to keyword search trie for remapping. Not remapped.");
        } else {
            peps.push_back(i);
        }
    }

    if (kwsprep(kws)) {
        g_log->error("CREATE",
                     "Error preparing keyword search trie for peptide-protein remapping. No remapping performed.");
        kwsfree(kws);
        return;
    }

    // reset stream to beginning of file
    if (m_fin.eof()) {
        m_fin.close();
        if (!myFileOpen(m_fin, m_fastaFileName, true)) {
            g_log->error("CREATE",
                         "Error opening the .fasta file \"" + m_fastaFileName + "\" for peptide-protein remapping.");
            g_log->crash();
            return;
        }
    } else {
        m_fin.seekg(0);
    }

    string id("");
    string fullDescr("");
    string seq("");

    while (parseNextProtein(id, fullDescr, seq)) {

        struct kwsmatch *matches;

        int num_found = kwsexec_multiple(kws, seq.c_str(), seq.length(), &matches);

        if (num_found < 0) continue;

        for (int j = 0; j < num_found; j++) {
            int index = matches[j].index;

            int offset = (int) (matches[j].offset[0]);
            int length = (int) (matches[j].size[0]);

            stringstream contextss;
            for (int pre = offset - 4; pre < offset; pre++) {
                if (pre < 0 || pre >= (int) seq.length()) {
                    contextss << '-';
                } else if (pre == 0 && seq[pre] == 'M') {
                    contextss << '-';
                } else {
                    contextss << seq[pre];
                }
            }
            contextss << '_';
            for (int post = offset + length; post < offset + length + 4; post++) {
                if (post >= (int) seq.length() || post < 0) {
                    contextss << '-';
                } else {
                    contextss << seq[post];
                }
            }

            if (refreshTrypticOnly) {
                string context = contextss.str();
                Peptide testPep(peps[index]->first, 1);
                testPep.prevAA = context[3];
                testPep.nextAA = context[5];
                if (testPep.NTT() != 2) {
                    continue;
                }
            }

            if (!(peps[index]->second)) {
                peps[index]->second = new vector<pair<string, string> >;
            }

            pair<string, string> p;
            p.first = id;
            p.second = contextss.str();

            peps[index]->second->push_back(p);
        }

        free(matches);
    }

#endif // #ifndef __LGPL__

}

