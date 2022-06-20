#include "SpectraSTPeptideLibIndex.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"

#include <sstream>
#include <stdlib.h>

#define NORMALCELLCOLOR   "#FFDDDD"
#define HEADERCELLCOLOR   "#42D4FD"

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

/* Class: SpectraSTPeptideLibIndex
 * 
 * Implements a library index on the peptide to facilitate retrieval by peptide. 
 * Note that this is not used by SpectraST for search, but it is useful for library manipulation
 */

using namespace std;

extern SpectraSTLog *g_log;

// constructor for creation
SpectraSTPeptideLibIndex::SpectraSTPeptideLibIndex(string idxFileName) :
        SpectraSTLibIndex(idxFileName, "Peptide"),
        m_map(),
        m_peptideSequenceCount(0),
        m_peptideIonCount(0),
        m_peptideSpectrumCount(0),
        m_chargeCount(7, 0),
        m_trypticCount(0),
        m_semitrypticCount(0),
        m_nontrypticCount(0),
        m_modsCount(),
        m_nonPeptideCount(0),
        m_nonPeptideIonCount(0),
        m_nonPeptideSpectrumCount(0),
        m_isUnique(true),
        m_prob9999(0),
        m_prob999(0),
        m_prob99(0),
        m_prob9(0),
        m_prob0(0),
        m_nreps20(0),
        m_nreps10(0),
        m_nreps4(0),
        m_nreps2(0),
        m_nreps1(0) {

    m_curPeptide = m_map.end(); // forbid iteration if it's created for inserting

}

// constructor for retrieval
SpectraSTPeptideLibIndex::SpectraSTPeptideLibIndex(string idxFileName, ifstream *libFinPtr, bool binaryLib) :
        SpectraSTLibIndex(idxFileName, libFinPtr, "Peptide", binaryLib),
        m_map(),
        m_peptideSequenceCount(0),
        m_peptideIonCount(0),
        m_peptideSpectrumCount(0),
        m_chargeCount(7, 0),
        m_trypticCount(0),
        m_semitrypticCount(0),
        m_nontrypticCount(0),
        m_modsCount(),
        m_nonPeptideCount(0),
        m_nonPeptideIonCount(0),
        m_nonPeptideSpectrumCount(0),
        m_isUnique(true),
        m_prob9999(0),
        m_prob999(0),
        m_prob99(0),
        m_prob9(0),
        m_prob0(0),
        m_nreps20(0),
        m_nreps10(0),
        m_nreps4(0),
        m_nreps2(0),
        m_nreps1(0) {


    // read the index into memory from the file
    readFromFile();
    m_curPeptide = m_map.begin();

}

// destructor
SpectraSTPeptideLibIndex::~SpectraSTPeptideLibIndex() {

}

// insertEntry - used to insert an entry to the index
void SpectraSTPeptideLibIndex::insertEntry(SpectraSTLibEntry *entry, fstream::off_type offset) {

    double prob = entry->getProb();

    if (prob > 0.9999) {
        m_prob9999++;
    } else if (prob > 0.999) {
        m_prob999++;
    } else if (prob > 0.99) {
        m_prob99++;
    } else if (prob > 0.9) {
        m_prob9++;
    } else {
        m_prob0++;
    }

    unsigned int nreps = entry->getNrepsUsed();

    if (nreps >= 20) {
        m_nreps20++;
    } else if (nreps >= 10) {
        m_nreps10++;
    } else if (nreps >= 4) {
        m_nreps4++;
    } else if (nreps >= 2) {
        m_nreps2++;
    } else {
        m_nreps1++;
    }

    m_entryCount++;

    Peptide *pep = entry->getPeptidePtr();

    if (!pep) {
        m_nonPeptideSpectrumCount++;

        string::size_type slashPos = 0;
        int charge = entry->getCharge();
        string name = entry->getName();
        if ((slashPos = name.rfind('/')) != string::npos) {
            name = name.substr(0, slashPos);
        }
        stringstream subkeyss;
        subkeyss << charge << '|' << '0' << '|';
        string frag(entry->getFragType());
        if (!(frag.empty())) {
            subkeyss << frag;
        }

        if (m_map.find(name) == m_map.end()) {
            m_nonPeptideCount++;
            m_nonPeptideIonCount++;
        } else if ((m_map[name]).find(subkeyss.str()) == (m_map[name]).end()) {
            m_nonPeptideIonCount++;
        } else {
            m_isUnique = false;
        }

        if (charge > 5) {
            m_chargeCount[6]++;
        } else if (charge >= 0) {
            m_chargeCount[charge]++;
        }

        ((m_map[name])[subkeyss.str()]).push_back(offset);

        return;
    }

    // peptides

    m_peptideSpectrumCount++;

    // create key for sub-map. The sub-key is of the format <charge>|<mod string in .msp style>
    string subkey = constructSubkey(entry);

    map<string, map<string, vector<fstream::off_type> > >::iterator found = m_map.find(pep->stripped);
    if (found == m_map.end()) {
        // new stripped peptide
        m_peptideSequenceCount++;
        m_peptideIonCount++;
    } else if ((m_map[pep->stripped]).find(subkey) == (m_map[pep->stripped]).end()) {
        // new peptide ion
        m_peptideIonCount++;
    } else {
        m_isUnique = false;
    }

    if (pep->charge > 5) {
        m_chargeCount[6]++;
    } else if (pep->charge >= 0) {
        m_chargeCount[pep->charge]++;
    }

    map<string, unsigned int> presentModTypes;
    if (pep->getAllPresentModTypes(presentModTypes)) {
        for (map<string, unsigned int>::iterator m = presentModTypes.begin(); m != presentModTypes.end(); m++) {
            string tag = m->first;
            if (m_modsCount.find(tag) != m_modsCount.end()) {
                m_modsCount[tag]++;
            } else {
                m_modsCount[tag] = 1;
            }
        }
    }

    if (pep->NTT() == 2) {
        m_trypticCount++;
    } else if (pep->NTT() == 1) {
        m_semitrypticCount++;
    } else {
        m_nontrypticCount++;
    }

    // insert
    ((m_map[pep->stripped])[subkey]).push_back(offset);


}

// writeToFile - writes the peptide index to file. NOTE that this function must be modified together with
// readFromFile if the file format is changed!
void SpectraSTPeptideLibIndex::writeToFile() {

    ofstream idxFout;

    if (!myFileOpen(idxFout, m_idxFileName)) {
        g_log->error("CREATE", "Cannot open PEPIDX file \"" + m_idxFileName + " for writing peptide index.");
        return;
    }

    idxFout << "### " << m_idxFileName << endl;
    idxFout << "### " << endl;

    printStats(idxFout, "### ");

    idxFout << "### ===" << endl;

    for (map<string, map<string, vector<fstream::off_type> > >::iterator i = m_map.begin(); i != m_map.end(); i++) {

        for (map<string, vector<fstream::off_type> >::iterator j = ((*i).second).begin();
             j != ((*i).second).end(); j++) {

            idxFout << (*i).first << '\t' << (*j).first << '\t';

            for (vector<fstream::off_type>::iterator k = ((*j).second).begin(); k != ((*j).second).end(); k++) {
                idxFout << (*k) << ' ';
            }

            idxFout << endl;
        }
    }


}

// readFromFile - reads the peptide index from file and creates the map object in memory. 
// NOTE that this function must be modified together with writeToFile if the file format is changed!
void SpectraSTPeptideLibIndex::readFromFile() {

    ifstream idxFin;
    if (!myFileOpen(idxFin, m_idxFileName)) {
        g_log->error("CREATE", "Cannot open PEPIDX file \"" + m_idxFileName + " for reading peptide index. Exiting.");
        g_log->crash();
    }

    string line("");

    char firstChar = (char) (idxFin.peek());
    if (firstChar == '#') {
        while (nextLine(idxFin, line, "### ==="));
    }

    string lastKey("");

    while (nextLine(idxFin, line, "", "")) {
        string::size_type pos = 0;

        string key = nextToken(line, pos, pos, "\t\r\n");
        string subkey = nextToken(line, pos, pos, "\t\r\n");
        int charge = 0;
        string mods("");
        string frag("");
        parseSubkey(subkey, charge, mods, frag);

        Peptide *pep = NULL;
        if (key[0] == '_') {
            // not a peptide
            m_nonPeptideIonCount++;

            if (key != lastKey) {
                m_nonPeptideCount++;
            }

        } else {

            pep = new Peptide(key, charge, mods);

            //    if (!pep->isGood()) {
            //      g_log->error("GENERAL", "Unrecognized modification: " + mods + " ? Please check that your usermods file is in the current directory.");
            //      g_log->crash();
            //    }

            m_peptideIonCount++;

            if (key != lastKey) {
                m_peptideSequenceCount++;
            }
        }

        unsigned int spectrumCount = 0;

        while (pos < line.length()) {
            fstream::off_type offset = strtoull((nextToken(line, pos, pos)).c_str(), NULL, 10);
            ((m_map[key])[subkey]).push_back(offset);

            spectrumCount++;

            pos++;
        }

        m_entryCount += spectrumCount;

        if (charge > 5) {
            m_chargeCount[6] += spectrumCount;
        } else if (charge >= 0) {
            m_chargeCount[charge] += spectrumCount;
        }

        if (key[0] == '_') {

            m_nonPeptideSpectrumCount += spectrumCount;

        } else {
            // is a peptide

            m_peptideSpectrumCount += spectrumCount;

            map<string, unsigned int> presentModTypes;
            if (pep->getAllPresentModTypes(presentModTypes)) {
                for (map<string, unsigned int>::iterator m = presentModTypes.begin(); m != presentModTypes.end(); m++) {
                    string tag = m->first;
                    if (m_modsCount.find(tag) != m_modsCount.end()) {
                        m_modsCount[tag] += spectrumCount;
                    } else {
                        m_modsCount[tag] = spectrumCount;
                    }
                }
            }

            if (pep->NTT() == 2) {
                m_trypticCount++;
            } else if (pep->NTT() == 1) {
                m_semitrypticCount++;
            } else {
                m_nontrypticCount++;
            }

            delete (pep);

        } // if (pep)

        if (((m_map[key])[subkey]).size() > 1) {
            m_isUnique = false;
        }

        lastKey = key;

    } // next line in .pepidx file

}

// retrieve - gets all the hits with a particular stripped peptide, charge or mods. If charge and/or mods is
// not given, all the entries with the same stripped peptide will be retrieved.
void SpectraSTPeptideLibIndex::retrieve(vector<SpectraSTLibEntry *> &hits, string peptide, int charge, string mods,
                                        string frag) {

    // find the stripped peptide
    map<string, map<string, vector<fstream::off_type> > >::iterator foundIndex = m_map.find(peptide);
    if ((foundIndex == m_map.end())) {
        // not found
        return;
    }

    // loops over all the entries with that stripped peptide, and determines which ones we want based on the
    // charge and the mods
    for (map<string, vector<fstream::off_type> >::iterator i = (*foundIndex).second.begin();
         i != (*foundIndex).second.end(); i++) {

        string subkey((*i).first);
        int subkeyCharge = 0;
        string subkeyMods("");
        string subkeyFrag("");
        parseSubkey(subkey, subkeyCharge, subkeyMods, subkeyFrag);

        if (/*frag.empty() ||*/ (frag == subkeyFrag)) { // disallow "any frag" retrieval
            // don't care about frag, or frag matched
            if ((charge == 0) || (charge == subkeyCharge)) {
                // don't care about charge, or charge matched
                if (mods.empty() || (mods == subkeyMods)) {
                    // don't care about mods, or mods matched
                    for (vector<fstream::off_type>::iterator j = (*i).second.begin(); j != (*i).second.end(); j++) {
                        m_libFinPtr->seekg(*j);
                        SpectraSTLibEntry *entry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
                        //					entry->readFromFile(*m_libFinPtr);

                        hits.push_back(entry);
                    }
                }
            }
        }
    }
}

void SpectraSTPeptideLibIndex::retrieve(vector<SpectraSTLibEntry *> &hits, string peptide, string subkey) {

    int charge = 0;
    string mods("");
    string frag("");
    parseSubkey(subkey, charge, mods, frag);

    retrieve(hits, peptide, charge, mods, frag);
}


// isInIndex - similar to retrieve(), but only returns whether or not an entry is found, and will not actually retrieve the entry. 
bool SpectraSTPeptideLibIndex::isInIndex(string peptide, int charge, string mods, string frag) {

    map<string, map<string, vector<fstream::off_type> > >::iterator foundIndex = m_map.find(peptide);
    if ((foundIndex == m_map.end())) {
        // not found
        return (false);
    }

    // loops over all the entries with that stripped peptide, and determines which ones we want based on the
    // charge and the mods
    for (map<string, vector<fstream::off_type> >::iterator i = (*foundIndex).second.begin();
         i != (*foundIndex).second.end(); i++) {

        string subkey((*i).first);
        int subkeyCharge = 0;
        string subkeyMods("");
        string subkeyFrag("");
        parseSubkey(subkey, subkeyCharge, subkeyMods, subkeyFrag);

        if (/*frag.empty() ||*/ (frag == subkeyFrag)) { // disallow "any frag" retrieval
            // don't care about frag, or frag matched
            if ((charge == 0) || (charge == subkeyCharge)) {
                // don't care about charge, or charge matched
                if (mods.empty() || (mods == subkeyMods)) {
                    // don't care about mods, or mods matched
                    return (true);
                }
            }
        }
    }
    return (false);
}

bool SpectraSTPeptideLibIndex::isInIndex(string peptide, string subkey) {

    int charge = 0;
    string mods("");
    string frag("");
    parseSubkey(subkey, charge, mods, frag);

    return (isInIndex(peptide, charge, mods, frag));
}

// nextPeptide - sequential access method. See SpectraSTSpLibImporter::import() for an example of how to use it
bool SpectraSTPeptideLibIndex::nextPeptide(string &peptide, vector<string> &subkeys) {

    peptide = "";
    subkeys.clear();

    if (m_curPeptide == m_map.end()) {
        return (false);
    }

    peptide = (*m_curPeptide).first;

    for (map<string, vector<fstream::off_type> >::iterator i = (*m_curPeptide).second.begin();
         i != (*m_curPeptide).second.end(); i++) {
        string subkey = (*i).first;

        subkeys.push_back(subkey);

    }

    m_curPeptide++;
    return (true);

}

void SpectraSTPeptideLibIndex::reset() {
    // cerr << " <<< resetting >>> " << endl;
    m_curPeptide = m_map.begin();
}

// hasSharedSequence - finds a peptide in the index with a shared sequence as pep.
bool SpectraSTPeptideLibIndex::hasSharedSequence(Peptide &pep, string &foundPeptide, string frag) {

    // first, check to see if there's another peptide ion with the same stripped peptide
    map<string, map<string, vector<fstream::off_type> > >::iterator found = m_map.find(pep.stripped);
    if (found != m_map.end()) {
        if (found->second.size() > 1) {
            // must have another peptide ion with same stripped peptide as this one
            foundPeptide = pep.stripped;
            return (true);
        } else {
            // check to make sure this only one peptide ion is not self
            stringstream subkeyss;
            subkeyss << pep.charge << '|' << pep.mspMods() << '|' << frag;
            if (subkeyss.str() != found->second.begin()->first) {
                foundPeptide = pep.stripped;
                return (true);
            }
        }
    }

    // next, go through the entire index to find any stripped peptide with shared sequence - quite inefficient
    for (map<string, map<string, vector<fstream::off_type> > >::iterator i = m_map.begin(); i != m_map.end(); i++) {
        if ((pep.stripped != i->first) &&
            ((pep.stripped.length() > i->first.length() && pep.stripped.find(i->first, 0) != string::npos) ||
             (pep.stripped.length() < i->first.length() && i->first.find(pep.stripped, 0) != string::npos))) {
            foundPeptide = i->first;
            return (true);
        }
    }
    return (false);
}

// isUniqueLibrary - returns TRUE if the library is "unique" and FALSE otherwise
bool SpectraSTPeptideLibIndex::isUniqueLibrary() {
    return (m_isUnique);
}

// printStats - prints out statistics
void SpectraSTPeptideLibIndex::printStats(ostream &out, string linePrefix) {

    out << linePrefix << "Total number of spectra in library: " << m_entryCount << endl;

    if (m_entryCount != m_peptideSpectrumCount) {
        out << linePrefix << "Total number of peptide spectra in library: " << m_peptideSpectrumCount << endl;
    }

    out << linePrefix << "Total number of distinct peptide ions in library: " << m_peptideIonCount << endl;
    out << linePrefix << "Total number of distinct stripped peptides in library: " << m_peptideSequenceCount << endl;

    if (m_nonPeptideSpectrumCount > 0) {
        out << linePrefix << "Total number of non-peptide spectra in library: " << m_nonPeptideSpectrumCount << endl;
        out << linePrefix << "Total number of distinct non-peptide ions in library: " << m_nonPeptideIonCount << endl;
        out << linePrefix << "Total number of distinct non-peptide entities in library: " << m_nonPeptideCount << endl;
    }

    out << linePrefix << endl;

    out << linePrefix << "CHARGE            ";
    for (unsigned int c = 1; c <= 5; c++) {
        out << "+" << c << ": " << m_chargeCount[c] << " ; ";
    }
    out << ">+5: " << m_chargeCount[6] << " ; ";
    out << "Unk: " << m_chargeCount[0] << endl;

    out << linePrefix << "TERMINI           Tryptic: " << m_trypticCount << " ; Semi-tryptic: " << m_semitrypticCount
        << " ; Non-tryptic: " << m_nontrypticCount << endl;
    out << linePrefix << "PROBABILITY       >0.9999: " << m_prob9999 << " ; 0.999-0.9999: " << m_prob999
        << " ; 0.99-0.999: " << m_prob99 << " ; 0.9-0.99: " << m_prob9 << " ; <0.9: " << m_prob0 << endl;
    out << linePrefix << "NREPS             20+: " << m_nreps20 << " ; 10-19: " << m_nreps10 << " ; 4-9: " << m_nreps4
        << " ; 2-3: " << m_nreps2 << " ; 1: " << m_nreps1 << endl;

    out << linePrefix << "MODIFICATIONS     ";
    char lastAA = '0';
    for (map<string, unsigned int>::iterator m = m_modsCount.begin(); m != m_modsCount.end(); m++) {
        if (m->first[0] == lastAA) {
            out << " ; ";
        } else if (lastAA != '0') {
            out << endl << linePrefix << "                  ";
        }
        out << m->first << ": " << m->second;
        lastAA = m->first[0];
    }

    if (lastAA == '0') {
        // no modification at all
        out << "None";
    }

    out << endl;

}

void SpectraSTPeptideLibIndex::getAllSequences(vector<string> &seqs) {

    seqs.clear();
    for (map<string, map<string, vector<fstream::off_type> > >::iterator i = m_map.begin(); i != m_map.end(); i++) {
        seqs.push_back(i->first);
    }
}

void SpectraSTPeptideLibIndex::parseSubkey(string &subkey, int &charge, string &mods, string &frag) {

    string::size_type dividerPos = subkey.find('|', 0);
    charge = atoi((subkey.substr(0, dividerPos)).c_str());

    string::size_type secondDividerPos = subkey.find('|', dividerPos + 1);

    mods = "";
    frag = "";
    if (secondDividerPos == string::npos) {
        mods = subkey.substr(dividerPos + 1);
    } else {
        // new style, with fragmentation type indicator
        mods = subkey.substr(dividerPos + 1, secondDividerPos - dividerPos - 1);
        if (secondDividerPos + 1 < subkey.length()) {
            frag = subkey.substr(secondDividerPos + 1);
        }
    }

}

string SpectraSTPeptideLibIndex::constructSubkey(SpectraSTLibEntry *entry) {

    stringstream subkey;
    Peptide *pep = entry->getPeptidePtr();
    subkey << pep->charge << '|' << pep->mspMods() << '|';

    string frag(entry->getFragType());
    if (!(frag.empty())) {
        subkey << frag;
    }

    return (subkey.str());

}
