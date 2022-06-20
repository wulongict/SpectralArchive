#ifndef SPECTRASTPEPTIDELIBINDEX_HPP_
#define SPECTRASTPEPTIDELIBINDEX_HPP_

#include "SpectraSTLibIndex.hpp"
#include <fstream>
#include <string>
#include <map>

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
 * Note that this is not used by SpectraST for search, but it is useful for library manipulation.
 */


using namespace std;

class SpectraSTPeptideLibIndex : public SpectraSTLibIndex {

public:
    SpectraSTPeptideLibIndex(string idxFileName);

    SpectraSTPeptideLibIndex(string idxFileName, ifstream *libFinPtr, bool binaryLib);

    virtual ~SpectraSTPeptideLibIndex();

    virtual void insertEntry(SpectraSTLibEntry *entry, fstream::off_type offset);

    // File I/O methods
    virtual void writeToFile();

    virtual void readFromFile();

    // Retrieval methods
    void
    retrieve(vector<SpectraSTLibEntry *> &hits, string peptide, int charge = 0, string mods = "", string frag = "");

    void retrieve(vector<SpectraSTLibEntry *> &hits, string peptide, string subkey);

    bool isInIndex(string peptide, int charge = 0, string mods = "", string frag = "");

    bool isInIndex(string peptide, string subkey);

    // static method for parsing subkey
    static void parseSubkey(string &subkey, int &charge, string &mods, string &frag);

    // Sequential access method
    bool nextPeptide(string &peptide, vector<string> &subkeys);

    void reset();

    // check if there is any peptide in the library has a shared sequence as pep
    bool hasSharedSequence(Peptide &pep, string &foundPeptide, string frag = "");

    // check if the library is a "unique" library, meaning each peptide ion is represented by only one spectrum
    bool isUniqueLibrary();

    // method to print the statistics
    void printStats(ostream &out, string linePrefix = "");

    // method to extract all distinct peptide sequences
    void getAllSequences(vector<string> &seqs);

    static string constructSubkey(SpectraSTLibEntry *entry);

private:

    // statistics

    unsigned int m_peptideSequenceCount;
    unsigned int m_peptideIonCount;
    unsigned int m_peptideSpectrumCount;

    vector<unsigned int> m_chargeCount;
    unsigned int m_trypticCount;
    unsigned int m_semitrypticCount;
    unsigned int m_nontrypticCount;

    map<string, unsigned int> m_modsCount;

    unsigned int m_nonPeptideCount;
    unsigned int m_nonPeptideIonCount;
    unsigned int m_nonPeptideSpectrumCount;

    unsigned int m_prob9999;
    unsigned int m_prob999;
    unsigned int m_prob99;
    unsigned int m_prob9;
    unsigned int m_prob0;

    unsigned int m_nreps20;
    unsigned int m_nreps10;
    unsigned int m_nreps4;
    unsigned int m_nreps2;
    unsigned int m_nreps1;

    bool m_isUnique;

    // the index
    map<string, map<string, vector<fstream::off_type> > > m_map;

    // an iterator to m_map
    map<string, map<string, vector<fstream::off_type> > >::iterator m_curPeptide;


};

#endif /*SPECTRASTPEPTIDELIBINDEX_HPP_*/
