#ifndef SPECTRASTLIB_HPP_
#define SPECTRASTLIB_HPP_

#include "SpectraSTLibEntry.hpp"
#include "SpectraSTMzLibIndex.hpp"
#include "SpectraSTPeptideLibIndex.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTCreateParams.hpp"
#include "FileUtils.hpp"
#include <string>
#include <vector>
#include <fstream>

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

/* Class: SpectraSTLib
 * 
 * The class that represents the library in memory. NOTE however that since the library
 * file .splib is huge, it is typically not read into memory in its entirety. Instead 
 * SpectraSTLib only keeps the index and the fstream objects in memory, and only goes 
 * to the .splib file to retrieve the entries when it is asked to, except for a few recently used entries cached 
 * in memory.
 * 
 */

using namespace std;


class SpectraSTLib {

public:
    SpectraSTLib(vector<string> &impFileNames, SpectraSTCreateParams *createParams);

    SpectraSTLib(string libFileName, SpectraSTSearchParams *searchParams, bool loadPeptideIndex = false);

    ~SpectraSTLib();

    // for create, import the library from file
    void insertEntry(SpectraSTLibEntry *entry);

    // for search, also used in sp lib importer, test is bad conflicting ID: isBadConflictingID()
    void retrieve(vector<SpectraSTLibEntry *> &hits, double lowMz, double highMz, bool shortAnnotation = true);

    SpectraSTPeptideLibIndex *getPeptideLibIndexPtr() { return (m_pepIndex); }


    // only used for importer
    void writePreamble(vector<string> &lines);



    // Fingerprint
    SpectraSTMzLibIndex *getMzLibIndexPtr() { return (m_mzIndex); }

    vector<vector<float> > &getFingerprint() { return (m_fingerprint); }
    // END Fingerprint




private:
    // never use functions -----------------
    void printDotTimeProfiles();
    // void printDotHistograms();

    string getLibFileName() { return (m_libFileName); }

    int getCount() { return (m_count); }

    void printStats();// not implemented
    // above functions are never used


    // m_mzIndex - points to the precursor m/z index object. This object is instantiated by SpectraSTLib and
    // IS the property of SpectraSTLib.
    SpectraSTMzLibIndex *m_mzIndex;

    // m_pepIndex - points to the peptide index object. This object is instantiated by SpectraSTLib and
    // IS the property of SpectraSTLib.
    SpectraSTPeptideLibIndex *m_pepIndex;

    // m_searchParams & m_createParams - points to the params object. The one corresponding to the
    // mode (Search/Create) will be instantiated by SpectraSTMain and passed into here (using the
    // respective constructor); the other will be set to NULL. These pointers are also used to keep track of the
    // which mode this SpectraSTLib object is created for.
    SpectraSTSearchParams *m_searchParams;
    SpectraSTCreateParams *m_createParams;

    // m_newLibId - counter used to assign a unique ID to each entry during Create mode
    unsigned int m_newLibId;

    // The file names. All these are the full file names that include the path and the extension
    string m_libFileName;   // the .splib file used for searching, or for holding the result of creating
    string m_txtFileName;   // if the .splib is binary, a separate text file .sptxt is also printed for human readability
    vector<string> m_impFileNames; // files to be imported and create library from
    string m_mzIdxFileName;  // the corresponding .spidx file
    string m_pepIdxFileName; // the corresponding .pepidx file
    FileName m_libFileNameStruct;

    // fstream objects for i/o
    ifstream m_libFin;
    ofstream m_libFout;
    ofstream m_txtFout;
    ofstream *m_mrmFout;
    ofstream *m_mgfFout;
    ofstream *m_PAIdentFout;

    // Counter
    unsigned int m_count;

    // Flag to indicate we cannot open a .sptxt file to write
    bool m_noSptxt;

    // Utility functions for initialization
    void initializeLibSearchMode(bool loadPeptideIndex);

    void initializeLibCreateMode();

    void extractDatabaseFileFromPreamble(bool binary);

    // Fingerprint
    vector<vector<float> > m_fingerprint;

    static bool sortFingerprintDsc(float a, float b);
    void calcLibFingerprint();
    // END Fingerprint


};

#endif /*SPECTRALIB_HPP_*/
