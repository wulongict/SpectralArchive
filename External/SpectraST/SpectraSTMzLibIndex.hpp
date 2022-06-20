#ifndef SPECTRASTMZLIBINDEX_HPP_
#define SPECTRASTMZLIBINDEX_HPP_

#include "SpectraSTLibEntry.hpp"
#include "SpectraSTLibIndex.hpp"
#include <iostream>
#include <vector>
#include <queue>
#include <string>


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

/* Class: SpectraSTMzLibIndex
 * 
 * Implements a library index on the precursor m/z value. 
 * 
 * NOTE ON CACHING: The recently retrieved entries are cached in memory for efficiency. The cache capacity
 * is specified by the cacheRange argument in the constructor for retrieval. See the retrieve() method for
 * more information.
 * 
 */

using namespace std;

class SpectraSTMzLibIndex : public SpectraSTLibIndex {


public:

    // constructor for creation
    SpectraSTMzLibIndex(string idxFileName);

    // constructor for retrieval
    SpectraSTMzLibIndex(string idxFileName, ifstream *libFinPtr, double cacheRange, bool binaryLib);

    virtual ~SpectraSTMzLibIndex();

    // Insertion method - creation mode only
    virtual void insertEntry(SpectraSTLibEntry *entry, fstream::off_type offset);

    // File I/O methods
    virtual void writeToFile();

    virtual void readFromFile();

    // Retrieval method
    void retrieve(vector<SpectraSTLibEntry *> &hits, double lowMz, double highMz, bool shortAnnotation = false);

    // Sequential access method - the returned SpectraSTLibEntry object becomes property of caller!
    SpectraSTLibEntry *nextEntry();

    bool nextFileOffset(fstream::off_type &offset);

    SpectraSTLibEntry *thisEntry();

    void reset();

    // Sequential access methods with sorting - the returned SpectraSTLibEntry object becomes property of caller!
    void sortEntriesByNreps();

    void sortEntriesBySN();

    // added by LWU
    void sortEntriesByOffsets();

    SpectraSTLibEntry *nextSortedEntry();

    bool nextSortedFileOffset(fstream::off_type &offset);

    SpectraSTLibEntry *thisSortedEntry();

    void print();

    void printDotTimeProfiles(ofstream &fout);
    // void printDotHistograms(ofstream& fout);

    static const double CACHE_ALL;


private:


    // m_offsets - The file offsets into the .splib file. Each bin is 1 m/z unit wide, and contains a list of file offsets that
    // points to all entries of that m/z value in the .splib files.
    vector<vector<fstream::off_type> > m_offsets;

    // m_cache - The cache, pointers to 1 m/z-unit-wide bins of SpectraSTLibEntry's.
    // The cache is allocated for the entire m/z range, but at any given time, only a fraction
    // of the bins are active. When a bin is active (in use), the pointer points to a vector of
    // SpectraSTLibEntry's. When it is inactive, the pointer is NULL.
    vector<vector<SpectraSTLibEntry *> *> m_cache;

    // m_cacheQueue - queue to remember which bin is created for the longest, which will be
    // the first to be inactivated when the cache goes over capacity.
    queue<unsigned int> m_cacheQueue;

    // m_cacheSize - the current number of active bins
    int m_cacheSize;

    // m_cacheCapacity - the maximum allowable number of active bins. If m_cacheSize gets bigger
    // than m_cacheCapacity, a bin of entries will be freed to maintain the cache size below capacity.
    int m_cacheCapacity;

    // m_curBin, m_curOffset - indices into the m_offsets 2-D array. For sequential access using nextEntry()
    unsigned int m_curBin;
    int m_curOffset;

    // m_sortedOffsets - vectors of sorted (by whatever criterion, Nreps is the only implemented so far) offsets
    vector<pair<fstream::off_type, double> > *m_sortedOffsets;

    // m_curSortedOffset - index into the m_sortedOffets vector. For sequential access using nextSortedEntry()
    int m_curSortedOffset;

    // Private methods
    void initialize();

    void freeCacheBin(unsigned int bin);

    unsigned int calcBinNumber(double mass);

    unsigned int calcBinMz(unsigned int binNum);

    static bool sortEntriesDesc(pair<fstream::off_type, double> a, pair<fstream::off_type, double> b);

};


#endif /*SPECTRAMZLIBINDEX_H_*/
