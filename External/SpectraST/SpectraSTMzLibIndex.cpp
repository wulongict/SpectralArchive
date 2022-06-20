#include "SpectraSTMzLibIndex.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
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

#include <algorithm> // for std::sort

extern SpectraSTLog *g_log;

// Pass this as argument to cacheRange in the constructor to force SpectraST to keep everything
// in memory.
const double SpectraSTMzLibIndex::CACHE_ALL = 999999.0;

// constructor for creating
SpectraSTMzLibIndex::SpectraSTMzLibIndex(string idxFileName) :
        SpectraSTLibIndex(idxFileName, "PrecursorMZ"),
        m_offsets(MAX_MZ - MIN_MZ + 1),
        m_cacheSize(0),
        m_cacheCapacity(0),
        m_cache(MAX_MZ - MIN_MZ + 1),
        m_cacheQueue(),
        m_curBin(0),
        m_curOffset(-1),
        m_sortedOffsets(NULL),
        m_curSortedOffset(-1) {

    initialize();

}

// constructor for retrieval. Note that the cache capacity is set to the ceiling of cacheRange.
// e.g. if cacheRange is 6.5 Th, a maximum of 7 bins can be in memory at any given time.
SpectraSTMzLibIndex::SpectraSTMzLibIndex(string idxFileName, ifstream *libFinPtr, double cacheRange, bool binaryLib) :
        SpectraSTLibIndex(idxFileName, libFinPtr, "PrecursorMZ", binaryLib),
        m_offsets(MAX_MZ - MIN_MZ + 1),
        m_cacheSize(0),
        m_cacheCapacity((int) cacheRange + 1),
        m_cache(MAX_MZ - MIN_MZ + 1),
        m_cacheQueue(),
        m_curBin(0),
        m_curOffset(-1),
        m_sortedOffsets(NULL),
        m_curSortedOffset(-1) {

    initialize();
    readFromFile();

}

// destructor - Since all the entries stored in the cache remain property of the index, they
// have to be deleted here.
SpectraSTMzLibIndex::~SpectraSTMzLibIndex() {

    for (vector<vector<SpectraSTLibEntry *> *>::iterator cBin = m_cache.begin(); cBin != m_cache.end(); cBin++) {
        if (*cBin) {
            for (vector<SpectraSTLibEntry *>::iterator entry = (*cBin)->begin(); entry != (*cBin)->end(); entry++) {
                delete (*entry);
            }
            delete (*cBin);
        }
    }

    if (m_sortedOffsets) delete (m_sortedOffsets);

}

// insertEntry - given a new entry and its file offset, adds it into the index
void SpectraSTMzLibIndex::insertEntry(SpectraSTLibEntry *entry, fstream::off_type offset) {

    // finds out which bin this entry should go into
    unsigned int binNum = calcBinNumber((unsigned int) (entry->getPrecursorMz()));

    // insert the file offset
    (m_offsets[binNum]).push_back(offset);

    m_entryCount++;

}

// writeToFile - write the entire index to a file. The .spidx simply contains lines of the form
// <m/z value>\t<offset 1> <offset 2> ... <offset n>
void SpectraSTMzLibIndex::writeToFile() {

    ofstream idxFout;
    if (!myFileOpen(idxFout, m_idxFileName)) {
        g_log->error("CREATE", "Cannot open SPIDX file \"" + m_idxFileName + "\" for writing index.");
        return;
    }

    for (unsigned int b = 0; b < (unsigned int) (m_offsets.size()); b++) {
        idxFout << calcBinMz(b) << '\t';
        for (vector<fstream::off_type>::iterator i = m_offsets[b].begin(); i != m_offsets[b].end(); i++) {
            idxFout << (*i) << " ";
        }
        idxFout << endl;
    }

}

// readFromFile - read the entire index into memory from a file. MUST BE modified together with writeToFile()
// to ensure sychronization of the .spidx file format.
void SpectraSTMzLibIndex::readFromFile() {

    ifstream idxFin;
    if (!myFileOpen(idxFin, m_idxFileName)) {
        g_log->error("SEARCH", "Cannot open SPIDX file \"" + m_idxFileName + "\" for reading index. Exiting.");
        g_log->crash();
        return;
    }

    initialize();

    string line;

    while (nextLine(idxFin, line, "")) {

        string::size_type pos = 0;

        // first token in the line is the m/z value of that bin
        string massStr = nextToken(line, 0, pos);
        unsigned int binNum = calcBinNumber(atof(massStr.c_str()));

        // all tokens after the first are file offsets
        while (pos < line.length()) {
            string offsetStr = nextToken(line, pos, pos);
            if (offsetStr != "") {
                fstream::off_type offset = strtoull(offsetStr.c_str(), NULL, 10);
                (m_offsets[binNum]).push_back(offset);
                m_entryCount++;
            }
        }
    }

}

// retrieve - Given the target m/z and a tolerance around it, retrieves all library entries
// that fall into the corresponding bins, and return them in the vector 'entries'. It first checks
// if the required bins are already in memory (being active in the cache); if not, it updates the cache (and deletes old cached item
// if the cache becomes full) first. The returned entries remain property of SpectraSTMzLibIndex, and 
// WILL BE RE-USED for different queries!
//
// NOTE: The target m/z tolerance is not exact; it is possible that entries outside of the range 
// [targetMz - targetMzTolerance, targetMz + targetMzTolerance] are retrieved. This is necessarily the result of
// having an index that puts entries into 1 m/z unit-wide bins. Within the same bin, some entries may lie within the 
// range and others may not, and the index has no way of knowing which ones. 
// The function, however, is guaranteed to retrieve at least all the entries within the range.
void
SpectraSTMzLibIndex::retrieve(vector<SpectraSTLibEntry *> &entries, double lowMz, double highMz, bool shortAnnotation) {

    unsigned int low = calcBinNumber(lowMz);
    unsigned int high = calcBinNumber(highMz);

    for (unsigned int b = low; b <= high; b++) {
        // add by long
//        cout << "[Info] looking for bin: " << b << endl;
        if (m_cache[b] == NULL) {
            // this bin is not yet in memory, read it

            m_cache[b] = new vector<SpectraSTLibEntry *>();
            for (vector<fstream::off_type>::iterator offset = (m_offsets[b]).begin();
                 offset != (m_offsets[b]).end(); offset++) {
                // add by long
//                cout << "[Info] read disk: " << b << endl;
                m_libFinPtr->seekg(*offset);

                SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib, shortAnnotation);

                newEntry->setLibFileOffset(*offset);
                m_cache[b]->push_back(newEntry);
            }
            m_cacheSize++;

            // add by long
//            cout << "[Info] looking for bin: " << b << endl;

            // get into the cache queue
            m_cacheQueue.push(b);

            // if cache is full, pop the cache queue and free the oldest active bin
            if (m_cacheSize > m_cacheCapacity) {
                // add by long
                //cout << "[Info] cache is full: " << b << endl;

                unsigned int old = m_cacheQueue.front();
                while (old >= low && old <= high) { // todo: there is a bug here. If I use this to get a lot of  spectra, there will be a problem
                    // still need this bin, pop it but re-push to end of queue
                    m_cacheQueue.pop();
                    m_cacheQueue.push(old);
                    old = m_cacheQueue.front();
                }
                // add by long
                //cout << "[Info] found a useless bin: " << b << endl;

                // finally find an old bin that's not needed, pop and free this one
                m_cacheQueue.pop();
                freeCacheBin(old);
            }
        }
        //cout << "[Info] got bin " << b << endl;

        // now we can retrieve
        for (vector<SpectraSTLibEntry *>::iterator entry = m_cache[b]->begin(); entry != m_cache[b]->end(); entry++) {
            entries.push_back(*entry);
        }
        //cout << "[Info] end of bin " << b << endl;

    }


}

// freeCacheBin - frees a bin of entries in the cache
void SpectraSTMzLibIndex::freeCacheBin(unsigned int bin) {

    // check to make sure it's active
    if (m_cache[bin]) {

        for (vector<SpectraSTLibEntry *>::iterator entry = (m_cache[bin])->begin();
             entry != (m_cache[bin])->end(); entry++) {
            // delete the entries themselves
            delete (*entry);
        }
        // delete the pointer to the vector of entries
        delete (m_cache[bin]);
        m_cache[bin] = NULL;

        m_cacheSize--;
    }
}

// nextEntry - sequential access of the m/z index. returns NULL if there's no more entry left.
SpectraSTLibEntry *SpectraSTMzLibIndex::nextEntry() {

    m_curOffset++;

    while (m_curOffset >= (int) (m_offsets[m_curBin].size())) {
        m_curBin++;
        m_curOffset = 0;
        if (m_curBin >= (unsigned int) m_offsets.size()) {
            return NULL;
        }
    }

    m_libFinPtr->seekg((m_offsets[m_curBin])[m_curOffset]);
    SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
    newEntry->setLibFileOffset((m_offsets[m_curBin])[m_curOffset]);

    return (newEntry);

}

bool SpectraSTMzLibIndex::nextFileOffset(fstream::off_type &offset) {

    m_curOffset++;

    while (m_curOffset >= (int) (m_offsets[m_curBin].size())) {
        m_curBin++;
        m_curOffset = 0;
        if (m_curBin >= (unsigned int) m_offsets.size()) {
            return (false);
        }
    }
    offset = (m_offsets[m_curBin])[m_curOffset];
    return (true);

}

SpectraSTLibEntry *SpectraSTMzLibIndex::thisEntry() {

    if (m_curBin >= (unsigned int) m_offsets.size() || m_curOffset >= (int) (m_offsets[m_curBin].size())) {
        return NULL;
    }

    m_libFinPtr->seekg((m_offsets[m_curBin])[m_curOffset]);
    SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
    newEntry->setLibFileOffset((m_offsets[m_curBin])[m_curOffset]);

    return (newEntry);
}

// sortEntriesByNreps - sorts the entries by number of replicates, so that nextSortedEntry will return the entries in descending
// order of number of replicates.
void SpectraSTMzLibIndex::sortEntriesByNreps() {

    reset();

    m_sortedOffsets = new vector<pair<fstream::off_type, double> >;

    SpectraSTLibEntry *entry = NULL;
    while ((entry = nextEntry())) {
        unsigned int nreps = entry->getNrepsUsed();
        double prob = entry->getProb();

        pair<fstream::off_type, double> p;
        p.first = entry->getLibFileOffset();
        p.second = prob + (double) nreps - 1.0; // aggregate nreps and prob for sorting

        m_sortedOffsets->push_back(p);
        delete (entry);
        entry = NULL;
    }

    sort(m_sortedOffsets->begin(), m_sortedOffsets->end(), SpectraSTMzLibIndex::sortEntriesDesc);


}

// sortEntriesBySN - sorts the entries by signal-to-noise ratio, so that nextSortedEntry will return the entries in descending
// order of S/N.
void SpectraSTMzLibIndex::sortEntriesBySN() {

    reset();

    m_sortedOffsets = new vector<pair<fstream::off_type, double> >;

    for (vector<vector<fstream::off_type> >::iterator bin = m_offsets.begin(); bin != m_offsets.end(); bin++) {

        vector<pair<fstream::off_type, double> > entriesInBin;

        for (vector<fstream::off_type>::iterator en = bin->begin(); en != bin->end(); en++) {

            m_libFinPtr->seekg(*en);
            SpectraSTLibEntry *entry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
            double sn = entry->getPeakList()->calcSignalToNoise();

            pair<fstream::off_type, double> p;
            p.first = (*en);
            p.second = sn;

            entriesInBin.push_back(p);
            delete (entry);

        }

        sort(entriesInBin.begin(), entriesInBin.end(), SpectraSTMzLibIndex::sortEntriesDesc);
        m_sortedOffsets->insert(m_sortedOffsets->end(), entriesInBin.begin(), entriesInBin.end());

    }

}

// sort Entries by offset; add by LWU on Mar 8 2020
void SpectraSTMzLibIndex::sortEntriesByOffsets() {
    reset();
    m_sortedOffsets = new vector<pair<fstream::off_type, double> >;

    for (vector<vector<fstream::off_type> >::iterator bin = m_offsets.begin(); bin != m_offsets.end(); bin++) {
        for (vector<fstream::off_type>::iterator en = bin->begin(); en != bin->end(); en++) {
            pair<fstream::off_type, double> p;
            p.first = (*en);
            p.second = -1.0*p.first;
            m_sortedOffsets->push_back(p);
        }
    }

    sort(m_sortedOffsets->begin(), m_sortedOffsets->end(), SpectraSTMzLibIndex::sortEntriesDesc);

}


// nextSortedEntry - returns the next sorted entry - for now only sorting by Nreps is implemented (see sortEntriesByNreps).
SpectraSTLibEntry *SpectraSTMzLibIndex::nextSortedEntry() {

    if (!m_sortedOffsets) return (nextEntry());

    m_curSortedOffset++;

    if (m_curSortedOffset >= (int) (m_sortedOffsets->size())) {
        return (NULL);
    } else {

        fstream::off_type offset = ((*m_sortedOffsets)[m_curSortedOffset]).first;
        m_libFinPtr->seekg(offset);
        SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
        newEntry->setLibFileOffset(offset);

        return (newEntry);
    }

}

bool SpectraSTMzLibIndex::nextSortedFileOffset(fstream::off_type &offset) {

    if (!m_sortedOffsets) {
        bool success = nextFileOffset(offset);
        return (success);
    }

    m_curSortedOffset++;

    if (m_curSortedOffset >= (int) (m_sortedOffsets->size())) {
        return (false);
    } else {
        offset = ((*m_sortedOffsets)[m_curSortedOffset]).first;
        return (true);
    }

}

SpectraSTLibEntry *SpectraSTMzLibIndex::thisSortedEntry() {

    if (!m_sortedOffsets) {
        return (thisEntry());
    }

    if (m_curSortedOffset >= (int) (m_sortedOffsets->size())) {
        return NULL;
    } else {
        fstream::off_type offset = ((*m_sortedOffsets)[m_curSortedOffset]).first;
        m_libFinPtr->seekg(offset);
        SpectraSTLibEntry *newEntry = new SpectraSTLibEntry(*m_libFinPtr, m_binaryLib);
        newEntry->setLibFileOffset(offset);
        return (newEntry);
    }
}

// reset - goes back to the start, used for the nextEntry and nextSortedEntry sequential readers.
void SpectraSTMzLibIndex::reset() {
    m_curBin = 0;
    m_curOffset = -1;
    m_curSortedOffset = -1;
}

// calcBinNumber - calculates the bin number that holds entries of a particular m/z value.
unsigned int SpectraSTMzLibIndex::calcBinNumber(double mz) {
    if (mz <= 0.0) return (0);
    unsigned int mzInt = (unsigned int) mz;
    if (mzInt < MIN_MZ) {
        return (0);
    }
    if (mzInt > MAX_MZ) {
        return ((unsigned int) (MAX_MZ - MIN_MZ));
    }
    return ((unsigned int) (mzInt - MIN_MZ));

}

// calcBinMz - calculates the m/z value (rounded down to the nearest integer) that corresponds to a particular bin number. 
unsigned int SpectraSTMzLibIndex::calcBinMz(unsigned int binNum) {
    return (binNum + MIN_MZ);
}

// initialize - simply clear the bins and the cache
void SpectraSTMzLibIndex::initialize() {

    for (vector<vector<fstream::off_type> >::iterator oBin = m_offsets.begin(); oBin != m_offsets.end(); oBin++) {
        (*oBin).clear();
    }

    for (vector<vector<SpectraSTLibEntry *> *>::iterator cBin = m_cache.begin(); cBin != m_cache.end(); cBin++) {
        (*cBin) = NULL;
    }
}

// sortEntriesDesc - comparator for descending order
bool SpectraSTMzLibIndex::sortEntriesDesc(pair<fstream::off_type, double> a, pair<fstream::off_type, double> b) {

    return (b.second < a.second);
}

// printDotTimeProfiles
void SpectraSTMzLibIndex::printDotTimeProfiles(ofstream &fout) {

    for (vector<vector<SpectraSTLibEntry *> *>::iterator bin = m_cache.begin(); bin != m_cache.end(); bin++) {
        if (*bin) {
            for (vector<SpectraSTLibEntry *>::iterator en = (*bin)->begin(); en != (*bin)->end(); en++) {
                (*en)->printDotTimeProfile(fout);
            }
        }
    }
}


