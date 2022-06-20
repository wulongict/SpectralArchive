#include "SpectraSTLibIndex.hpp"
#include "SpectraSTMzLibIndex.hpp"
#include "SpectraSTPeptideLibIndex.hpp"
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

/* Class: SpectraSTLibIndex
 * 
 * Implements a library index on the precursor m/z value. Can be subclassed for different implementation
 * of the same index, or for a different index (on another attribute of the entry, for example).
 * 
 * NOTE ON CACHING: This implementation of the index allows caching of the recently retrieved entries in
 * memory. If the list of query spectra is sorted by its precursor m/z prior to searching, this
 * allows the search to proceed orders-of-magnitude faster, since the bottleneck of the search lies in
 * accessing the file for the candidate entries. Please see the cache() function for more information on the implementation. 
 * 
 */

using namespace std;

// constructor for creating
SpectraSTLibIndex::SpectraSTLibIndex(string idxFileName, string idxType) :
        m_idxFileName(idxFileName),
        m_libFinPtr(NULL),
        m_idxType(idxType),
        m_retrievalMode(false),
        m_binaryLib(false),
        m_entryCount(0) {

}

// constructor for retrieval
SpectraSTLibIndex::SpectraSTLibIndex(string idxFileName, ifstream *libFinPtr, string idxType, bool binaryLib) :
        m_idxFileName(idxFileName),
        m_libFinPtr(libFinPtr),
        m_idxType(idxType),
        m_retrievalMode(true),
        m_binaryLib(binaryLib),
        m_entryCount(0) {

}

// destructor - Since all the entries stored in the cache remain property of the index, they
// have to be deleted here.
SpectraSTLibIndex::~SpectraSTLibIndex() {}
