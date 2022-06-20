#ifndef SPECTRASTLIBINDEX_HPP_
#define SPECTRASTLIBINDEX_HPP_

#include "SpectraSTLibEntry.hpp"
#include <iostream>
#include <vector>
#include <list>
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

/* Class: SpectraSTLibIndex
 * 
 * Implements a library index. Can be subclassed for different indexes (on another attribute of the entry, for example).
 * 
 */

using namespace std;

class SpectraSTLibIndex {


public:

    SpectraSTLibIndex(string idxFileName, string idxType);

    SpectraSTLibIndex(string idxFileName, ifstream *libFinPtr, string idxType, bool binaryLib);


    virtual ~SpectraSTLibIndex();

    virtual void insertEntry(SpectraSTLibEntry *entry, fstream::off_type offset) = 0;

    // File I/O methods
    virtual void writeToFile() = 0;

    virtual void readFromFile() = 0;

    string getIdxType() { return (m_idxType); }

    unsigned int getEntryCount() { return (m_entryCount); }

protected:

    // m_libFinPtr - The ifstream of the .splib file. Used to retrieve entries.
    ifstream *m_libFinPtr;

    // m_retrievalMode - TRUE if the index is instantiated for retrieval, FALSE if the index is
    // instantiated for creation.
    bool m_retrievalMode;

    // m_idxFileName - the file name of the index
    string m_idxFileName;

    // m_idxType - the index type (m/z, peptide, etc)
    string m_idxType;

    // m_binaryLib - whether the library is in binary format
    bool m_binaryLib;

    // m_entryCount
    unsigned int m_entryCount;
};


#endif /*SPECTRALIBINDEX_H_*/
