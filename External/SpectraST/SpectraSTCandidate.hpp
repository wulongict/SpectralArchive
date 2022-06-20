#ifndef SPECTRASTCANDIDATE_HPP_
#define SPECTRASTCANDIDATE_HPP_

#include "SpectraSTLibEntry.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTSimScores.hpp"
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

/* Class: SpectraSTCandidate
 * 
 * A container class which holds a candidate entry retrieved from the library. 
 * In addition to the actual library entry (which it references), it is also used to hold
 * scoring information after the spectrum is compared to the query spectrum.
 * 
 * 
 */


class SpectraSTCandidate {

public:
    SpectraSTCandidate(SpectraSTLibEntry *entry, SpectraSTSearchParams &params);

    ~SpectraSTCandidate();

    // Simple accessor methods
    SpectraSTLibEntry *getEntry() { return m_entry; };

    double getSortKey() { return m_sortKey; }

    SpectraSTSimScores &getSimScoresRef() { return m_simScores; }


    // Output method that prints the peak list to a file
    void printSpectrum(string spectrumFileName);

    // Checks if the candidate's F value is above the threshold specified in the search parameters
    bool passTopHitFvalThreshold();

    bool passLowerHitsFvalThreshold();


    // Methods for sorting
    void setSortKey(double key) { m_sortKey = key; }

    static bool sortPtrsDesc(SpectraSTCandidate *a, SpectraSTCandidate *b);

    // Methods for recording similarity scores in library entry
    // void recordDotInLibraryEntry();

private:

    // m_entry - Pointer to the library entry object
    // IMPORTANT NOTE:
    // Since the library entries are cached (see SpectraSTMzLibIndex), the library entry object is NOT the property of SpectraSTCandidate. It remains
    // the property of SpectraSTMzLibIndex, and WILL BE RE-USED by other searches (That's precisely the point of caching)!!
    SpectraSTLibEntry *m_entry;

    SpectraSTSearchParams &m_params;

    SpectraSTSimScores m_simScores;

    // m_sortKey - stores the value using which the candidates will be sorted to create the output list.
    // This can be one of the fields of m_simScores, or it can be anything.
    double m_sortKey;


};

#endif /*SPECTRASTCANDIDATE_HPP_*/
