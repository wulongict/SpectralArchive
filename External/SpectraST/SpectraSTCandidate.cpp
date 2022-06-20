#include "SpectraSTCandidate.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <sstream>
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

/* Class: SpectraSTCandidate
 * 
 * A container class which holds a candidate entry retrieved from the library. 
 * In addition to the actual library entry (which it references), it is also used to hold
 * scoring information after the spectrum is compared to the query spectrum.
 * 
 * 
 */

using namespace std;

extern SpectraSTLog *g_log;

// Constructor
SpectraSTCandidate::SpectraSTCandidate(SpectraSTLibEntry *entry, SpectraSTSearchParams &params) :
        m_entry(entry),
        m_params(params) {

}

// Destructor
SpectraSTCandidate::~SpectraSTCandidate() {}


// passTopHitFvalThreshold - Checks if the candidate's F value is above the threshold specified in the search parameters	
bool SpectraSTCandidate::passTopHitFvalThreshold() {
    if (m_params.hitListTopHitFvalThreshold < 0.0000001) {
        return (true);
    } else {
        return (m_simScores.fval >= m_params.hitListTopHitFvalThreshold);
    }
}

// passLowerHitsFvalThreshold - Checks if the candidate's F value is above the threshold specified in the search parameters	
bool SpectraSTCandidate::passLowerHitsFvalThreshold() {
    if (m_params.hitListLowerHitsFvalThreshold < 0.0000001) {
        return (true);
    } else {
        return (m_simScores.fval >= m_params.hitListLowerHitsFvalThreshold);
    }
}


// sortPtrsDesc - comparison method (passed to the sort function)	
bool SpectraSTCandidate::sortPtrsDesc(SpectraSTCandidate *a, SpectraSTCandidate *b) {

    return (a->m_sortKey > b->m_sortKey);

}

// printSpectrum - writes the library entry to the file spectrumFileName
void SpectraSTCandidate::printSpectrum(string spectrumFileName) {

    ofstream fout;
    if (!myFileOpen(fout, spectrumFileName)) {
        g_log->error("SEARCH",
                     "Cannot open file \"" + spectrumFileName + "\" for printing matched library spectrum. Skipped.");
        return;
    }
    m_entry->writeToFile(fout);

}

/*
void SpectraSTCandidate::recordDotInLibraryEntry() {
  
  m_entry->recordDot(m_simScores.dot);
  
}
*/
