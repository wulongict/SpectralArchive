#include "SpectraSTSearchTaskStats.hpp"
#include "SpectraSTLog.hpp"
#include <sstream>

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

/* Class: SpectraSTSearchTaskStats
 * 
 * Class to manage search task stats
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTSearchTaskStats::SpectraSTSearchTaskStats() :
        m_numTopHits(6, 0),
        m_numBadTopHits(6, 0),
        m_numGoodTopHits(6, 0),
        m_numDecoyTopHits(6, 0),
        m_numDecoyBadTopHits(6, 0),
        m_numDecoyGoodTopHits(6, 0),
        m_numLowerHits(11, 0),
        m_numDecoyLowerHits(11, 0),
        m_numTotalLowerHits(0),
        m_numTotalDecoyLowerHits(0) {
//  m_numTotalLowerSingletonHits(0), 
//  m_numTotalDecoyLowerSingletonHits(0) {
}

// destructor
SpectraSTSearchTaskStats::~SpectraSTSearchTaskStats() {

}

// processSearch - take the search results and store some interesting stats. TO BE EXPANDED
void SpectraSTSearchTaskStats::processSearch(SpectraSTSearch *s) {

    if (s->isNoMatch()) {
        return;
    }

    SpectraSTLibEntry *entry = s->m_candidates[0]->getEntry();
    unsigned int numPeaks = entry->getPeakList()->getNumPeaks();
    unsigned int numAssignedPeaks = entry->getPeakList()->getNumAssignedPeaks();
    double fracAssigned = (double) (numAssignedPeaks) / (double) (numPeaks);

    int ch = entry->getCharge() - 1;
    //unsigned int ch = (unsigned int)(fracAssigned * 5 - 0.01) ;
    if (ch < 0) ch = 0;
    if (ch > 4) ch = 4; // treat all 5+ charges in same bin

    if (s->isDecoy()) {
        map<string, int>::iterator found = m_nonDecoyImpure.find(entry->getName());
        if (found != m_nonDecoyImpure.end()) {
            found->second++;
        } else {
            m_nonDecoyImpure[entry->getName()] = 1;
        }
    }


    m_numTopHits[ch]++;
    if (s->isLikelyGood()) {
        m_numGoodTopHits[ch]++;
    } else if (s->isLikelyBad()) {
        m_numBadTopHits[ch]++;
    }

    if (s->isDecoy(1)) {
        m_numDecoyTopHits[ch]++;
        if (s->isLikelyGood()) {
            m_numDecoyGoodTopHits[ch]++;
        } else if (s->isLikelyBad()) {
            m_numDecoyBadTopHits[ch]++;
        }
    }

    for (unsigned int rank = 2; rank <= 10 && rank <= (unsigned int) (s->m_candidates.size()); rank++) {
        bool isSingleton = s->isSingleton(rank);
        m_numTotalLowerHits++;
//    if (isSingleton) {
//      m_numTotalLowerSingletonHits++;
//    }
        m_numLowerHits[rank]++;
        if (s->isDecoy(rank)) {
            m_numDecoyLowerHits[rank]++;
            m_numTotalDecoyLowerHits++;
//      if (isSingleton) {
//	m_numTotalDecoyLowerSingletonHits++;
//      }
        }
    }

}

// logStats - write the stats in the log file. TO BE EXPANDED
void SpectraSTSearchTaskStats::logStats() {

    stringstream bss;

    bss << "Breakdown: ";
    bss.precision(0);
    for (int ch = 0; ch <= 4; ch++) {
        bss << '+' << ch + 1 << " = " << m_numTopHits[ch] << " ; ";
    }
    g_log->log("SEARCH STATS", bss.str());

    // decoy analysis
    if (m_numTotalDecoyLowerHits > 0) { // if there is any decoys at all
        stringstream decoyss;
        decoyss.precision(4);
        decoyss << "Decoy analysis:";

        for (unsigned int rank = 2; rank <= 10; rank++) {
            if (m_numLowerHits[rank] > 0) {
                decoyss << " D" << rank << " = " << fixed
                        << ((double) m_numDecoyLowerHits[rank] / (double) m_numLowerHits[rank]) << ";";
            }
        }

        for (int ch = 0; ch <= 4; ch++) {
            if (m_numBadTopHits[ch] > 0) {
                decoyss << " D1bad(+" << ch + 1 << ") = " << fixed
                        << ((double) m_numDecoyBadTopHits[ch] / (double) m_numBadTopHits[ch]) << ";";
            }
            if (m_numGoodTopHits[ch] > 0) {
                decoyss << " D1good(+" << ch + 1 << ") = " << fixed
                        << ((double) m_numDecoyGoodTopHits[ch] / (double) m_numGoodTopHits[ch]) << ";";
            }
            if (m_numTopHits[ch] > 0) {
                decoyss << " D1(+" << ch + 1 << ") = " << fixed
                        << ((double) m_numDecoyTopHits[ch] / (double) m_numTopHits[ch]) << ";";
            }
        }

        if (m_numTotalLowerHits > 0) {
            decoyss << " D2-10 = " << fixed << ((double) m_numTotalDecoyLowerHits / (double) m_numTotalLowerHits)
                    << ";";
            //     if (m_numTotalLowerSingletonHits > 0) {
//	decoyss << " D2-10(SINGLETON) = " << fixed << ((double)m_numTotalDecoyLowerSingletonHits / (double)m_numTotalLowerSingletonHits) << ";";
//      }
        }

        g_log->log("SEARCH STATS", decoyss.str());

    }

    for (map<string, int>::iterator i = m_nonDecoyImpure.begin(); i != m_nonDecoyImpure.end(); i++) {
        if (i->second > 5) {
            stringstream badDecoyss;
            badDecoyss << "Frequent decoy hit: " << i->first << " (" << i->second << " times)";
            g_log->log("SEARCH STATS", badDecoyss.str());
        }
    }

}
    
    
