#ifndef SPECTRASTSEARCHTASKSTATS_HPP
#define SPECTRASTSEARCHTASKSTATS_HPP

#include "SpectraSTSearch.hpp"
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
 * Class to manage the statistics of a search task 
 */


using namespace std;

class SpectraSTSearchTaskStats {
public:
    SpectraSTSearchTaskStats();

    ~SpectraSTSearchTaskStats();

    void processSearch(SpectraSTSearch *s);

    void logStats();

protected:

    vector<unsigned int> m_numTopHits;
    vector<unsigned int> m_numBadTopHits;
    vector<unsigned int> m_numGoodTopHits;

    vector<unsigned int> m_numDecoyTopHits;
    vector<unsigned int> m_numDecoyBadTopHits;
    vector<unsigned int> m_numDecoyGoodTopHits;

    vector<unsigned int> m_numLowerHits;
    vector<unsigned int> m_numDecoyLowerHits;

    unsigned int m_numTotalLowerHits;
    unsigned int m_numTotalDecoyLowerHits;
    // unsigned int m_numTotalLowerSingletonHits;
    // unsigned int m_numTotalDecoyLowerSingletonHits;

    map<string, int> m_nonDecoyImpure;

};

#endif
