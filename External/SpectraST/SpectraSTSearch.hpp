#ifndef SPECTRASTSEARCH_HPP_
#define SPECTRASTSEARCH_HPP_

#include "SpectraSTLib.hpp"
#include "SpectraSTCandidate.hpp"
#include "SpectraSTSearchOutput.hpp"
#include "SpectraSTQuery.hpp"
// #include "SpectraSTSearchTaskStats.hpp"

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

/* Class: SpectraSTSearch
 * 
 * Implements the search for one query spectrum.
 *  
 */

class SpectraSTSearch {

public:

    friend class SpectraSTSearchTaskStats;

    friend class SpectraSTMzXMLSearchTask;

    friend class SpectraSTMgfSearchTask;

    SpectraSTSearch(SpectraSTQuery *query,
                    SpectraSTSearchParams &params,
                    SpectraSTSearchOutput *output);

    virtual ~SpectraSTSearch();

    void search(SpectraSTLib *lib);

    void print();

    SpectraSTLibEntry *getTopHit() { return (m_candidates.size() > 0 ? m_candidates[0]->getEntry() : NULL); }

    bool isLikelyGood();

    bool isLikelyBad();

    bool isNoMatch() { return (m_candidates.size() == 0); }

    bool isDecoy(unsigned int rank = 1);

    bool isSingleton(unsigned int rank = 1);


private:

    SpectraSTQuery *m_query;

    // the search params
    SpectraSTSearchParams &m_params;

    // the candidates
    vector<SpectraSTCandidate *> m_candidates;

    // the output object responsible for printing the search results
    SpectraSTSearchOutput *m_output;

//  void calcDeltaSimpleDots();
//  void calcHitsStats();

    bool isWithinPrecursorTolerance(SpectraSTLibEntry *entry);

    void detectHomologs();

    void finalizeScoreUsePValue();

    void finalizeScoreLinearCombination(unsigned int numHits);

    static double upperIncompleteGamma(double x, double a);

    static double calcAndersonDarlingScore(vector<double> &allPValues);


};

#endif /*SPECTRASTSEARCH_HPP_*/
