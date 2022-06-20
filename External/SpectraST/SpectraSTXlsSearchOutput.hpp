#ifndef SPECTRASTXLSSEARCHOUTPUT_HPP_
#define SPECTRASTXLSSEARCHOUTPUT_HPP_

#include "SpectraSTSearchOutput.hpp"


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

/* Class: SpectraSTXlsSearchOutput
 * 
 * Outputter to .xls format
 *    
 */


class SpectraSTXlsSearchOutput : public SpectraSTSearchOutput {

public:
    SpectraSTXlsSearchOutput(string outFullFileName, string inputFileExt, SpectraSTSearchParams &searchParams,
                             bool printHeader = true);

    virtual ~SpectraSTXlsSearchOutput();

    virtual void printHeader();

    virtual void printStartQuery(string query, double precursorMz, int assumedCharge, double retentionTime);

    virtual void printEndQuery(string query) {}

    virtual void printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry, SpectraSTSimScores &simScores);

    virtual void printFooter() {}

    virtual void printAbortedQuery(string query, string message);

private:
    bool m_printHeader;
};

#endif /*SPECTRASTXLSSEARCHOUTPUT_HPP_*/
