#ifndef SPECTRASTSEARCHOUTPUT_HPP_
#define SPECTRASTSEARCHOUTPUT_HPP_

#include "SpectraSTSimScores.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTLibEntry.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "cramp.hpp"
#endif

#include <string>
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

/* Class: SpectraSTSearchOutput
 * 
 * Abstract base class for an outputter. Subclass for different output formats.
 *  
 */



using namespace std;

class SpectraSTSearchOutput {

public:
    SpectraSTSearchOutput(string outputFileName, string searchFileExt, SpectraSTSearchParams &searchParams);

    virtual ~SpectraSTSearchOutput();

    virtual void setInstrInfo(rampInstrumentInfo *instrInfo) { m_instrInfo = instrInfo; }

    // the actual printing methods. subclasses need to implement their own for desired behavior
    virtual void printHeader() = 0;

    virtual void printStartQuery(string query, double precursorMz, int assumedCharge, double retentionTime) = 0;

    virtual void printEndQuery(string query) = 0;

    virtual void
    printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry, SpectraSTSimScores &simScores) = 0;

    virtual void printFooter() = 0;

    virtual void printAbortedQuery(string query, string message) = 0;

    void openFile();

    void closeFile();

    string getOutputFileName();

    string getOutputPath();

    static SpectraSTSearchOutput *
    createSpectraSTSearchOutput(string searchFileName, SpectraSTSearchParams &searchParams);

protected:

    ofstream *m_fout;
    string m_outputFileName;
    string m_outputPath;
    string m_searchFileExt;
    SpectraSTSearchParams &m_searchParams;
    rampInstrumentInfo *m_instrInfo;
    string m_query;

};

#endif /*SPECTRASTSEARCHOUTPUT_HPP_*/
