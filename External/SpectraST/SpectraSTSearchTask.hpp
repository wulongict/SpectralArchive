#ifndef SPECTRASTSEARCHTASK_HPP_
#define SPECTRASTSEARCHTASK_HPP_

#include "SpectraSTLib.hpp"
#include "SpectraSTSearchOutput.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTSearchTaskStats.hpp"
#include <vector>
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

/* Class: SpectraSTSearchTask
 * 
 * Abstract base class that implements one search task.
 * Subclass to handle different search file formats. 
 */


using namespace std;

class SpectraSTSearchTask {

public:
    SpectraSTSearchTask(vector<string> &searchFileFullNames, SpectraSTSearchParams &params, SpectraSTLib *lib);

    virtual ~SpectraSTSearchTask();

    virtual void preSearch();

    virtual void postSearch();

    // abstract method; each subclass should implement its own for desired file reading behavior
    virtual void search() = 0;

    unsigned int getSearchCount() { return m_searchCount; }

    static string getQueryStr(string prefix, int scanNum, int charge);

    static SpectraSTSearchTask *
    createSpectraSTSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params, SpectraSTLib *lib);


protected:

    // the names of all the search files in this task
    vector<string> m_searchFileNames;

    // pointer to the library object - NOT a property of this class
    SpectraSTLib *m_lib;

    // the names of all the output files (one per search file)
    //  vector<string> m_outputFileNames;

    // pointers to the output objects (one per search file) - these ARE properties of this class
    vector<SpectraSTSearchOutput *> m_outputs;

    // reference to search params object
    SpectraSTSearchParams &m_params;

    // object to manage stats
    SpectraSTSearchTaskStats m_searchTaskStats;

    // a vector of all the selected queries
    vector<string> m_selectedList;

    // methods to manage selected lists
    void readSelectedListFile();

    bool isInSelectedList(string name);

    // counters and flags
    bool m_searchAll;
    unsigned int m_searchCount;


};

#endif /*SPECTRASTSEARCHTASK_HPP_*/
