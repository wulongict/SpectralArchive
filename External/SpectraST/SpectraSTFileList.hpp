#ifndef SPECTRASTFILELIST_HPP
#define SPECTRASTFILELIST_HPP

#include "SpectraSTCreateParams.hpp"
#include "SpectraSTSearchParams.hpp"

#include <string>
#include <vector>

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

/* Class: SpectraSTFileList
 * 
 * The class that represents a list of files as input to SpectraST. Useful when the number of files are so great that
 * runs over the command line, and for queuing up SpectraST jobs.
 * 
 */

using namespace std;

class SpectraSTFileList {

public:
    // general purpose file list
    SpectraSTFileList(vector<string> &listFileNames);

    // search file list
    SpectraSTFileList(vector<string> &listFileNames, SpectraSTCreateParams &createParams);

    // create file list
    SpectraSTFileList(vector<string> &listFileNames, SpectraSTSearchParams &searchParams);

    ~SpectraSTFileList();

    static int isFileList(vector<string> &fileNames);

    unsigned int getSearchCount();

private:

    // the list of .list files!
    vector<string> m_listFileNames;

    // current mode and default mode (Search or Create)
    char m_mode;
    char m_defaultMode;

    // counter
    unsigned int m_curListFile;

    // the files currently processing
    vector<string> m_currentFiles;

    // the params objects, current and default
    SpectraSTSearchParams m_defaultSearchParams;
    SpectraSTCreateParams m_defaultCreateParams;
    SpectraSTSearchParams m_searchParams;
    SpectraSTCreateParams m_createParams;

    // search count
    unsigned int m_searchCount;

    // private methods
    void parse();

    bool setOptions(string optionLine);

    void issueTask();

};

#endif
