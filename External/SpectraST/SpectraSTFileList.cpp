#include "SpectraSTFileList.hpp"
#include "SpectraSTLib.hpp"
#include "SpectraSTSearchTask.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <iostream>
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

/* Class: SpectraSTFileList
 * 
 * The class that represents a list of files as input to SpectraST. Useful when the number of files are so great that
 * runs over the command line, and for queuing up SpectraST jobs.
 * 
 */

using namespace std;

extern bool g_quiet;
extern SpectraSTLog *g_log;

// Constructor - general lists assuming no knowledge of mode. Will know the mode as the .list files are parsed.
SpectraSTFileList::SpectraSTFileList(vector<string> &listFileNames) :
        m_listFileNames(),
        m_mode('u'),
        m_defaultMode('u'),
        m_curListFile(0),
        m_currentFiles(),
        m_searchParams(),
        m_defaultSearchParams(),
        m_createParams(),
        m_defaultCreateParams(),
        m_searchCount(0) {

    for (vector<string>::iterator i = listFileNames.begin(); i != listFileNames.end(); i++) {
        m_listFileNames.push_back(*i);
        parse();
        m_curListFile++;
    }
}

// Constructor - list for create operation.
SpectraSTFileList::SpectraSTFileList(vector<string> &listFileNames, SpectraSTCreateParams &createParams) :
        m_listFileNames(),
        m_mode('c'),
        m_defaultMode('c'),
        m_curListFile(0),
        m_currentFiles(),
        m_searchParams(),
        m_defaultSearchParams(),
        m_createParams(createParams),
        m_defaultCreateParams(createParams),
        m_searchCount(0) {

    for (vector<string>::iterator i = listFileNames.begin(); i != listFileNames.end(); i++) {
        m_listFileNames.push_back(*i);
        parse();
        m_curListFile++;
    }
}

// Constructor - list for search operation.
SpectraSTFileList::SpectraSTFileList(vector<string> &listFileNames, SpectraSTSearchParams &searchParams) :
        m_listFileNames(),
        m_mode('s'),
        m_defaultMode('s'),
        m_curListFile(0),
        m_currentFiles(),
        m_searchParams(searchParams),
        m_defaultSearchParams(searchParams),
        m_createParams(),
        m_defaultCreateParams(),
        m_searchCount(0) {

    for (vector<string>::iterator i = listFileNames.begin(); i != listFileNames.end(); i++) {
        m_listFileNames.push_back(*i);
        parse();
        m_curListFile++;
    }

}

// Destructor
SpectraSTFileList::~SpectraSTFileList() {
}

// parse - goes down the .list file, set options and issue tasks
// A .list file has the following format: each line is a file (of appropriate type). Empty lines and lines starting with # are ignored.
// Lines starting with a '?' separates one task from another, and anything following the '?' are command-line options for the task
// running on the files after that line. The tasks will be run sequentially, and barring serious errors, it will try to abandon a task
// that produces an error and move onto the next one (instead of truncating the whole queue). Options specified after the '?' override // those at the command-line. Search and Create jobs can be mixed in the same .list file.
//
// Example:
// # let's first create the library
// ? -cNmyLib
// foo.msp
// bar.msp
//
// # now we can search
// ? -sF -sLmyLib.splib
// 0001.mzXML
// 0002.mzXML
// 0003.mzXML
//
// # search another batch with another option
// ? -sF -sM4.0 -sLmyLib.splib
// 0004.mzXML
// 0005.mzXML

void SpectraSTFileList::parse() {

    bool validOptions = true;

    // open the .list file
    ifstream fin;
    if (!myFileOpen(fin, m_listFileNames[m_curListFile])) {
        g_log->error("GENERAL", "Cannot open LIST file \"" + m_listFileNames[m_curListFile] + "\". File skipped.");
        return;
    }

    if (!g_quiet) {
        cout << "Processing LIST file \"" << m_listFileNames[m_curListFile] << "\"..." << endl;
    }

    string line("");
    while (nextLine(fin, line)) {
        line = crop(line);

        if (line.empty() || line[0] == '#') {
            // empty or comment, ignore
            continue;
        }

        if (line[0] == '?') {
            // at a task separator line

            // issue current task
            if (!m_currentFiles.empty() && validOptions) {
                issueTask();
            }

            // start a new task
            m_currentFiles.clear();
            validOptions = setOptions(line);

        } else {

            // add the line to the current file list
            string::size_type dummy = 0;
            m_currentFiles.push_back(
                    nextToken(line, 0, dummy, "#"));  // read only until '#' in case there is trailing comments
        }
    }

    // issue the last task
    if (!m_currentFiles.empty() && validOptions) {
        issueTask();
    }

}

// setOptions - set the options based on the optionLine
bool SpectraSTFileList::setOptions(string optionLine) {

    // first, reset to defaults, which are SpectraST defaults and those that are specified on the command-line.
    m_searchParams = m_defaultSearchParams;
    m_createParams = m_defaultCreateParams;

    m_mode = m_defaultMode;

    // remove '?' and trailing comments
    string::size_type dummy = 0;
    optionLine = nextToken(optionLine, 1, dummy, "#");

    string::size_type pos = 0;
    string option("");
    while (!(option = nextToken(optionLine, pos, pos, " \"\t\r\n", " \t\r\n")).empty()) {
        if (optionLine[pos] == '\"') {
            pos++; // into the quoted string
            string inQuote = nextToken(optionLine, pos, pos, "\"");
            option += inQuote;
            pos++; // into the space following the end quote
        }


        if (option[0] != '-' || (option.size() > 2 && option[1] != 's' && option[1] != 'c')) {
            // doesn't look like an option, ignore
            continue;
        }
        if (m_mode != 'u' && m_mode != option[1]) {
            // inconsistent mode
            g_log->error("GENERAL", "Inconsistent mode (Search or Create) in " + optionLine + " in list file \"" +
                                    m_listFileNames[m_curListFile] + "\". Skipping task.");
            return (false);
        }
        if (option[1] == 's') {
            m_mode = 's';
            m_searchParams.addOption(option.substr(2).c_str());
        } else if (option[1] == 'c') {
            m_mode = 'c';
            m_createParams.addOption(option.substr(2).c_str());
        }
        pos++;
    }

    return (true);

}

// issueTask - runs the search or create tasks
void SpectraSTFileList::issueTask() {

    // create string of output files for console messages
    stringstream filess;
    filess << '\"' << m_currentFiles[0] << '\"';
    if (m_currentFiles.size() > 1) {
        filess << "..\"" << m_currentFiles[m_currentFiles.size() - 1] << '\"';
    }

    if (m_mode == 's') {

        if (!g_quiet) {
            cout << "> Starting SEARCH task on files " << filess.str() << "..." << endl;
        }

        m_searchParams.finalizeOptions();

        if (m_searchParams.libraryFile.empty()) {
            // no library file set!!
            g_log->error("SEARCH",
                         "No library file specified. Cannot proceed with search on files " + filess.str() + ".");
            return;
        }

        if (!g_quiet) {
            cout << "Library File loaded: \"" << m_searchParams.libraryFile << "\"." << endl;
        }

        SpectraSTLib *lib = new SpectraSTLib(m_searchParams.libraryFile, &m_searchParams);
        SpectraSTSearchTask *searchTask = SpectraSTSearchTask::createSpectraSTSearchTask(m_currentFiles, m_searchParams,
                                                                                         lib);

        if (searchTask) {
            searchTask->preSearch();
            searchTask->search();
            searchTask->postSearch();
            m_searchCount += searchTask->getSearchCount();
            delete (searchTask);
        }

        delete (lib);

    } else if (m_mode == 'c') {

        if (!g_quiet) {
            cout << "> Starting CREATE task on files " << filess.str() << "..." << endl;
        }

        m_createParams.finalizeOptions();

        SpectraSTLib *lib = new SpectraSTLib(m_currentFiles, &m_createParams);
        delete (lib);

    } else {
        g_log->error("GENERAL",
                     "Ambiguous mode (Search or Create) for input files \"" + filess.str() + "\". Skipping task.");
        return;
    }

}

// isFileList - checks to see if the vector of fileNames are all .list files
// returns 1 if all files are .list files, 0 if all files are NOT .list files, and -1 if there are both .list and non .list files.
int SpectraSTFileList::isFileList(vector<string> &fileNames) {

    int isList = -1;
    for (vector<string>::iterator f = fileNames.begin(); f != fileNames.end(); f++) {

        FileName fn;
        parseFileName(*f, fn);
        if (fn.ext == ".list") {
            if (isList == 0) {
                // inconsistent
                return (-1);
            } else {
                isList = 1;
            }
        } else {
            if (isList == 1) {
                // inconsistent
                return (-1);
            } else {
                isList = 0;
            }
        }
    }

    return (isList);
}

// getSearchCount - returns the search count.
unsigned int SpectraSTFileList::getSearchCount() {
    return (m_searchCount);
}
