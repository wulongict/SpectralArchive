#include "SpectraSTSearchTask.hpp"
#include "SpectraSTMspSearchTask.hpp"
#include "SpectraSTMzXMLSearchTask.hpp"
#include "SpectraSTDtaSearchTask.hpp"
#include "SpectraSTMgfSearchTask.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <string>
#include <iostream>
#include <fstream>
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

/* Class: SpectraSTSearchTask
 * 
 * Abstract base class that implements one search task.
 * Subclass to handle different search file formats. 
 */

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor - all subclasses must call this constructor
SpectraSTSearchTask::SpectraSTSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                         SpectraSTLib *lib) :
        m_searchFileNames(searchFileNames),
        m_params(params),
        m_lib(lib),
        m_outputs(),
        m_searchCount(0),
        m_searchTaskStats(),
        m_selectedList(),
        m_searchAll(true) {

    for (vector<string>::iterator f = searchFileNames.begin(); f != searchFileNames.end(); f++) {

        SpectraSTSearchOutput *output = SpectraSTSearchOutput::createSpectraSTSearchOutput(*f, m_params);
        m_outputs.push_back(output);

    }

    // if a selected list is specified (i.e. search only a subset of the queries), read it now.
    if (!(m_params.filterSelectedListFileName.empty())) {
        readSelectedListFile();
        m_searchAll = false;
    }

}

// destructor
SpectraSTSearchTask::~SpectraSTSearchTask() {

    /*

    // tgz the saved spectra if asked
    if (m_params.saveSpectra && m_params.tgzSavedSpectra) {
      for (vector<string>::iterator f = m_outputFileNames.begin(); f != m_outputFileNames.end(); f++) {
        FileName fn;
        parseFileName(*f, fn);

        myTgz(fn.path + fn.name + ".match.tgz", fn.name + ".match/", true);
        removeDir(fn.path + fn.name + ".match/");

        myTgz(fn.path + fn.name + ".query.tgz", fn.name + ".query/", true);
        removeDir(fn.path + fn.name + ".query/");
      }
    }

    */

    if (g_verbose) {
        cout << endl;
    }

    // print out finishing messages
    if (!g_quiet) {

        if (m_searchFileNames.size() == 1) {
            cout << "Finished searching \"" << m_searchFileNames[0] << "\" (" << m_searchCount << " spectra searched.)"
                 << endl;
            cout << "Output written to \"" << m_outputs[0]->getOutputFileName() << "\"." << endl;

            /*
            if (m_params.saveSpectra && !m_params.tgzSavedSpectra) {
          FileName fn;
          parseFileName(m_outputFileNames[0], fn);
          cout << "Query spectra saved to directory \"" << fn.path << fn.name << ".query/\"." << endl;
          cout << "Matched library spectra saved to directory \"" << fn.path << fn.name << ".match/\"." << endl;
            } else if (m_params.saveSpectra && m_params.tgzSavedSpectra) {
          FileName fn;
          parseFileName(m_outputFileNames[0], fn);
          cout << "Query spectra saved to tgz file \"" << fn.path << fn.name << ".query.tgz\"." << endl;
          cout << "Matched library spectra saved to tgz file \"" << fn.path << fn.name << ".match.tgz\"." << endl;
            }
            */

        } else {
            // multiple files
            cout << "Finished searching \"" << m_searchFileNames[0] << "\" ... \"";
            cout << m_searchFileNames[m_searchFileNames.size() - 1] << "\" (" << m_searchCount << " spectra searched.)"
                 << endl;
            cout << "Output written to \"" << m_outputs[0]->getOutputFileName() << "\" ... ";
            cout << m_outputs[m_outputs.size() - 1]->getOutputFileName() << "\"." << endl;

            /*
            if (m_params.saveSpectra && !m_params.tgzSavedSpectra) {
          FileName fn1;
          parseFileName(m_outputFileNames[0], fn1);
          FileName fn2;
          parseFileName(m_outputFileNames[m_outputFileNames.size() - 1], fn2);
          cout << "Query spectra saved to directories \"" << fn1.path << fn1.name << ".query/\" ... \"";
          cout << fn2.path << fn2.name << ".query/\"." << endl;
          cout << "Matched library spectra saved to directories \"" << fn1.path << fn1.name << ".match/\" ... \"";
          cout << fn2.path << fn2.name << ".match/\"." << endl;

            } else if (m_params.saveSpectra && m_params.tgzSavedSpectra) {
          FileName fn1;
          parseFileName(m_outputFileNames[0], fn1);
          FileName fn2;
          parseFileName(m_outputFileNames[m_outputFileNames.size() - 1], fn2);
          cout << "Query spectra saved to tgz files \"" << fn1.path << fn1.name << ".query.tgz\" ... \"";
          cout << fn2.path << fn2.name << ".query.tgz\"." << endl;
          cout << "Matched library spectra saved to tgz files \"" << fn1.path << fn1.name << ".match.tgz\" ... \"";
          cout << fn2.path << fn2.name << ".match.tgz\"." << endl;
            }
            */
        }


    }

    // delete the output objects
    for (vector<SpectraSTSearchOutput *>::iterator op = m_outputs.begin(); op != m_outputs.end(); op++) {
        if (*op) {
            delete (*op);
        }
    }


}

// preSearch - called before search() is called. if any groundwork needs to be laid before any search,
// this is where to do it.
void SpectraSTSearchTask::preSearch() {

    // not doing anything
}

// preSearch - called after search() is called. if any finishing touch needs to be done after any search,
// this is where to do it.
void SpectraSTSearchTask::postSearch() {

    // not doing anything
}

// readSelectedListFile - reads in the selected queries. should be common for any search file format.
void SpectraSTSearchTask::readSelectedListFile() {

    if (!g_quiet) {
        cout << "Selected query list file loaded: \"" << m_params.filterSelectedListFileName << "\"." << endl;
    }

    ifstream fin;
    if (!myFileOpen(fin, m_params.filterSelectedListFileName)) {
        g_log->error("SEARCH", "Cannot open query list file \"" + m_params.filterSelectedListFileName +
                               "\". No spectrum will be searched.");
        return;
    }

    string line;
    while (nextLine(fin, line, "", "")) {
        string::size_type pos;
        string s = nextToken(line, 0, pos);
        if (!s.empty()) {
            m_selectedList.push_back(s);
            if (g_verbose) {
                cout << "\tQuery selected: " << s << endl;
            }
        }
    }
}

// isInSelectedList - tries to find the query name in the selected list and returns true if found.
bool SpectraSTSearchTask::isInSelectedList(string name) {

    vector<string>::iterator i;
    for (i = m_selectedList.begin(); i != m_selectedList.end(); i++) {

        if (name == *i) {
            return (true);
        }
    }
    return (false);
}


// createSpectraSTSearchTask - factory method to create the proper search task object for different search file formats. Modify
// if a subclass is implemented!
SpectraSTSearchTask *
SpectraSTSearchTask::createSpectraSTSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                               SpectraSTLib *lib) {

    if (searchFileNames.empty()) return (NULL);

    // need to ensure all the search files have same type. check the extension of the first search file
    string firstExt("");
    char *firstRampExt = rampValidFileType(searchFileNames[0].c_str());
    if (!firstRampExt) {
        string::size_type firstDummy = getExtension(searchFileNames[0], firstExt);
    } else {
        // if RAMP thinks this is valid, it could be mzXML or mzData. but we make no
        // distinction between them, since they will be read by RAMP anyway
        firstExt = ".mzXML";
    }

    vector<string> goodSearchFileNames;
    goodSearchFileNames.push_back(searchFileNames[0]);

    // mixed file formats not allowed, so we check each search file to make sure it has the same extension as the first one
    // if any of them isn't, they are simply skipped in the search and won't affect the rest of the search files
    if (searchFileNames.size() > 1) {
        for (vector<string>::iterator i = searchFileNames.begin() + 1; i != searchFileNames.end(); i++) {
            string ext("");
            char *rampExt = rampValidFileType((*i).c_str());
            if (!rampExt) {
                string::size_type dummy = getExtension(*i, ext);
            } else {
                ext = ".mzXML"; // allowing mixed formats that are all RAMP-compatible
            }
            if (ext != firstExt) {
                g_log->error("SEARCH", "Different search file format for \"" + *i + "\". File skipped.");
            } else {
                goodSearchFileNames.push_back(*i);
            }
        }
    }

    // create the search task objects
    if (firstExt == ".msp" || firstExt == ".sptxt" || firstExt == ".MSP" || firstExt == ".SPTXT") {
        return (new SpectraSTMspSearchTask(goodSearchFileNames, params, lib));
    } else if (firstExt == ".mzXML") { // Note that this stands for any of the RAMP valid types (see above)
        return (new SpectraSTMzXMLSearchTask(goodSearchFileNames, params, lib));
    } else if (firstExt == ".dta" || firstExt == ".DTA") {
        return (new SpectraSTDtaSearchTask(goodSearchFileNames, params, lib));

    } else if (firstExt == ".mgf" || firstExt == ".MGF") {
        return (new SpectraSTMgfSearchTask(goodSearchFileNames, params, lib));

    } else {
        return (NULL);
    }
}




