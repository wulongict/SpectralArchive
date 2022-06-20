#include "SpectraSTSearchOutput.hpp"
#include "SpectraSTTxtSearchOutput.hpp"
#include "SpectraSTXlsSearchOutput.hpp"
#include "SpectraSTPepXMLSearchOutput.hpp"
#include "SpectraSTHtmlSearchOutput.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <iostream>



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

extern SpectraSTLog *g_log;

// constructor
SpectraSTSearchOutput::SpectraSTSearchOutput(string outputFileName, string searchFileExt,
                                             SpectraSTSearchParams &searchParams) :
        m_outputFileName(outputFileName),
        m_fout(NULL),
        m_searchFileExt(searchFileExt),
        m_searchParams(searchParams),
        m_instrInfo(NULL),
        m_query("") {

    getPath(outputFileName, m_outputPath);

    // don't open the files yet, wait until it is needed -- when OpenFile() is called.

}

// destructor
SpectraSTSearchOutput::~SpectraSTSearchOutput() {

    if (m_instrInfo) {
        delete (m_instrInfo);
    }
    if (m_fout) {
        delete (m_fout);
    }


}

// openFile - opens the file for output
void SpectraSTSearchOutput::openFile() {

    if (m_fout) return; // already open
    m_fout = new ofstream();
    if (!myFileOpen(*m_fout, m_outputFileName)) {
        g_log->error("SEARCH", "Cannot open file \"" + m_outputFileName + "\" for writing search results. Exiting.");
        g_log->crash();
    }

}

// closeFile - closes the file for output
void SpectraSTSearchOutput::closeFile() {

    if (m_fout) {
        delete (m_fout);
        m_fout = NULL;
    }
}

string SpectraSTSearchOutput::getOutputFileName() {
    return (m_outputFileName);
}

string SpectraSTSearchOutput::getOutputPath() {
    return (m_outputPath);
}

// createSpectraSTSearchOutput - "factory" to create proper output object for different formats. Modify this if you implement
// a new outputter!
SpectraSTSearchOutput *
SpectraSTSearchOutput::createSpectraSTSearchOutput(string searchFileName, SpectraSTSearchParams &searchParams) {

    FileName fn;
    parseFileName(searchFileName, fn);

    string outputDirectory(fn.path);
    if (!(searchParams.outputDirectory.empty())) {
        outputDirectory = searchParams.outputDirectory;
    }

    string ext(searchParams.outputExtension);
    if (ext == "nxls") ext = "xls";
    if (ext == "ntxt") ext = "txt";

    string outputFileName(outputDirectory + fn.name + "." + ext);
    makeFullPath(outputFileName);

    /*
    if (searchParams.saveSpectra) {
      makeDir(outputDirectory + fn.name + ".match/");
      makeDir(outputDirectory + fn.name + ".query/");
    }
    */

    if (searchParams.outputExtension == "txt") {
        return (new SpectraSTTxtSearchOutput(outputFileName, fn.ext, searchParams));
    } else if (searchParams.outputExtension == "ntxt") {
        return (new SpectraSTTxtSearchOutput(outputFileName, fn.ext, searchParams, false));
    } else if (searchParams.outputExtension == "xls") {
        return (new SpectraSTXlsSearchOutput(outputFileName, fn.ext, searchParams));
    } else if (searchParams.outputExtension == "nxls") {
        return (new SpectraSTXlsSearchOutput(outputFileName, fn.ext, searchParams, false));
    } else if (searchParams.outputExtension == "xml" || searchParams.outputExtension == "pepXML" ||
               searchParams.outputExtension == "pep.xml") {
        return (new SpectraSTPepXMLSearchOutput(outputFileName, fn.ext, searchParams, fn.path));
    } else if (searchParams.outputExtension == "html") {
        return (new SpectraSTHtmlSearchOutput(outputFileName, fn.ext, searchParams));
    } else {
        // unknown output type. should never happen.
        return (new SpectraSTPepXMLSearchOutput(outputFileName, fn.ext, searchParams, fn.path));
    }

    return (new SpectraSTPepXMLSearchOutput(outputFileName, fn.ext, searchParams, fn.path));
}
