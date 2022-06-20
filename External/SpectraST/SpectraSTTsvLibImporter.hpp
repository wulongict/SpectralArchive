#ifndef SPECTRASTTSVLIBIMPORTER_HPP_
#define SPECTRASTTSVLIBIMPORTER_HPP_

#include "SpectraSTLibImporter.hpp"
#include "SpectraSTPeakList.hpp"
#include "Peptide.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "cramp.hpp"
#endif

#include <map>
#include <string>

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/


/* Class: SpectraSTTsvLibImporter
 * 
 * Implements a library importer for a tab-delimited format. 
 * Intended for cases where the identifications are provided in some simple text format,
 * for example, as a list of manually validated results. 
 *
 * Each line of the file must be of the format:
 * <mzXML/mzData/mzML file><TAB><scan num><TAB><identification><TAB><probability><TAB><scores>
 * 
 * where <mzXML/mzData/mzML file> is the relative or absolute path to the file containing the spectrum
 * <scan num> is the scan number of the spectrum
 * <identification> is the identification of that spectrum, in the format R.ABC[160]EFGHIK.L/3.
 * <scores> is of the format "xc=3.0203,dc=0.5600,pb=0.9999,xx=1234" (this field is optional)
 * 
 * Anything in a line after a '#' will be ignored.
 */

using namespace std;

class SpectraSTTsvLibImporter : public SpectraSTLibImporter {

public:

    SpectraSTTsvLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTTsvLibImporter();

    virtual void import();


private:

    // keeps track of the number of mzXML files open
    int m_numMzXMLOpen;

    // the dataset name (usually specified in m_params.datasetName)
    string m_datasetName;

    // hash to map an mzXML file name to its cRamp object
    map<string, cRamp *> m_mzXMLFiles;

    // hash containing all queries to be imported (query name => (path to mzXML file, pointer to SpectraSTLibEntry))
    map<string, pair<vector<string>, SpectraSTLibEntry *> > m_queries;

    // unsigned int m_numImported;
    unsigned int m_numSkipped;
    unsigned int m_numPassedProbCutoff;

    void readFromFile(string &impFileName);

    bool loadSpectrum(string query, string path, SpectraSTLibEntry *entry, string altPath = "");

    bool processSearchHit(string &mzXMLFile, int scanNum, string &peptide, string &searchEngine,
                          string &instrumentType, string &path, double &prob, string &scores, string &protein);


};

#endif /*SPECTRASTTSVLIBIMPORTER_HPP_*/
