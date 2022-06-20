#ifndef SPECTRASTMZXMLLIBIMPORTER_HPP_
#define SPECTRASTMZXMLLIBIMPORTER_HPP_

#include "SpectraSTLibImporter.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "cramp.hpp"
#endif

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

/* Class: SpectraSTMzXMLLibImporter
 * 
 * Implements a library importer for the .mzXML file format. 
 * 
 */


class SpectraSTMzXMLLibImporter : public SpectraSTLibImporter {

public:

    SpectraSTMzXMLLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTMzXMLLibImporter();

    virtual void import();

private:

    unsigned int m_numScansInFile;
    unsigned int m_numMissingInFile;
    unsigned int m_numMS1InFile;
    unsigned int m_numFailedFilterInFile;
    unsigned int m_numImportedInFile;
    unsigned int m_numConsensusInFile;
    unsigned int m_numBadConsensusInFile;

    string m_datasetName;

    map<double, vector<SpectraSTLibEntry *> *> m_clusters;

    void readFromFile(string &impFileName);

    void importOne(cRamp *cramp, rampScanInfo *scanInfo, string &prefix);

    void formConsensusEntry(vector<SpectraSTLibEntry *> *cluster);

};

#endif /*SPECTRASTMZXMLLIBIMPORTER_HPP_*/
