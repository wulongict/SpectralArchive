#ifndef SPECTRASTPEPXMLLIBIMPORTER_HPP_
#define SPECTRASTPEPXMLLIBIMPORTER_HPP_


#include "SpectraSTLibImporter.hpp"
#include "SpectraSTPeakList.hpp"
#include "Peptide.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "Parsers/mzParser/cramp.hpp"
#endif

#include <map>
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


/* Class: SpectraSTPepXMLLibImporter
 * 
 * Implements a library importer for the .pepXML file format. Takes a pepXML file,
 * finds all queries with peptide identification over a certain confidence threshold, goes
 * back to the .mzXML file to retrieve the spectrum, and creates a library of these spectra.
 * 
 * Note that this importer will import all duplicate spectra of the same peptide.
 * 
 * If multiple instances of the same spectra are present (e.g. as separate queries for different charge state
 * in the same file, or as queries in different searches), the importer will resolve the conflict by
 * only taking the one with the highest probability.
 * 
 * 
 * 
 */

using namespace std;

class SpectraSTPepXMLLibImporter : public SpectraSTLibImporter {

public:

    SpectraSTPepXMLLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTPepXMLLibImporter();

    virtual void import();

    static bool parseQuery(string &query, string &baseName, int &firstScanNum, int &lastScanNum);


    static double linearRegress(vector<double> &vX, vector<double> &vY, pair<double, double> &coeff);

    static int calcLinearRegressionResiduals(vector<double> &vX, vector<double> &vY, pair<double, double> &coeff,
                                             vector<double> &residuals);


private:

    // keeps track of the number of mzXML files open
    int m_numMzXMLOpen;

    // the dataset name (usually specified in m_params.datasetName)
    string m_datasetName;

    // hash to map an mzXML file name to its cRamp object
    map<string, cRamp *> m_mzXMLFiles;

    // hash containing all queries to be imported (query name => (path to mzXML file, pointer to SpectraSTLibEntry))
    map<string, pair<vector<string>, SpectraSTLibEntry *> > m_queries;

    unsigned int m_numQuery;
    unsigned int m_numQueryWithPeptideProphetProb;
    unsigned int m_numQueryWithiProphetProb;
    unsigned int m_numQueryWithPercolatorProb;
    unsigned int m_numQueryPassedProbCutoff;
    unsigned int m_numSkipped;
    double m_probCutoff;

    map<string, pair<double, vector<double> > > *m_rtLandmarks;
    map<string, pair<double, double> > m_normalizeRTCoeff;
    map<string, map<double, double> > m_normalizeRTPivots;

    void readFromFile(string &impFileName);

    bool processSearchHit(string &searchHit, string &query, int charge, string &rtstr, string &searchEngine,
                          string &instrumentType, string &fragType, string &path, string &altPath,
                          string &experimentLabel, string &imbstr);

    bool loadSpectrum(string query, string path, vector<SpectraSTLibEntry *> &entry, string altPath = "");

    // this is deprecated
    void setStaticMods(string aas, string masses, string variables);

    void addModificationsToPeptide(string &peptide, string &modInfo);

    bool readOneScan(cRamp *cramp, int scanNum, SpectraSTLibEntry *entry, string &baseName, string &fullFileName,
                     bool silent = false, double minMz = 0.0, double maxMz = 999999.9);

    cRamp *openCRamp(string &baseName, string &path, string &altPath);

    void extractBracket(string &query, string &path, string &altPath);

    double getProbCutoffFromFDR(double desiredFDR, string &errorsStr, string &minProbsStr, double &actualFDR);


    double normalizeRT(string &baseName, double rRT);

    double calculateIRTLinearRegression(double rRT, pair<double, double> &coeff);

    double calculateIRTLinearInterpolation(double rRT, map<double, double> &pivots);

    void readRTLandmarksFromFile();

    void populateRTLandmarks();

    bool doNormalizeRTLinearRegression(pair<double, double> &coeff, string &baseName);

    bool doNormalizeRTLinearInterpolation(map<double, double> &pivots, string &baseName);

};

#endif /*SPECTRASTMSPLIBIMPORTER_HPP_*/
