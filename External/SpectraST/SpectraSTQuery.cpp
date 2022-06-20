#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <sstream>
#include <cmath>
#include <cstdlib>
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

/* Class: SpectraSTQuery
 * 
 * Manages a query.
 *  
 *  
 */


using namespace std;

extern SpectraSTLog *g_log;


// constructor
SpectraSTQuery::SpectraSTQuery(string name, double precursorMz, int charge, string comments,
                               SpectraSTPeakList *peakList) :
        m_name(name),
        m_startScanNum(0),
        m_endScanNum(0),
        m_prefix(""),
        m_precursorMz(precursorMz),
        m_defaultCharge(charge),
        m_comments(comments),
        m_defaultPeakList(peakList),
        m_retentionTime(-1.0),
        m_precursorMzRange(0.0),
        m_precursors(NULL) {

    // See if name is TPP-style. If so, populate the m_prefix, m_startScanNum and m_endScanNum fields
    string prefix("");
    int startScanNum = 0;
    int endScanNum = 0;
    int parsedCharge = 0;
    if (parseQueryTPPStyle(name, prefix, startScanNum, endScanNum, parsedCharge)) {
        m_prefix = prefix;
        m_startScanNum = startScanNum;
        m_endScanNum = endScanNum;
    }

    if (!m_defaultPeakList) {
        m_defaultPeakList = new SpectraSTPeakList(precursorMz, m_defaultCharge);
    }
    if (charge != 0) {
        m_peakLists[charge] = m_defaultPeakList;
    }
}


// constructor - DDA
SpectraSTQuery::SpectraSTQuery(string prefix, int startScanNum, int endScanNum, double precursorMz, int charge,
                               string comments, SpectraSTPeakList *peakList) :
        m_name(constructQueryTPPStyle(prefix, startScanNum, endScanNum, charge)),
        m_prefix(prefix),
        m_startScanNum(startScanNum),
        m_endScanNum(endScanNum),
        m_precursorMz(precursorMz),
        m_defaultCharge(charge),
        m_comments(comments),
        m_defaultPeakList(peakList),
        m_retentionTime(-1.0),
        m_precursorMzRange(0.0),
        m_precursors(NULL) {

    if (!m_defaultPeakList) {
        m_defaultPeakList = new SpectraSTPeakList(precursorMz, charge);
    }
    if (charge != 0) {
        m_peakLists[charge] = m_defaultPeakList;
    }
}

// constructor for DIA/SWATH query (multiple precursors)
SpectraSTQuery::SpectraSTQuery(string prefix, int startScanNum, int endScanNum, double minPrecursorMz,
                               double maxPrecursorMz, vector<double> *precursors, string comments,
                               SpectraSTPeakList *peakList) :
        m_name(constructQueryTPPStyle(prefix, startScanNum, endScanNum, 0)),
        m_prefix(prefix),
        m_startScanNum(startScanNum),
        m_endScanNum(endScanNum),
        m_precursors(precursors),
        m_defaultCharge(0),
        m_comments(comments),
        m_defaultPeakList(peakList),
        m_retentionTime(-1.0) {

    m_precursorMz = (minPrecursorMz + maxPrecursorMz) / 2.0;
    m_precursorMzRange = fabs(maxPrecursorMz - minPrecursorMz);

    if (!m_defaultPeakList) {
        m_defaultPeakList = new SpectraSTPeakList(m_precursorMz, 0);
    }

}

// destructor
SpectraSTQuery::~SpectraSTQuery() {

    if (m_defaultPeakList) delete (m_defaultPeakList);

    if (m_precursors) delete (m_precursors);

    /*
    cout << "Deleting " << m_peakLists.size() << endl;

    for (map<int, SpectraSTPeakList*>::iterator i = m_peakLists.begin(); i != m_peakLists.end(); i++) {
      cout << "  " << i->first << endl;
      if (i->second) {
        cout << "  " << "DEL" << endl;
        delete (i->second);
        i->second = NULL;
      }
    }
    */
}

// printQuerySpectrum - outputs the query spectrum (as an .msp file) for later spectrum plotting
void SpectraSTQuery::printQuerySpectrum(string querySpectrumFileName) {

    ofstream specFout;
    if (!myFileOpen(specFout, querySpectrumFileName)) {
        g_log->error("SEARCH", "Cannot open \"" + querySpectrumFileName + "\" for printing query spectrum.");
        return;
    }

    specFout << "Name: " << m_name << endl;

    if (m_defaultCharge != 0) {
        specFout << "MW: " << m_precursorMz * m_defaultCharge << endl;
    }
    specFout << "PrecursorMZ: " << m_precursorMz << endl;
    specFout << "Comment: " << m_comments << endl;
    m_defaultPeakList->writeToFile(specFout);
    specFout << endl;

}

// constructQueryName - construct a query name of the format <prefix>.<scanNum>.<scanNum>.<charge>
// the <scanNum> fields will be padded with preceding zeros up to 5 characters long
string SpectraSTQuery::constructQueryTPPStyle(string prefix, int startScanNum, int endScanNum, int charge) {

    stringstream ss;
    ss << prefix << '.';
    ss.width(5);
    ss.fill('0');
    ss << right << startScanNum;
    ss << '.';
    ss.width(5);
    ss.fill('0');
    ss << right << endScanNum;
    ss << '.';
    ss << right << charge;
    return (ss.str());

}

bool
SpectraSTQuery::parseQueryTPPStyle(string &query, string &prefix, int &firstScanNum, int &lastScanNum, int &charge) {

    string::size_type lastDotPos = query.rfind('.');
    if (lastDotPos == string::npos) return (false);

    string::size_type secondLastDotPos = query.rfind('.', lastDotPos - 1);  // .<firstScanNum>.<lastScanNum>
    if (secondLastDotPos == string::npos) return (false);

    string::size_type thirdLastDotPos = query.rfind('.', secondLastDotPos - 1);
    if (thirdLastDotPos == string::npos) return (false);

    charge = atoi(query.substr(lastDotPos + 1).c_str());
    firstScanNum = atoi((query.substr(thirdLastDotPos + 1, secondLastDotPos - thirdLastDotPos - 1)).c_str());
    lastScanNum = atoi((query.substr(secondLastDotPos + 1, lastDotPos - secondLastDotPos - 1)).c_str());
    prefix = query.substr(0, thirdLastDotPos);

    return (true);
}

SpectraSTPeakList *SpectraSTQuery::getPeakList(int charge) {
    if (charge == 0) {
        return (m_defaultPeakList);
    } else {
        map<int, SpectraSTPeakList *>::iterator found = m_peakLists.find(charge);
        if (found != m_peakLists.end()) {
            return (found->second);
        } else {
            return (m_defaultPeakList);
        }
    }
}

void SpectraSTQuery::addPossibleCharge(int charge) {

    if (m_defaultCharge == 0) m_defaultCharge = charge;

    map<int, SpectraSTPeakList *>::iterator found = m_peakLists.find(charge);
    if (found != m_peakLists.end()) {
        return;
    } else {
        m_peakLists[charge] = m_defaultPeakList; // don't copy -- assume same binning/scaling for all charges

        //    m_peakLists[charge] = new SpectraSTPeakList(*m_defaultPeakList);
//    if (m_defaultPeakList->getFragType() == "ETD") {
//      m_peakLists[charge]->setParentCharge(charge, true);
//    } else {
        m_peakLists[charge]->setParentCharge(charge, false);
//    }
    }
}

bool SpectraSTQuery::isPossibleCharge(int charge) {
    return (charge == 0 || m_defaultCharge == 0 || m_peakLists.find(charge) != m_peakLists.end());
}


bool SpectraSTQuery::matchPrecursorMz(double mz, double tolerance, int charge, unsigned int numIsotopes) {

    if (charge == 0) charge = 1; // just assume +1 if library entry has no charge
    if (tolerance >= 0.5) numIsotopes = 0; // no need to consider isotope error if the tolerance is already so big

    if (!m_precursors) {
        // don't have list of precursors, try whole range

        double minMz = m_precursorMz - 0.5 * m_precursorMzRange - tolerance;
        double maxMz = m_precursorMz + 0.5 * m_precursorMzRange + tolerance;

        for (unsigned int n = 0; n <= numIsotopes; n++) {
            double testMz = mz + (double) n / (double) charge;
            if (testMz >= minMz && testMz <= maxMz) {
                return (true);
            }
        }
        return (false);
    }

    // consider list of precursors
    for (vector<double>::iterator i = m_precursors->begin(); i != m_precursors->end(); i++) {

        double minMz = (*i) - tolerance;
        double maxMz = (*i) + tolerance;

        for (unsigned int n = 0; n <= numIsotopes; n++) {
            double testMz = mz + (double) n / (double) charge;
            if (testMz >= minMz && testMz <= maxMz) {
                return (true);
            }
        }
    }

    return (false);

}
