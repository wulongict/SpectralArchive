#ifndef SPECTRASTQUERY_HPP_
#define SPECTRASTQUERY_HPP_

#include "SpectraSTPeakList.hpp"


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
 */

class SpectraSTQuery {

public:

    // Generic
    SpectraSTQuery(string name, double precursorMz, int charge, string comments, SpectraSTPeakList *peakList);

    // DDA query from mzXML
    SpectraSTQuery(string prefix, int startScanNum, int endScanNum, double precursorMz, int charge, string comments,
                   SpectraSTPeakList *peakList);

    // DIA query from mzXML
    SpectraSTQuery(string prefix, int startScanNum, int endScanNum, double minPrecursorMz, double maxPrecursorMz,
                   vector<double> *precursors, string comments, SpectraSTPeakList *peakList);

    virtual ~SpectraSTQuery();

    string getName() { return m_name; }

    double getPrecursorMz() { return m_precursorMz; }

    vector<double> *getPrecursors() { return m_precursors; }

    int getScanNum() { return m_startScanNum; }

    int getStartScanNum() { return m_startScanNum; }

    int getEndScanNum() { return m_startScanNum; }

    string getPrefix() { return m_prefix; }

    int getDefaultCharge() { return m_defaultCharge; }

    string getComments() { return m_comments; }

    SpectraSTPeakList *getPeakList(int charge = 0);

    double getRetentionTime() { return m_retentionTime; }

    void setRetentionTime(double rt) { if (rt >= 0.0) m_retentionTime = rt; }

    void printQuerySpectrum(string querySpectrumFileName);

    void addPossibleCharge(int charge);

    bool isPossibleCharge(int charge);

    bool matchPrecursorMz(double mz, double tolerance, int charge, unsigned int numIsotopes = 0);

    static string constructQueryTPPStyle(string prefix, int startScanNum, int endScanNum, int charge);

    static bool parseQueryTPPStyle(string &query, string &prefix, int &firstScanNum, int &lastScanNum, int &charge);

private:

    // fields
    string m_name;
    string m_prefix;
    int m_startScanNum;
    int m_endScanNum;

    double m_precursorMz;
    int m_defaultCharge;
    string m_comments;

    double m_retentionTime;
    double m_precursorMzRange;

    vector<double> *m_precursors;

    map<int, SpectraSTPeakList *> m_peakLists;

    // The peak list is owned by this query
    SpectraSTPeakList *m_defaultPeakList;

};

#endif /*SPECTRASTQUERY_HPP_*/
