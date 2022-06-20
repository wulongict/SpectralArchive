#ifndef SPECTRASTLOG_HPP
#define SPECTRASTLOG_HPP

#include <string>
#include <fstream>
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

/* Class: SpectraSTLog
 * 
 * Manages a log file that takes care of error/warning reporting and various logging. 
 * 
 */


using namespace std;

class SpectraSTLog {
public:
    SpectraSTLog(string logFileName);

    void log(string msg);

    void log(string tag, string msg);

    void log(string tag, string msg, string stamp);

    void error(string tag, string msg);

    void warning(string tag, string msg);

    void crash();

    unsigned int getNumError() { return (m_numError); }

    unsigned int getNumWarning() { return (m_numWarning); }

    void printErrors();

    void printWarnings();

    ~SpectraSTLog();

private:
    ofstream m_fout;
    string m_logFileName;
    bool m_good;
    unsigned int m_numError;
    unsigned int m_numWarning;
    vector<string> m_errors;
    vector<string> m_warnings;
};

#endif
