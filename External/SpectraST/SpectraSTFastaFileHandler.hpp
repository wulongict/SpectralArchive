#ifndef SPECTRASTFASTAFILEHANDLER_HPP_
#define SPECTRASTFASTAFILEHANDLER_HPP_

#include <vector>
#include <map>
#include <string>
#include <fstream>

/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.08 


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


/* Class: SpectraSTFastaFileHandler
 * 
 * Manages a .fasta protein sequence file.
 * 
 * MAJOR WORK IN PROGRESS
 */

using namespace std;

class SpectraSTFastaFileHandler {

public:

    SpectraSTFastaFileHandler(string fastaFileName);

    virtual ~SpectraSTFastaFileHandler();

    void enableRandomAccess();

    bool findProtein(string &id, string &fullDescr, string &seq);

    bool nextProtein(string &id, string &fullDescr, string &seq);

    void refresh(map<string, vector<pair<string, string> > *> &pepMappings, bool refreshTrypticOnly = false);

private:

    ifstream m_fin;
    string m_fastaFileName;
    map<string, fstream::off_type> *m_offsets;
    fstream::off_type m_nextProteinOffset;

    bool parseNextProtein(string &id, string &fullDescr, string &seq);

};

#endif /*SPECTRASTFASTAFILEHANDLER_HPP_*/
