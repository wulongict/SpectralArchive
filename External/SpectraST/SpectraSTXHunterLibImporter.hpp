#ifndef SPECTRASTXHUNTERLIBIMPORTER_HPP_
#define SPECTRASTXHUNTERLIBIMPORTER_HPP_

#include "SpectraSTLibImporter.hpp"

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

/* Class: SpectraSTXHunterLibImporter
 * 
 * Implements a library importer for the .hlf file format (used by X!Hunter).
 * Note that there is no guarantee that X!Hunter libraries will work well with SpectraST!
 * 
 */


class SpectraSTXHunterLibImporter : public SpectraSTLibImporter {

public:

    SpectraSTXHunterLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTXHunterLibImporter();

    virtual void import();


private:

    void readFromFile(string &impFileName);

    unsigned int m_numSkipped;

};

#endif /*SPECTRASTMSPLIBIMPORTER_HPP_*/
