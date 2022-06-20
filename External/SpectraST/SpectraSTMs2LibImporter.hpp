#ifndef SPECTRASTMS2LIBIMPORTER_HPP_
#define SPECTRASTMS2LIBIMPORTER_HPP_

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

/* Class: SpectraSTMs2LibImporter
 * 
 * Implements a library importer for the .ms2 file format (used by BiblioSpec).
 * Note that there is no guarantee that BiblioSpec libraries will work well with SpectraST!
 * 
 */


class SpectraSTMs2LibImporter : public SpectraSTLibImporter {

public:

    SpectraSTMs2LibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTMs2LibImporter();

    virtual void import();


private:

    void readFromFile(string &impFileName);

};

#endif /*SPECTRASTMS2LIBIMPORTER_HPP_*/
