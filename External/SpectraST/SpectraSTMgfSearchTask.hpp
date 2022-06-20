#ifndef SPECTRASTMGFSEARCHTASK_HPP_
#define SPECTRASTMGFSEARCHTASK_HPP_

#include "SpectraSTSearchTask.hpp"

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

/* Class: SpectraSTMgfSearchTask
 *
 * Subclass of SpectraSTSearchTask that handles .mgf file formats
 *
 */

class SpectraSTMgfSearchTask : public SpectraSTSearchTask {

public:
    SpectraSTMgfSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params, SpectraSTLib *lib);

    virtual ~SpectraSTMgfSearchTask();

    virtual void search();

private:
    void searchOneFile(unsigned int fileIndex);


};

#endif /*SPECTRASTMGFSEARCHTASK_HPP_*/

