#ifndef PROGRESSCOUNT_HPP_
#define PROGRESSCOUNT_HPP_

#include <string>


/*

Library       : ProgressCount
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

Institute for Systems Biology, hereby disclaims all copyright interest 
in Spectrast written by Henry Lam.

*/

/* Class: ProgressCount
 * 
 * A helper class that keeps track of the progress of some looping procedure.
 * To use, instantiate a ProgressCount object before the loop (provide the number of iterations if known.)
 * Then call start() with a message to get things started. During the loop, call increment() once every iteration.
 * Call done() when the loop is complete.
 */


using namespace std;

class ProgressCount {

public:
    ProgressCount(bool display, int step, int total = 0);

    virtual ~ProgressCount();

    void start(string msg);

    void increment();

    void done();

    int count() { return (m_count); }

private:

    bool m_display; // whether to display progress to cout or not
    int m_count; // counter
    int m_progress; // the next milestone
    int m_total; // the total number of increment() calls expected
    int m_step; // the number of increment() calls between milestones


};

#endif /*PROGRESSCOUNT_HPP_*/
