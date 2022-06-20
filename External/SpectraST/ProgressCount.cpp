#include "ProgressCount.hpp"
#include <iostream>


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
 * A helper class that keeps track of the progress of some looping procedure
 */



// Constructor - display is whether or not to display the progress, step is the number
// of increment() calls before a progress report is displayed, total is the expected total
// number of increment() calls. If total != 0, a percentage of progress will be displayed -
// the step in this case refers to the number of 10% before a progress report; if
// the total number is not known, so only the absolute number of increment() calls so far is reported. 
ProgressCount::ProgressCount(bool display, int step, int total) :
        m_display(display),
        m_count(0),
        m_progress(0),
        m_total(total),
        m_step(step) {

    if (m_total != 0 && m_step > 5) {
        // won't allow step > 50%
        m_step = 1;
    }

    m_progress = m_step;
}

// start - starts the counter. The msg will preceed the actual progress counting, which
// is of the format <msg>...10%...20%...30%... etc.
void ProgressCount::start(string msg) {
    if (m_display) {
        cout << msg << "...";
        cout.flush();
    }
}

ProgressCount::~ProgressCount() {
}

// increment - increment the counter, check if it reaches a milestone for reporting
void ProgressCount::increment() {

    m_count++;
    if (m_total != 0) {
        // know the total, report percentage
        if (m_count > m_progress * 0.1 * m_total) {
            if (m_display) {
                cout << m_progress * 10 << "%...";
                cout.flush();
            }
            m_progress += m_step;
        }
    } else {
        // don't know the total, report absolute number
        if (m_count > m_progress) {
            if (m_display) {
                cout << m_progress << "...";
                cout.flush();
            }
            m_progress += m_step;
        }
    }
}

// done - finishing, print "DONE!" message.
void ProgressCount::done() {
    if (m_display) {
        cout << "DONE!" << endl;
    }
}
