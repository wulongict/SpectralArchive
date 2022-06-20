#include "Predicate.hpp"
#include <math.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#ifndef STANDALONE_LINUX
#include "common/sysdepend.h" // define round() for msvc
#endif

/*
 
Library       : Predicate
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

/* Class: Predicate
 * 
 * Manages a predicate that needs to be repeated tested.
 */


using namespace std;


// Constructor for a string predicate
Predicate::Predicate(string refName, string oper, string refValue, bool defaultBool) :
        m_refName(refName),
        m_oper(oper),
        m_type('S'),
        m_refValueS(refValue),
        m_refValueI(0),
        m_refValueD(0.0),
        m_default(defaultBool) {

}

// Constructor for a int predicate
Predicate::Predicate(string refName, string oper, int refValue, bool defaultBool) :
        m_refName(refName),
        m_oper(oper),
        m_type('I'),
        m_refValueS(""),
        m_refValueI(refValue),
        m_refValueD(0.0),
        m_default(defaultBool) {

}

// Constructor for a double predicate
Predicate::Predicate(string refName, string oper, double refValue, bool defaultBool) :
        m_refName(refName),
        m_oper(oper),
        m_type('D'),
        m_refValueS(""),
        m_refValueI(0),
        m_refValueD(refValue),
        m_default(defaultBool) {

}

// Destructor
Predicate::~Predicate() {

}

// compRef - compare to reference value
bool Predicate::compRef(string value) {
    if (m_type != 'S') {
        return (m_default);
    }

    if (m_oper == "==") {
        return (value == m_refValueS);
    } else if (m_oper == "!=") {
        return (value != m_refValueS);
    } else if (m_oper == "<=") {
        return (value <= m_refValueS);
    } else if (m_oper == ">=") {
        return (value >= m_refValueS);
    } else if (m_oper == "<") {
        return (value < m_refValueS);
    } else if (m_oper == ">") {
        return (value > m_refValueS);
    } else if (m_oper == "=~") {
        return (value.find(m_refValueS, 0) != string::npos);

        // TODO: may want to have real regular expression support later
    } else if (m_oper == "!~") {
        return (value.find(m_refValueS, 0) == string::npos);
    } else {
        return (m_default);
    }
}

// compRef - compare to reference value
bool Predicate::compRef(int value) {
    if (m_type != 'I') {
        return (m_default);
    }

    if (m_oper == "==") {
        return (value == m_refValueI);
    } else if (m_oper == "!=") {
        return (value != m_refValueI);
    } else if (m_oper == "<=") {
        return (value <= m_refValueI);
    } else if (m_oper == ">=") {
        return (value >= m_refValueI);
    } else if (m_oper == "<") {
        return (value < m_refValueI);
    } else if (m_oper == ">") {
        return (value > m_refValueI);
    } else if (m_oper == "=~") {
        return (m_default);
    } else if (m_oper == "!~") {
        return (m_default);
    } else {
        return (m_default);
    }
}

// compRef - compare to reference value
bool Predicate::compRef(double value) {
    if (m_type != 'D') {
        return (m_default);
    }
    if (m_oper == "==") {
        return (value == m_refValueD);
    } else if (m_oper == "!=") {
        return (value != m_refValueD);
    } else if (m_oper == "<=") {
        return (value <= m_refValueD);
    } else if (m_oper == ">=") {
        return (value >= m_refValueD);
    } else if (m_oper == "<") {
        return (value < m_refValueD);
    } else if (m_oper == ">") {
        return (value > m_refValueD);
    } else if (m_oper == "=~") {
        return (m_default);
    } else if (m_oper == "!~") {
        return (m_default);
    } else {
        return (m_default);
    }
}

// comp - compare two values
bool Predicate::comp(string value1, string value2) {
    if (m_type != 'S') {
        return (m_default);
    }

    if (m_oper == "==") {
        return (value1 == value2);
    } else if (m_oper == "!=") {
        return (value1 != value2);
    } else if (m_oper == "<=") {
        return (value1 <= value2);
    } else if (m_oper == ">=") {
        return (value1 >= value2);
    } else if (m_oper == "<") {
        return (value1 < value2);
    } else if (m_oper == ">") {
        return (value1 > value2);
    } else if (m_oper == "=~") {
        return (value1.find(value2, 0) != string::npos);
    } else if (m_oper == "!~") {
        return (value1.find(value2, 0) == string::npos);
    } else {
        return (m_default);
    }
}

// comp - compare two values
bool Predicate::comp(int value1, int value2) {
    if (m_type != 'I') {
        return (m_default);
    }

    if (m_oper == "==") {
        return (value1 == value2);
    } else if (m_oper == "!=") {
        return (value1 != value2);
    } else if (m_oper == "<=") {
        return (value1 <= value2);
    } else if (m_oper == ">=") {
        return (value1 >= value2);
    } else if (m_oper == "<") {
        return (value1 < value2);
    } else if (m_oper == ">") {
        return (value1 > value2);
    } else if (m_oper == "=~") {
        return (m_default);
    } else if (m_oper == "!~") {
        return (m_default);
    } else {
        return (m_default);
    }
}

// comp - compare two values
bool Predicate::comp(double value1, double value2) {
    if (m_type != 'D') {
        return (m_default);
    }

    if (m_oper == "==") {
        return (value1 == value2);
    } else if (m_oper == "!=") {
        return (value1 != value2);
    } else if (m_oper == "<=") {
        return (value1 <= value2);
    } else if (m_oper == ">=") {
        return (value1 >= value2);
    } else if (m_oper == "<") {
        return (value1 < value2);
    } else if (m_oper == ">") {
        return (value1 > value2);
    } else if (m_oper == "=~") {
        return (m_default);
    } else if (m_oper == "!~") {
        return (m_default);
    } else {
        return (m_default);
    }
}

// changeType - attempt to change the predicate to new type.
void Predicate::changeType(char newType) {
    if (newType == m_type) {
        return;
    }

    if (newType == 'D' && m_type == 'I') {
        // int to double, simply do a type cast
        m_refValueD = (double) m_refValueI;
    } else if (newType == 'I' && m_type == 'D') {
        // double to int, do a rounding
        m_refValueI = (int) (round(m_refValueD));
    } else if (newType == 'S' && m_type == 'D') {
        // double to string, convert double to string representation
        stringstream ss;
        ss << m_refValueD;
        m_refValueS = ss.str();
    } else if (newType == 'S' && m_type == 'I') {
        // int to string, convert double to string representation
        stringstream ss;
        ss << m_refValueI;
        m_refValueS = ss.str();
    } else if (newType == 'D' && m_type == 'S') {
        // string to double, do a atof()
        m_refValueD = atof(m_refValueS.c_str());
    } else if (newType == 'I' && m_type == 'S') {
        // string to int, do a atoi()
        m_refValueD = atoi(m_refValueS.c_str());
    }
    m_type = newType;
}

