#ifndef PREDICATE_HPP
#define PREDICATE_HPP

#include <string>

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

class Predicate {

public:
    Predicate(string refName, string oper, string refValue, bool defaultBool = true);

    Predicate(string refName, string oper, int refValue, bool defaultBool = true);

    Predicate(string refName, string oper, double refValue, bool defaultBool = true);

    ~Predicate();

    string getRefName() { return (m_refName); }

    char getType() { return (m_type); }


    // compare to reference value
    bool compRef(string value);

    bool compRef(int value);

    bool compRef(double value);

    // compare two values
    bool comp(string value1, string value2);

    bool comp(int value1, int value2);

    bool comp(double value1, double value2);

    // change the type (int, double, string) to compare
    void changeType(char newType);

private:
    string m_oper; // operand, can be ==, !=, >, >=, <, <=, =~, !~

    char m_type;

    string m_refName;
    string m_refValueS;
    int m_refValueI;
    double m_refValueD;

    bool m_default; // default return value for failed comparison

};

#endif
