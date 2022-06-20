#ifndef XMLWALKER_HPP_
#define XMLWALKER_HPP_

#include <vector>
#include <map>
#include <string>
#include <fstream>


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

/* Class: XMLWalker
 * 
 * A simplistic parser of XML files -- only works with pepXML files generated at ISB.
 * It needs the whitespaces to be in the right places -- which is most commonly true, but not guaranteed.
 * 
 * Assumptions made (there could be more):
 * - Nothing can be outside angle brackets (the pepXML schema also has this requirement anyway)
 * - No line (\n seperated) can have more than one element (<foo attr="value"><bar/></foo> is not okay on a line)
 * - If an element spans multiple lines, it must have open and close tags; <foo/> type elements cannot span more than one line.
 * - Important "mile post" open and end tags are at the start of the line.
 * 
 * Can also be used to parse XML segments in a string buffer. This can be used for precise parsing of nested structures, in which
 * an inner element is first parsed from an outer element, followed by parsing of the inner element as a string, etc. 
 *
 * The basic idea is the following. The user supplies a list of attributes in a hash of (attr => value), that he/she want to extract 
 * out of the file. 
 * Then walk() is called repeatedly, which will parse the XML files and detect any desired attr along the way.
 * walk() will return TRUE as it hits "mile posts" -- which are also specified by the user. 
 * The user can then examine the hash to get the values (if any found).
 * walk() will return FALSE if it can't find any "mile post" left, at which point the parsing is done.
 * The user can specify attributes in context, for instance, "name@employer" asks the parser to only find any "name" fields within
 * an employer element, not just any "name" fields. The user can also grab the whole context by supplying "@employer", which will
 * cause the parser to return everything within an employer element as a string.
 *
 * See SpectraSTPepXMLLibImporter.cpp for how to use it.
 */

using namespace std;

class XMLWalker {
public:

    // For XML file
    XMLWalker(string xmlFileName, map<string, string> *fields = NULL, bool recordAll = false, string divider = "|");

    // For XML segment in a string buffer
    XMLWalker(istream *bufss, map<string, string> *fields = NULL, bool recordAll = false, string divider = "|");

    // Destructor
    virtual ~XMLWalker();

    void setFields(map<string, string> *fields, bool recordAll = false, string divider = "|");

    // signals good opening of file
    bool good();

    bool skip(string untilLinePrefix);

    bool walk(string startLinePrefix, string endLinePrefix);

    string current();

    void stopAt(vector<string> &linePrefixes);

private:

    // the XML file name
    string m_xmlFileName;

    // pairs of (attr, value). The user supplies the desired attrs, the parser will try to fill in the values as it walks.
    map<string, string> *m_fields;

    // stores the contexts that the user wants to enforce for each field
    map<string, int> m_contexts;
    string m_curLine;
    string m_buffer;

    bool m_recordAll;
    string m_divider;

    istream *m_fin;
    bool m_good;

    void setContexts();

    vector<string> m_stopLinePrefixes;

};

#endif /*XMLWALKER_HPP_*/

