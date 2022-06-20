#include "XMLWalker.hpp"
#include "FileUtils.hpp"
#include <iostream>


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

// Constructor for XML file input
// fields: a (attr => value) hash, the user should supply the attr's, and the parser will fill in the value's.
// recordAll: set to TRUE if you want multiple instances of the found attr to be returned as a concatenated string 
//            if FALSE, it will only return the first-seen instance.
// divider: a string to delineate multiple instances, only applies when recordAll is set to TRUE.
XMLWalker::XMLWalker(string xmlFileName, map<string, string> *fields, bool recordAll, string divider) :
        m_xmlFileName(xmlFileName),
        m_fields(fields),
        m_recordAll(recordAll),
        m_divider(divider),
        m_fin(NULL),
        m_good(true),
        m_curLine(""),
        m_stopLinePrefixes() {

    m_fin = new ifstream();
    if (!myFileOpen(*((ifstream *) m_fin), m_xmlFileName)) {
        m_good = false;
        return;
    }

    setContexts();
}

// Constructor for string input
// fields: a (attr => value) hash, the user should supply the attr's, and the parser will fill in the value's.
// recordAll: set to TRUE if you want multiple instances of the found attr to be returned as a concatenated string 
//            if FALSE, it will only return the first-seen instance.
// divider: a string to delineate multiple instances, only applies when recordAll is set to TRUE.
XMLWalker::XMLWalker(istream *bufss, map<string, string> *fields, bool recordAll, string divider) :
        m_xmlFileName(""),
        m_fields(fields),
        m_recordAll(recordAll),
        m_divider(divider),
        m_fin(bufss),
        m_good(true),
        m_curLine(""),
        m_stopLinePrefixes() {

    setContexts();
}

// Destructor
XMLWalker::~XMLWalker() {
    if (m_fin) {
        delete (m_fin);
    }
}

// good - returns TRUE if the XML file opens okay
bool XMLWalker::good() {
    return (m_good);
}

// current - returns the current line
string XMLWalker::current() {
    return (m_curLine);
}

// skip - jumps until it reaches a line prefixed with untilLinePrefix. No parsing is done in the process.
bool XMLWalker::skip(string untilLinePrefix) {

    if (!m_good) {
        return (false);
    }

    while (nextLine(*m_fin, m_curLine, untilLinePrefix));
    if (m_curLine == "_EOF_") {
        return (false);
    }
    return (true);

}

// setFields - set the (attr => value) hash.
// fields: a (attr => value) hash, the user should supply the attr's, and the parser will fill in the value's.
// recordAll: set to TRUE if you want multiple instances of the found attr to be returned as a concatenated string 
//            if FALSE, it will only return the first-seen instance.
// divider: a string to delineate multiple instances, only applies when recordAll is set to TRUE.
void XMLWalker::setFields(map<string, string> *fields, bool recordAll, string divider) {
    m_fields = fields;
    m_recordAll = recordAll;
    m_divider = divider;
    setContexts();
}

void XMLWalker::stopAt(vector<string> &stopLinePrefixes) {
    m_stopLinePrefixes = stopLinePrefixes;
}

// setContexts - parses out the (attr => value) and stores any contexts that we need to monitor
void XMLWalker::setContexts() {

    // m_contexts is a hash of (context => status). context is usually an element name.
    // status = 0 if we are not currently inside that element
    // status = 1 if we are currently inside the context, and this context will be done when we go to the next line
    //            (for those <foo/> type elements)
    // status = 2 if we are currently inside the context, and will remain inside until we hit the close tag

    m_contexts.clear();

    map<string, string>::iterator f;
    for (f = m_fields->begin(); f != m_fields->end(); f++) {
        string::size_type atPos = 0;
        if ((atPos = (*f).first.find('@', 0)) != string::npos) {
            // context specified
            string newContext = (*f).first.substr(atPos + 1);
            m_contexts[newContext] = 0;
        }
    }
}

// walk - the main parsing method. go find a line prefixed with startLinePrefix, start there, and return if
// we hit a line prefixed with endLinePrefix.
bool XMLWalker::walk(string startLinePrefix, string endLinePrefix) {

    if (!m_good) {
        return (false);
    }

    // clear the (attr => value) hash
    if (m_fields) {
        map<string, string>::iterator i;
        for (i = m_fields->begin(); i != m_fields->end(); i++) {
            (*i).second = "_NOT_FOUND_";
        }
    }

    // jumps until we hit startLinePrefix
    if (m_curLine.compare(0, startLinePrefix.length(), startLinePrefix) != 0) {
        while (nextLineCropped(*m_fin, m_curLine, startLinePrefix)) {

            // check if we hit a stopAt line
            for (vector<string>::iterator slpfx = m_stopLinePrefixes.begin();
                 slpfx != m_stopLinePrefixes.end(); slpfx++) {
                if (m_curLine.substr(0, slpfx->length()) == *slpfx) {
                    return (false);
                }
            }
        }
        if (m_curLine == "_EOF_") {
            return (false);
        }
    }

    do {

        // check if we hit a stopAt line
        for (vector<string>::iterator slpfx = m_stopLinePrefixes.begin(); slpfx != m_stopLinePrefixes.end(); slpfx++) {
            if (crop(m_curLine).substr(0, slpfx->length()) == *slpfx) {
                return (false);
            }
        }

        if (!m_fields) {
            continue;
        }

        // monitor the contexts
        map<string, int>::iterator cx;
        for (cx = m_contexts.begin(); cx != m_contexts.end(); cx++) {

            string::size_type openContext = m_curLine.find("<" + cx->first, 0);
            string::size_type closeContext = m_curLine.find("</" + cx->first + ">", 0);

            if (openContext != string::npos && closeContext != string::npos) {
                cx->second = 1; // one line element, like <foo>....</foo>
                continue;
            }

            if (closeContext != string::npos) {
                cx->second = 0;   // closing this context since we just saw a close tag, </foo>
                continue;
            }
            if (openContext == string::npos) {
                continue;
            }

            string::size_type oneLineClose = m_curLine.find("/>", openContext + cx->first.size() + 1);
            if (oneLineClose != string::npos) {
                cx->second = 1;   // this is a one-line element, like <foo a1=v1 a2=v2/>
                continue;
            }

            string::size_type closeBracket = m_curLine.find(">", openContext + cx->first.size() + 1);
            if (closeBracket != string::npos) {
                cx->second = 2;  // this is a multi-line, possibly. open the context and wait for the close context.
            } else {
                // an open element that drags more than one line... bad
                // ignore this element together...
            }
        }

        // go through each attr in m_fields and tries to find the values.
        map<string, string>::iterator f;
        for (f = m_fields->begin(); f != m_fields->end(); f++) {

            string field = (*f).first;
            string::size_type atPos = (*f).first.find('@', 0);
            if (atPos != string::npos) {
                // context specified
                string context = (*f).first.substr(atPos + 1);
                if (m_contexts[context] == 0) {
                    // not in context! forget about this field
                    continue;
                } else {
                    // in context!
                    field = (*f).first.substr(0, atPos);
                    if (field.empty()) {
                        // entire context wanted
                        if (f->second == "_NOT_FOUND_") {
                            f->second = m_curLine + "\n";
                        } else {
                            f->second += (m_curLine + "\n");
                        }
                    }
                }
            }

            // we're not doing recordAll, already found this attr before, no need to find it again.
            if (!m_recordAll && (*f).second != "_NOT_FOUND_" && !field.empty()) {
                continue;
            }

            // we allow two types of attr=>value associations
            // 1. as an attribute within an element, e.g. <foo attr="value">
            string attrStr = ' ' + field + '=';
            string::size_type attrStrLen = attrStr.length();
            // 2. as a (name, value) pair of attributes, e.g. <foo name="attr" value="value">
            string paramStr = "\"" + field + "\"";
            string::size_type paramStrLen = paramStr.length();


            string::size_type dummy = 0;
            string::size_type start = 0;

            if (!field.empty() && (start = m_curLine.find(attrStr, 0)) != string::npos) {
                //this is an attribute, what's following the '=' is the value we want
                string value = nextToken(m_curLine, start + attrStrLen + 1, dummy, "\t\r\n\"");

                if (!m_recordAll) {
                    (*f).second = value;
                } else if ((*f).second == "_NOT_FOUND_") {
                    (*f).second = value;
                } else {
                    (*f).second += m_divider + value;
                }

            } else if (!field.empty() && (start = m_curLine.find(paramStr, 0)) != string::npos) {
                string pre("");
                if (start >= 5)
                    pre.assign(m_curLine, start - 5, 5);
                string post("");
                string::size_type quotePos = 0;
                post = nextToken(m_curLine, start + paramStrLen, quotePos, "\"\t\r\n");
                if (pre == "name=" && post == "value=") {

                    // this is a parameter, what's following "value=" is the value we want
                    string value = nextToken(m_curLine, quotePos + 1, dummy, "\"\t\r\n");
                    if (!m_recordAll) {
                        (*f).second = value;
                    } else if ((*f).second == "_NOT_FOUND_") {
                        (*f).second = value;
                    } else {
                        (*f).second += m_divider + value;
                    }
                }
            }
        }



        // close all one-line contexts
        for (cx = m_contexts.begin(); cx != m_contexts.end(); cx++) {
            if ((*cx).second == 1) {
                (*cx).second = 0;
            }
        }


    } while (nextLineCropped(*m_fin, m_curLine, endLinePrefix));

    if (m_curLine == "_EOF_") {
        return (false);
    } else {
        return (true);
    }

}



