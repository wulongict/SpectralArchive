#include "FileUtils.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "cramp.hpp"
#endif

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

/*

Library       : FileUtils
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

/* FileUtils
 * 
 * Some utility functions for managing files. THIS IS MEANT TO BE AN INTERFACE OF SPECTRAST AND THE FILE SYSTEM.
 * ALL FILE-RELATED OPERATIONS OF SPECTRAST ARE DONE THROUGH THIS CLASS, SO THIS IS WHERE ALL FILE-SYSTEM SPECIFIC
 * IMPLEMENTATIONS NEED TO BE. THE REST OF SPECTRAST IS IGNORANT OF THE FILE SYSTEM.
 *    
 * THIS IS ALSO WHERE SPECTRAST GETS ALL THE TPP-WIDE DEFAULTS (through constants.h). NOWHERE ELSE IN SPECTRAST
 * NEEDS TO KNOW ABOUT THE REST OF TPP.
 */

// FILE HANDLING METHODS

// myFileOpen - opens a file for input, checks whether successful
bool myFileOpen(ifstream &fin, string fullFileName, bool binary) {

    if (binary) {
        fin.open(fullFileName.c_str(), ios::binary | ios::in);
    } else {
        fin.open(fullFileName.c_str());
    }
    if (!fin.good()) {
        return (false);
    } else {
        return (true);
    }
    return (true);
}

// myFileOpen - opens a file for output, checks whether successful
bool myFileOpen(ofstream &fout, string fullFileName, bool binary) {

    if (binary) {
        fout.open(fullFileName.c_str(), ios::binary | ios::out);
    } else {
        fout.open(fullFileName.c_str());
    }
    if (!fout.good()) {
        return (false);
    } else {
        return (true);
    }
    return (true);
}

// removeFile - deletes a file from the file system -- this is UNIX, change for Windows 
void removeFile(string fullFileName) {
    string rmCommand("rm -f ");
    rmCommand += fullFileName;
    int ret = system(rmCommand.c_str());
}

// removeDir - deletes an entire directory from the file system -- this is UNIX, change for Windows 
void removeDir(string dir) {
    string rmCommand("rm -rf ");
    rmCommand += dir;
    int ret = system(rmCommand.c_str());
}

// makeDir - creates a directory -- this is UNIX, change for Windows
void makeDir(string dir) {
    string mkdirCommand("mkdir -p ");
    mkdirCommand += dir;
    int ret = system(mkdirCommand.c_str());
}

// myTgz - uses tar to compress files to tgz's or extract files from tgz's
// if compressing, 'filesToOpen' is required. if extracting, it can be omitted to extract everything in that
// tgz. Note also that this function will "cd" to that directory first so that the filenames in the tgz will
// be relative (easier to move the tgz around).
void myTgz(string tgzFullFileName, string filesToOpen, bool compress) {

    FileName fn;
    parseFileName(tgzFullFileName, fn);

    string cdCommand("");

    // if path is given, cd to the directory first
    if (!fn.path.empty()) {
        cdCommand += "cd " + fn.path + ";";
    }

    string tgzCommand(cdCommand + " tar ");

    if (compress) {
        tgzCommand += "czf ";
    } else {
        tgzCommand += "xzf ";
    }
    tgzCommand += fn.name + fn.ext + ' ';

    if (compress && filesToOpen.empty()) {
        // don't know what to compress. just abort
        cerr << "No specification of what files to compress into a .tgz. No action taken." << endl;
        return;
    }

    tgzCommand += filesToOpen;
#ifdef STANDALONE_LINUX
    int ret = system(tgzCommand.c_str());
#else
    verified_system(tgzCommand.c_str()); // like system(), but handles multiple commmands in win32
#endif
}

// myGnuplot - calls gnuplot by system call - can specify the full path to the binary, or use the TPP default in constants.h
void myGnuplot(string plotFileName, string gnuplotBinary) {
    stringstream plotCommand;
    if (gnuplotBinary.empty()) {
        plotCommand << GNUPLOT_BINARY << ' ' << plotFileName;
    } else {
        plotCommand << gnuplotBinary << ' ' << plotFileName;
    }
#ifdef STANDALONE_LINUX
    int ret = system(plotCommand.str().c_str());
#else
    verified_system(plotCommand.str().c_str());
#endif
}

// getCGIPath - just retrieves the TPP default cgi-bin path in constants.h
string getCGIPath() {

#ifdef STANDALONE_LINUX
    return (CGI_BIN);
#else
    string tpp(getCGIFullBin());
    return (tpp);
#endif
}


// FILE NAME HANDLING METHODS

// parseFileName - Takes a full file name, i.e. <path>/<fileName>.<extension> and parses out 
// the pieces.
void parseFileName(string fullFileName, FileName &fn) {

    string::size_type dotPos = getExtension(fullFileName, fn.ext);
    string::size_type slashPos = getPath(fullFileName, fn.path);
    if (slashPos == string::npos) {
        fn.name = fullFileName.substr(0, dotPos);
    } else {
        fn.name = fullFileName.substr(slashPos + 1, dotPos - slashPos - 1);
    }
    //cout << "DEBUG:FILEUTILS:PARSEFILENAME>>" << fullFileName << endl;
    //cout << "DEBUG:FILEUTILS:PARSEFILENAME>>" << fn.name << endl;
}

// getPath - just gets the path from the fullFileName. return the position of the slash separating
// the path and the file name.
string::size_type getPath(string fullFileName, string &path) {
    // find the last slash '/' in the string; the path is assumed to the everything to the left
    // of this slash. The path WILL CONTAIN THIS SLASH!
    string::size_type len = fullFileName.length();
#ifdef STANDALONE_LINUX
    string::size_type pos = fullFileName.rfind('/', len);
#else
    string::size_type pos = (string::size_type)(findRightmostPathSeperator(fullFileName)); // return pos to rightmost / or \, or npos.
#endif

    if (pos == string::npos) {
        // no slash whatsoever, hence no path
        path = "";
        return (string::npos);
    } else {
        path = fullFileName.substr(0, pos + 1);
        return (pos);
    }
    return (string::npos);
}

// getExtension - just gets the extension from the fullFileName. return the position of the dot separating
// the extension and the file name.
string::size_type getExtension(string fullFileName, string &ext) {
    // find the last dot '.' in the string; the extension is assumed to the everything to the right
    // of this dot. The ext WILL CONTAIN THIS DOT!

    string::size_type len = fullFileName.length();
#ifdef STANDALONE_LINUX
    string::size_type slashPos = fullFileName.rfind('/', len);
#else
    string::size_type slashPos = (string::size_type)(findRightmostPathSeperator(fullFileName)); // return pos to rightmost / or \, or npos.
#endif

    // check special case of ramp-compatible formats (esp. to deal with .gz files)
    char *rampExt = rampValidFileType(fullFileName.c_str());
    if (rampExt) {
        ext = rampExt;
        return (len - ext.length());
    }

    string::size_type pos = fullFileName.rfind('.', len);

    if (slashPos != string::npos && pos != string::npos && slashPos > pos) {
        // dot before a slash, this means the dot is part of the path, not the extension of the file
        // hence, no extension.
        ext = "";
        return (len);
    }

    if (pos == string::npos) {
        // no dot whatsoever, hence no extension
        ext = "";
        return (len);
    } else {
        // check special case of .pep.xml
        if (pos >= 4 && fullFileName.substr(pos - 4, 4) == ".pep" && fullFileName.substr(pos) == ".xml") {
            pos -= 4;
        }
        ext = fullFileName.substr(pos);
        return (pos);
    }
    return (len);
}

// LINE PARSERS

// nextLine - a smarter getline that will remove the trailing carriage return '\r'. This is
// necessary since many text files have lines terminated by \r\n rather than \n, and the C++
// library function getline() stubbornly keeps the '\r' as part of the string.
bool nextLine(istream &fin, string &line) {

    if (fin.eof()) {
        line = "_EOF_";
        return (false);
    }

    getline(fin, line);
    if (line.empty()) return (true);

    string::size_type last = line.length() - 1;
    if (line[last] == '\r') {
        line.erase(last);
    }
    return (true);
}


// nextLine - smart line-by-line file parser. puts the next line of an input file stream in 'line',
// Flags (by returning false) if it gets a line starting with 'flagLinePrefix'. 
// 'flagLinePrefix' being empty means flag at a blank line.
// line will contain the string "_EOF_" if it reaches the end of the file. Will also return false in those
// cases. 
bool nextLine(istream &fin, string &line, string flagLinePrefix) {

    if (!nextLine(fin, line)) {
        return (false);
    }

    if ((flagLinePrefix.empty() && line.empty()) ||
        (!flagLinePrefix.empty() && line.compare(0, flagLinePrefix.length(), flagLinePrefix) == 0)) {
        return (false);
    } else {
        return (true);
    }
    return (false);

}

// nextLineCropped - same as nextLine except that the line is cropped
bool nextLineCropped(istream &fin, string &line, string flagLinePrefix) {

    if (!nextLine(fin, line)) {
        return (false);
    }

    line = crop(line);

    if ((flagLinePrefix.empty() && line.empty()) ||
        (!flagLinePrefix.empty() && line.compare(0, flagLinePrefix.length(), flagLinePrefix) == 0)) {
        return (false);
    } else {
        return (true);
    }
    return (false);
}

// nextLine - smart line-by-line file parser. puts the next line of an input file stream in 'line',
// skipping over any line starting with 'skipoverLinePrefix'. Also 
// flags (by returning false) if it gets a line starting with 'flagLinePrefix'. 
// 'skipoverLinePrefix' being empty means skip over all blank lines
// 'flagLinePrefix' being empty means flag at a blank line 
// 'line' will contain the string "_EOF_" if it reaches the end of the file (after skipping over
// every line it's supposed to skip over). Will also return false in those cases. 
bool nextLine(istream &fin, string &line, string flagLinePrefix, string skipoverLinePrefix) {

    line = "";
    while (nextLine(fin, line)) {
        if ((skipoverLinePrefix.empty() && line.empty()) ||
            (!skipoverLinePrefix.empty() && line.compare(0, skipoverLinePrefix.length(), skipoverLinePrefix) == 0)) {

            // skip over

        } else {
            break;
        }
    }
    if (line == "_EOF_") {
        return (false);
    }
    if ((flagLinePrefix.empty() && line.empty()) ||
        (!flagLinePrefix.empty() && line.compare(0, flagLinePrefix.length(), flagLinePrefix) == 0)) {
        return (false);
    } else {
        return (true);
    }
    return (false);
}


// nextLineCropped - same as nextLine except that the line is cropped
bool nextLineCropped(istream &fin, string &line, string flagLinePrefix, string skipoverLinePrefix) {

    line = "";
    while (nextLine(fin, line)) {
        line = crop(line);
        if ((skipoverLinePrefix.empty() && line.empty()) ||
            (!skipoverLinePrefix.empty() && line.compare(0, skipoverLinePrefix.length(), skipoverLinePrefix) == 0)) {

            // skip over

        } else {
            break;
        }
    }
    if (line == "_EOF_") {
        return (false);
    }
    if ((flagLinePrefix.empty() && line.empty()) ||
        (!flagLinePrefix.empty() && line.compare(0, flagLinePrefix.length(), flagLinePrefix) == 0)) {
        return (false);
    } else {
        return (true);
    }
    return (false);


}

// crop - crops the leading and trailing blank spaces from a string
string crop(string s) {

    string::size_type firstPos = s.find_first_not_of(" \t\r\n");
    string::size_type lastPos = s.find_last_not_of(" \t\r\n");
    if (firstPos == string::npos || lastPos == string::npos) {
        return ("");
    }

    return (s.substr(firstPos, lastPos - firstPos + 1));

}


// TOKENIZER

// nextToken - returns the next token from the string. Here's how it works, it starts scanning from the
// position 'from', skips over all characters in the string 'skipover', then reads in the token until
// it reaches any character in the string 'delim'. 'tokenEnd' will point to the position AFTER the end
// of the token (i.e. where the delimiter is).  
// 
// Typical usage: consecutive calls of nextToken(s, pos, pos, delim, skipover) will return the tokens in
// the string one by one, starting at position 'pos'.
string
nextToken(string s, string::size_type from, string::size_type &tokenEnd, const char *delim, const char *skipover) {

    if (from == string::npos || from >= s.length()) {
        // already got nothing
        tokenEnd = s.length();
        return ("");
    }

    // skips over all characters in the string 'skipover'
    string::size_type tokenStart = s.find_first_not_of(skipover, from);

    if (tokenStart == string::npos) {
        // whoa, everything from 'from' onwards are characters to be skipped over.
        // That is, no token exists. For this special case, tokenEnd is set to the end of the string.
        tokenEnd = s.length();
        return ("");
    }

    // reads until it reaches a delimiter
    tokenEnd = s.find_first_of(delim, tokenStart);
    if (tokenEnd == string::npos) {
        // token is all the way to the end of s
        tokenEnd = s.length();
        return (s.substr(tokenStart));
    }

    return (s.substr(tokenStart, tokenEnd - tokenStart));

}


vector<string> *split(string s, const char *delim, const char *skipover) {

    vector<string> *vs = new vector<string>();
    string::size_type pos = 0;

    while (pos < s.length()) {
        string token = nextToken(s, pos, pos, delim, skipover);
        vs->push_back(token);
    }
    return (vs);
}
