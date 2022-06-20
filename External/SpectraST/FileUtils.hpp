#ifndef FILEUTILS_HPP_
#define FILEUTILS_HPP_

#ifdef __MINGW__
#define WINDOWS_CYGWIN
#endif

#ifdef STANDALONE_LINUX

#include "SpectraST_util.h" // building standalone
#include "SpectraST_constants.h" // building standalone

#ifndef GNUPLOT_BINARY
#ifdef WINDOWS_CYGWIN
#define GNUPLOT_BINARY "wgnuplot.exe"
#else
#define GNUPLOT_BINARY "gnuplot"
#endif
#endif

#else
#include "common/util.h"
#include "common/constants.h" // the TPP defaults
#endif

#include <string>
#include <fstream>
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




using namespace std;

// a struct for file names
typedef struct _fileName {
    string path;
    string name;
    string ext;
} FileName;

// FILE HANDLING METHODS

// opening files, use the right version for input or output
bool myFileOpen(ifstream &fin, string fullFileName, bool binary = false);

bool myFileOpen(ofstream &fout, string fullFileName, bool binary = false);

// deleting files
void removeFile(string fileName);

void removeDir(string dir);

// making new directory
void makeDir(string dir);

// using tar to create/extract from tgz's
void myTgz(string tgzFullFileName, string filesToOpen, bool compress = false);

// calling gnuplot to plot
void myGnuplot(string plotFileName, string gnuBinary = "");

// getting the default cgi-bin path
string getCGIPath();



// FILE NAME HANDLING METHODS

// parse full file names to get the path, extension, and the file name
void parseFileName(string fullFileName, FileName &fn);

string::size_type getPath(string fullFileName, string &path);

string::size_type getExtension(string fullFileName, string &ext);


// LINE PARSERS

bool nextLine(istream &fin, string &line);

bool nextLine(istream &fin, string &line, string flagLinePrefix);

bool nextLine(istream &fin, string &line, string flagLinePrefix, string skipoverLinePrefix);

bool nextLineCropped(istream &fin, string &line, string flagLinePrefix);

bool nextLineCropped(istream &fin, string &line, string flagLinePrefix, string skipoverLinePrefix);

string crop(string s);

// TOKENIZER

string nextToken(string s, string::size_type from, string::size_type &tokenEnd, const char *delim = " \t\r\n",
                 const char *skipover = " \t\r\n");

vector<string> *split(string s, const char *delim = " \t\r\n", const char *skipover = " \t\r\n");

#endif /*FILEUTILS_HPP_*/
