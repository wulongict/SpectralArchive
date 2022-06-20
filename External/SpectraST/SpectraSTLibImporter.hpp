#ifndef SPECTRASTLIBIMPORTER_HPP_
#define SPECTRASTLIBIMPORTER_HPP_

#include "SpectraSTLib.hpp"
#include "SpectraSTCreateParams.hpp"
#include "SpectraSTLibEntry.hpp"
#include "Peptide.hpp"
#include "Predicate.hpp"
#include <fstream>
#include <vector>
#include <string>


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

/* Class: SpectraSTLibImporter
 * 
 * Abstract base class that implements a library importer. 
 * 
 * To import a new file format, subclass from this class, implement the constructor and the readFromFile method, and modify
 * the static instantiation method createSpectraSTLibImporter to call your own constructor. 
 * 
 */


using namespace std;


class SpectraSTLibImporter {

public:

    SpectraSTLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual ~SpectraSTLibImporter();

    virtual void import() = 0; // abstract method - must be implemented by subclasses

    virtual void printStats();

    // "factory" method to create objects of subclasses
    static SpectraSTLibImporter *
    createSpectraSTLibImporter(vector<string> &impFileNames, SpectraSTLib *lib, SpectraSTCreateParams &params);

    virtual bool passAllFilters(SpectraSTLibEntry *entry);

    bool satisfyFilterCriteria(SpectraSTLibEntry *entry);

    bool isInProbTable(SpectraSTLibEntry *entry, bool updateProb = false);

    bool isInProteinList(SpectraSTLibEntry *entry);

    bool isAllDecoyProteins(SpectraSTLibEntry *entry);

    virtual string constructOutputFileName();

    string getOutputFileName() { return (m_outputFileName); }

    void printProteinList();

protected:

    vector<string> m_preamble;

    vector<string> m_impFileNames;

    SpectraSTCreateParams &m_params;

    string m_outputFileName;

    // m_lib - pointer to the library object. This is NOT the property of SpectraSTLibImporter!
    SpectraSTLib *m_lib;

    // m_filter - pointer to a Predicate object which stores what we are filtering for
    vector<Predicate *> *m_filters;
    string m_filterStr;
    char m_filterLogic;
    unsigned int m_count;

    // m_probTable - pointer to a hash of (peptide ion) => prob
    // the string is of the format ABC[330]DEFGHIK/2 (i.e. the format of m_name used in SpectraSTLibEntry)
    map<string, double> *m_probTable;

    // m_proteinList - pointer to a hash of (protein) => allowable # of entries remaining
    map<string, int> *m_proteinList;

    // helper methods
    void parseFilterCriteria();

    void readProbTable();

    void readProteinList();

    Peptide *createPeptide(string peptide, int charge, string modStr, string spectrum, string fileType);

    bool insertOneEntry(SpectraSTLibEntry *entry, string fileType);

    void setDeamidatedNXST(Peptide *pep);

};

#endif /*SPECTRASTLIBIMPORTER_HPP_*/
