#include "SpectraSTLibImporter.hpp"
#include "SpectraSTMspLibImporter.hpp"
#include "SpectraSTXHunterLibImporter.hpp"
#include "SpectraSTPepXMLLibImporter.hpp"
#include "SpectraSTSpLibImporter.hpp"
#include "SpectraSTMs2LibImporter.hpp"
#include "SpectraSTTsvLibImporter.hpp"
#include "SpectraSTMzXMLLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Predicate.hpp"

#include <sstream>
#include <iostream>
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

/* Class: SpectraSTLibImporter
 * 
 * Abstract base class that implements a library importer. 
 * 
 * To import a new file format, subclass from this class, implement the constructor and the readFromFile method, and modify
 * the static instantiation method createSpectraSTLibImporter to call your own constructor. 
 * 
 */

extern SpectraSTLog *g_log;

// Constructor. If you subclass, REMEMBER to call this constructor of the parent class in your constructor!
SpectraSTLibImporter::SpectraSTLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                           SpectraSTCreateParams &params) :
        m_impFileNames(impFileNames),
        m_lib(lib),
        m_params(params),
        m_count(0),
        m_preamble(),
        m_outputFileName(""),
        m_filters(NULL),
        m_filterStr(""),
        m_filterLogic('\0'),
        m_probTable(NULL),
        m_proteinList(NULL) {

    if (lib) {

        // create the output library file name
        if (params.outputFileName.empty()) {
            m_outputFileName = constructOutputFileName();
        } else {
            m_outputFileName = params.outputFileName + ".splib";
        }

        // create the Predicate object if filterCriteria is specified
        if (!m_params.filterCriteria.empty()) {
            parseFilterCriteria();
        }

        // read the probability table if it is specified
        if (!m_params.useProbTable.empty()) {
            readProbTable();
        }

        // read the protein list if it is specified
        if (!m_params.useProteinList.empty()) {
            readProteinList();
        }

    }
}

// parseFilterCriteria - parses the filter criterion string.
// The criterion string is of the form <attr> <op> <value> [<logic> <attr> <op> <value]... where <attr> is any of the
// field of the entry or any of the fields in the comments. <logic> can be & (and) or | (or), but only one
// can be used throughout in the criterion. (If there are both & and |, & will be assumed throughout.) 
// <op> can be != (not equal), == (equal), <= (less than or equal), >= (greater than or equal)
// > (greater than), < (smaller than), =~ (contains as a substring) and !~ (does not contain)
// <value> can be numerical or string. It will figure out what it should be (hopefully).
void SpectraSTLibImporter::parseFilterCriteria() {

    m_filters = new vector<Predicate *>;

    string::size_type critStartPos = 0;
    string::size_type logicPos = 0;

    while (critStartPos < m_params.filterCriteria.length()) {

        logicPos = m_params.filterCriteria.find_first_of("&|", critStartPos);
        if (logicPos != string::npos) {
            if (m_filterLogic != '\0') {
                if (m_filterLogic != m_params.filterCriteria[logicPos]) {
                    g_log->error("CREATE",
                                 "FILTER: Mixing of AND and OR logics in filter criteria. AND logic is assumed throughout.");
                    m_filterLogic = '&';
                }
            } else {
                m_filterLogic = m_params.filterCriteria[logicPos];
            }
        } else {
            logicPos = m_params.filterCriteria.length();
            if (m_filterLogic == '\0') m_filterLogic = '&';
        }

        string criterion = m_params.filterCriteria.substr(critStartPos, logicPos - critStartPos);
        critStartPos = logicPos + 1;

        string attr("");
        string value("");
        string opStr("");
        string::size_type opPos = 0;
        if ((opPos = criterion.find("!=", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = "!=";
        } else if ((opPos = criterion.find(">=", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = ">=";
        } else if ((opPos = criterion.find("<=", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = "<=";
        } else if ((opPos = criterion.find("==", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = "==";
        } else if ((opPos = criterion.find(">", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 1));
            opStr = ">";
        } else if ((opPos = criterion.find("<", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 1));
            opStr = "<";
        } else if ((opPos = criterion.find("=~", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = "=~";
        } else if ((opPos = criterion.find("!~", 0)) != string::npos) {
            attr = crop(criterion.substr(0, opPos));
            value = crop(criterion.substr(opPos + 2));
            opStr = "!~";
        }

        if (opStr.empty() || attr.empty() || value.empty()) {
            continue;
        }

        // Some common ones; for these we know what the type (int, double or string) the value should be
        if (attr == "Name") {
            m_filters->push_back(new Predicate(attr, opStr, value));
        } else if (attr == "MW") {
            m_filters->push_back(new Predicate(attr, opStr, atof(value.c_str())));
        } else if (attr == "PrecursorMZ") {
            m_filters->push_back(new Predicate(attr, opStr, atof(value.c_str())));
        } else if (attr == "LibID") {
            m_filters->push_back(new Predicate(attr, opStr, atoi(value.c_str())));
        } else if (attr == "Charge") {
            m_filters->push_back(new Predicate(attr, opStr, atoi(value.c_str())));
        } else if (attr == "Mods") {
            m_filters->push_back(new Predicate(attr, opStr, value));
        } else if (attr == "Status") {
            m_filters->push_back(new Predicate(attr, opStr, value));
        } else if (attr == "FullName") {
            m_filters->push_back(new Predicate(attr, opStr, value));
        } else if (attr == "NumPeaks") {
            m_filters->push_back(new Predicate(attr, opStr, atoi(value.c_str())));
        } else if (attr == "FragType") {
            m_filters->push_back(new Predicate(attr, opStr, value));
        } else {
            // Not the common ones, the type is unknown. Will guess based on the type of the value supplied
            if (value.find_first_not_of("0123456789.") == string::npos) { // assume it's a double
                m_filters->push_back(new Predicate(attr, opStr, atof(value.c_str())));
            } else { // assume it's a string
                m_filters->push_back(new Predicate(attr, opStr, value));
            }
        }

    }

    // put together a string to describe what we are filtering for
    stringstream fss;
    fss << "FILTER for " << m_params.filterCriteria;
    g_log->log("CREATE", "Apply " + fss.str());

    if (!m_filterStr.empty()) {
        m_filterStr += "; ";
    }
    m_filterStr += fss.str();


}

// destructor
SpectraSTLibImporter::~SpectraSTLibImporter() {
    if (m_filters) {
        for (vector<Predicate *>::iterator i = m_filters->begin(); i != m_filters->end(); i++) {
            delete (*i);
        }
        delete m_filters;
    }
    if (m_probTable) {
        delete m_probTable;
    }
    if (m_proteinList) {
        delete m_proteinList;
    }
}

// printStats - nothing here.
void SpectraSTLibImporter::printStats() {}

// passAllFilters - calls each of the filtering methods
bool SpectraSTLibImporter::passAllFilters(SpectraSTLibEntry *entry) {

    return (isInProbTable(entry, true) &&
            isInProteinList(entry) &&
            satisfyFilterCriteria(entry) &&
            !isAllDecoyProteins(entry) &&
            (!(entry->getPeptidePtr()) || entry->getPeptidePtr()->isGood())
    );
}

// satisfyFilterCriteria - checks if the entry satisfies the criteria given by m_params.filterCriteria. See 
// the method parseFilterCriteria for more information. 
bool SpectraSTLibImporter::satisfyFilterCriteria(SpectraSTLibEntry *entry) {

    if (!m_filters) {
        return (true);
    }

    for (vector<Predicate *>::iterator pr = m_filters->begin(); pr != m_filters->end(); pr++) {

        Predicate *filter = (*pr);
        string attr = filter->getRefName();
        string commentValue("");

        if (m_filterLogic == '&') {

            if (attr == "Name") {
                if (!filter->compRef(entry->getName())) return (false);
            } else if (attr == "MW") {
                if (!filter->compRef(entry->getMw())) return (false);
            } else if (attr == "PrecursorMZ") {
                if (!filter->compRef(entry->getPrecursorMz())) return (false);
            } else if (attr == "LibID") {
                if (!filter->compRef((int) (entry->getLibId()))) return (false);
            } else if (attr == "Charge") {
                if (!filter->compRef((int) (entry->getCharge()))) return (false);
            } else if (attr == "Mods") {
                if (!filter->compRef(entry->getMods())) return (false);
            } else if (attr == "Status") {
                if (!filter->compRef(entry->getStatus())) return (false);
            } else if (attr == "FullName") {
                if (!filter->compRef(entry->getFullName())) return (false);
            } else if (attr == "NumPeaks") {
                if (!filter->compRef((int) (entry->getPeakList()->getNumPeaks()))) return (false);
            } else if (attr == "FragType") {
                if (!filter->compRef(entry->getFragType())) return (false);
            } else if (entry->getOneComment(attr, commentValue)) {
                string::size_type pos = 0;

                if (filter->getType() == 'D') {
                    if (!filter->compRef(atof(nextToken(commentValue, pos, pos, "/^\t\r\n").c_str()))) return (false);
                } else if (filter->getType() == 'I') {
                    if (!filter->compRef(atoi(nextToken(commentValue, pos, pos, "/^\t\r\n").c_str()))) return (false);
                } else {
                    if (!filter->compRef(commentValue)) return (false);
                }

            } else {
                // can't find it in comments, don't do anything - wait for the other predicates to evaluate
            }

        } else if (m_filterLogic == '|') {

            if (attr == "Name") {
                if (filter->compRef(entry->getName())) return (true);
            } else if (attr == "MW") {
                if (filter->compRef(entry->getMw())) return (true);
            } else if (attr == "PrecursorMZ") {
                if (filter->compRef(entry->getPrecursorMz())) return (true);
            } else if (attr == "LibID") {
                if (filter->compRef((int) (entry->getLibId()))) return (true);
            } else if (attr == "Charge") {
                if (filter->compRef((int) (entry->getCharge()))) return (true);
            } else if (attr == "Mods") {
                if (filter->compRef(entry->getMods())) return (true);
            } else if (attr == "Status") {
                if (filter->compRef(entry->getStatus())) return (true);
            } else if (attr == "FullName") {
                if (filter->compRef(entry->getFullName())) return (true);
            } else if (attr == "NumPeaks") {
                if (filter->compRef((int) (entry->getPeakList()->getNumPeaks()))) return (true);
            } else if (attr == "FragType") {
                if (!filter->compRef(entry->getFragType())) return (true);
            } else if (entry->getOneComment(attr, commentValue)) {
                string::size_type pos = 0;

                if (filter->getType() == 'D') {
                    if (filter->compRef(atof(nextToken(commentValue, pos, pos, "/^\t\r\n").c_str()))) return (true);
                } else if (filter->getType() == 'I') {
                    if (filter->compRef(atoi(nextToken(commentValue, pos, pos, "/^\t\r\n").c_str()))) return (true);
                } else {
                    if (filter->compRef(commentValue)) return (true);
                }

            } else {
                // can't find it in comments
            }
        }
    }

    if (m_filterLogic == '&') {
        return (true);
    } else if (m_filterLogic == '|') {
        return (false);
    } else {
        return (true);
    }

    return (true);
}

// readProbTable - reads a probability table from file. A prob table has lines of the format <peptide ion> <white space> <new prob>, 
// where <peptide ion> is the entry name in "interact" format: e.g. AC[160]DEFGHIK/2, and <new prob> is optional. Only library entries
// whose peptide ion (= entry name) matches any in the prob table will be retained. Optionally, if <new prob> is specified, the library
// entry's Prob= field in the comments will be updated to that alongside the peptide ion in the prob table.
void SpectraSTLibImporter::readProbTable() {

    ifstream fin;
    if (!myFileOpen(fin, m_params.useProbTable)) {
        g_log->error("CREATE", "Cannot read probability table. No probability adjustment/filtering will be performed.");
        return;
    }

    m_probTable = new map<string, double>;

    string line("");
    while (nextLine(fin, line, "_EOF_", "")) {
        string::size_type pos = 0;
        string pep = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
        if (pep.empty() || pep[0] == '#') {
            continue;
        }
        string probStr = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
        double prob = -1.0;
        if (!probStr.empty()) {
            prob = atof(probStr.c_str());
        }

        map<string, double>::iterator found = m_probTable->find(pep);
        if (found != m_probTable->end()) {
            if (prob > found->second) found->second = prob;
        } else {
            (*m_probTable)[pep] = prob;
        }
    }

    // put together a string to describe filtering by prob table
    stringstream fss;
    fss << "FILTER for entries in probability table \"" << m_params.useProbTable << "\"";
    g_log->log("CREATE", "Apply " + fss.str());

    if (!m_filterStr.empty()) {
        m_filterStr += "; ";
    }
    m_filterStr += fss.str();


}

// isInProbTable - returns true if the entry's peptide ion is in the prob table. Optionally, the peptide ion's associated
// probability can be updated to the one listed alongside the peptide ion in the prob table. This is a way to modify
// the probabilities after the initial import from PepXML, e.g. due to better validation scheme.
bool SpectraSTLibImporter::isInProbTable(SpectraSTLibEntry *entry, bool updateProb) {

    if (!m_probTable) {
        // no prob table, logic here is always return true, so that the entry won't be filtered
        return (true);
    }


    map<string, double>::iterator found = m_probTable->find(entry->getName());

    if (found == m_probTable->end()) {
        // also consider case when the prob table contains peptide only, without the charge
        string seqWithoutCharge("");

        if (entry->getPeptidePtr()) {
            seqWithoutCharge = entry->getPeptidePtr()->interactStyle();
        } else {
            string::size_type slashPos = 0;
            string name = entry->getName();
            if ((slashPos = name.rfind('/')) != string::npos) {
                seqWithoutCharge = name.substr(0, slashPos);
            } else {
                seqWithoutCharge = name;
            }
        }

        found = m_probTable->find(seqWithoutCharge);
    }

    /* Feature to extract entries with specified raw spectra. Commented out for efficiency.
    if (found == m_probTable->end()) {
      // also consider case when the prob table contains raw spectrum strings
      string rawSpectra = entry->getRawSpectra();
      string::size_type commaPos = 0;
      do {
        string spectrum = nextToken(rawSpectra, commaPos, commaPos, ",\t\r\n", ",");
      } while (!spectrum.empty() && (found = m_probTable->find(spectrum)) == m_probTable->end());
    }
    */

    // still not found
    if (found == m_probTable->end()) {
        return (false);
    }

    // found in prob table, now update prob if present
    if (updateProb && found->second >= 0.0) {
        stringstream probss;
        probss.precision(4);
        probss << fixed << found->second;
        entry->setOneComment("Prob", probss.str());
    }

    return (true);

}

// readProteinList - reads a protein list from a file. A protein list consists of lines in the format <protein> <whitespace> <cap>,
// where <cap> is the maximum number of peptide ions to be retained for that protein. If there are more than <cap> peptide ions
// associated with that protein, the entries will be added in a first-come-first-serve basis. More usefully, the library should
// first be sorted by the number of replicates (-cAN option), then the peptide ions with more replicates (more proteotypic) will be 
// retained ahead of those with fewer. Entries without any of its associated protein(s) in the protein list will be discarded.
void SpectraSTLibImporter::readProteinList() {

    ifstream fin;
    if (!myFileOpen(fin, m_params.useProteinList)) {
        g_log->error("CREATE", "Cannot read protein list. No filtering by protein will be performed.");
        return;
    }

    m_proteinList = new map<string, int>;

    string line("");
    while (nextLine(fin, line, "_EOF_", "")) {
        string::size_type pos = 0;
        string protein = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
        if (protein.empty() || protein[0] == '#') {
            continue;
        }
        string capStr = nextToken(line, pos, pos, " \t\r\n", " \t\r\n");
        int cap = 999999;
        if (!capStr.empty() && capStr[0] >= '0' && capStr[0] <= '9') {
            cap = atoi(capStr.c_str());
        }

        map<string, int>::iterator found = m_proteinList->find(protein);

        if (found != m_proteinList->end()) {
            if (found->second > cap) found->second = cap;
        } else {
            (*m_proteinList)[protein] = cap;
        }
    }

    // put together a string to describe filtering by prob table
    stringstream fss;
    fss << "FILTER for entries in protein list \"" << m_params.useProteinList << "\"";
    g_log->log("CREATE", "Apply " + fss.str());

    if (!m_filterStr.empty()) {
        m_filterStr += "; ";
    }
    m_filterStr += fss.str();


}

// isInProteinList - returns true if the entry is associated with any of the proteins in the protein list
bool SpectraSTLibImporter::isInProteinList(SpectraSTLibEntry *entry) {

    if (!m_proteinList) {
        return (true);
    }

    string proteins("");
    if (entry->getOneComment("Protein", proteins)) {

        string::size_type slashPos = proteins.find('/', 0);
        string protein("");

        if (slashPos == string::npos) {
            // single protein format
            map<string, int>::iterator found = m_proteinList->find(proteins);
            if (found != m_proteinList->end() && found->second > 0) {
                found->second--;
                return (true);
            } else {
                return (false);
            }
        } else {
            while (!((protein = nextToken(proteins, slashPos + 1, slashPos, "/\t\r\n")).empty())) {
                map<string, int>::iterator found = m_proteinList->find(protein);
                if (found != m_proteinList->end() && found->second > 0) {
                    found->second--;
                    return (true);
                }
            }
            return (false);
        }
    }
    return (true);
}

// printProteinList - outputs the protein list to cerr (debugging mostly)
void SpectraSTLibImporter::printProteinList() {

    if (!m_proteinList) return;

    for (map<string, int>::iterator i = m_proteinList->begin(); i != m_proteinList->end(); i++) {
        cerr << i->first << '\t' << i->second << endl;
    }
}

bool SpectraSTLibImporter::isAllDecoyProteins(SpectraSTLibEntry *entry) {

    if (m_params.removeDecoyProteins.empty()) {
        return (false);
    }

    vector<string> proteins;

    entry->getAllProteins(proteins);

    if (proteins.size() == 0) return (false);

    bool isAllDecoys = true;

    for (vector<string>::iterator i = proteins.begin(); i != proteins.end(); i++) {
        if (i->substr(0, m_params.removeDecoyProteins.length()) != m_params.removeDecoyProteins) {
            isAllDecoys = false;
            break;
        }
    }

    return (isAllDecoys);

}

// constructOutputFileName - if the output file name is not specified, this method can be used to piece together an
// output file name that is somewhat descriptive of what's happening
string SpectraSTLibImporter::constructOutputFileName() {

    stringstream ss;

    FileName fn0;
    parseFileName(m_impFileNames[0], fn0);
    ss << fn0.path << fn0.name;

    char oper = 'U';  // by default it's a UNION operation joining all imported files together

    if (m_impFileNames.size() == 1) {
        // nothing to do

    } else if (m_impFileNames.size() < 4) {
        // build a chain name if less than 4 files
        for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
            FileName fn;
            parseFileName(*i, fn);
            ss << '_' << oper << '_' << fn.name;
        }

    } else {
        // if more than 3 files, use the first name + operator + "plus"
        ss << '_' << oper << "_plus";
    }

    if (fn0.ext == ".splib" && m_impFileNames.size() == 1) {
        // single file. add a suffix so that it won't overwrite the file in a splib import
        ss << "_new";
    }

    ss << ".splib";

    return (ss.str());

}


// createSpectraSTLibImporter - this is the "factory method", a static instantiation method to return the proper
// importer object. If you subclass, REMEMBER to modify this method to call your constructor given the right extension!
SpectraSTLibImporter *SpectraSTLibImporter::createSpectraSTLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                                       SpectraSTCreateParams &params) {

    if (impFileNames.empty()) return (NULL);

    // type determined by first file
    FileName fn;
    parseFileName(impFileNames[0], fn);

    // check if all files have some extension
    for (vector<string>::iterator i = impFileNames.begin() + 1; i != impFileNames.end(); i++) {
        string ext;
        string::size_type dummy = getExtension(*i, ext);
        if (ext != fn.ext) {
            g_log->error("CREATE", "Libraries to be imported have mixed formats. Exiting.");
            return (NULL);
        }
    }

    if (fn.ext == ".msp") {
        return (new SpectraSTMspLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".hlf") {
        return (new SpectraSTXHunterLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".pepXML" || fn.ext == ".pep.xml" || fn.ext == ".xml") {
        return (new SpectraSTPepXMLLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".splib") {
        return (new SpectraSTSpLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".ms2") {
        return (new SpectraSTMs2LibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".tsv" || fn.ext == ".xls") {
        return (new SpectraSTTsvLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".mzXML" || fn.ext == ".mzML") {
        return (new SpectraSTMzXMLLibImporter(impFileNames, lib, params));
    } else if (fn.ext == ".sptxt") {
        // treat as .msp unless there is a corresponding .splib in the directory
        if (impFileNames.size() > 1) {
            g_log->error("CREATE", "Trying to import more than one .sptxt (as .msp) files. Exiting.");
            return (NULL);
        }
        string splibFile = fn.path + fn.name + ".splib";
        ifstream splibFileFin(splibFile.c_str());
        if (splibFileFin.good()) {
            // splib file already exists
            g_log->error("CREATE", "Trying to import a .sptxt (as .msp) but the corresponding .splib exists. Exiting.");
            splibFileFin.close();
            return (NULL);
        } else {
            // splib file not there, allow import and treat as .msp
            return (new SpectraSTMspLibImporter(impFileNames, lib, params));
        }

    } else {
        g_log->error("CREATE", "Libraries to be imported have illegal formats. Exiting.");
        return (NULL);
    }
    return (NULL);
}

Peptide *
SpectraSTLibImporter::createPeptide(string peptide, int charge, string modStr, string spectrum, string fileType) {

    Peptide *pep = new Peptide(peptide, charge, modStr);

    if (pep->hasUnknownMod) {
        g_log->error(fileType + " IMPORT",
                     "Peptide ID has unknown modification: \"" + peptide + "\". Skipped spectrum " + spectrum + " .");
        delete (pep);
        return (NULL);
    }

    if (pep->illegalPeptideStr) {
        g_log->error(fileType + " IMPORT",
                     "Illegal peptide ID string: \"" + peptide + "\". Skipped spectrum " + spectrum + " .");
        delete (pep);
        return (NULL);
    }

    unsigned int naa = pep->NAA();
    if (naa < m_params.minimumNumAAToInclude) {
        g_log->log(fileType + " IMPORT", "Peptide ID too short. Skipped spectrum " + spectrum + " .");
        delete (pep);
        return (NULL);
    }

    return (pep);
}

bool SpectraSTLibImporter::insertOneEntry(SpectraSTLibEntry *entry, string fileType) {

    SpectraSTPeakList *peakList = entry->getPeakList();

    if (!(m_params.keepRawIntensities)) {
        peakList->normalizeTo(10000, m_params.rawSpectraMaxDynamicRange);
    }

    if (m_params.centroidPeaks) peakList->centroid("TOF");

    if (peakList->getNumPeaks() < m_params.minimumNumPeaksToInclude) {
        g_log->log(fileType + " IMPORT", "Too few peaks in spectrum. Skipped entry \"" + entry->getName() + "\".");
        return (false);
    }

    if (!(m_params.setFragmentation.empty())) {
        // if user specifies something, always overwrite
        entry->setFragType(m_params.setFragmentation);
    } else {
        if (entry->getFragType().empty()) {
            string fragType("CID");
            entry->setFragType(fragType); // default to CID -- better than leaving it blank!
        }
    }

    if (passAllFilters(entry)) {

        if (!(m_params.skipRawAnnotation)) entry->annotatePeaks(true, false);

        m_lib->insertEntry(entry);
        return (true);

    } else {

        g_log->log(fileType + " IMPORT",
                   "Entry did not pass user-defined filter. Skipped entry \"" + entry->getName() + "\".");
        return (false);
    }

}

// setDeamidatedNXST - sets any N in NXST motifs to be Deamidated (for glycocaptured data)
void SpectraSTLibImporter::setDeamidatedNXST(Peptide *pep) {

    string::size_type Npos = 0;
    map<int, string>::iterator foundDeamidation;
    while (Npos < pep->stripped.length() - 2 && (Npos = pep->stripped.find('N', Npos)) != string::npos) {
        if (!pep->isModsSet) {
            break;
        }
        foundDeamidation = pep->mods.find((int) Npos);
        if (foundDeamidation == pep->mods.end()) {
            if (Npos + 2 < pep->stripped.length() &&
                (pep->stripped[Npos + 2] == 'S' || pep->stripped[Npos + 2] == 'T')) {
                // not already modified, and an S or T two sites later
                pep->mods[(int) Npos] = "Deamidated";
            }
        } else if (foundDeamidation->second == "Deamidated") {
            if (Npos + 2 < pep->stripped.length() &&
                (pep->stripped[Npos + 2] != 'S' && pep->stripped[Npos + 2] != 'T')) {
                // modified at the wrong site
                pep->mods.erase(foundDeamidation);
            }
        }
        Npos++;
    }

    if (pep->stripped[pep->stripped.length() - 2] == 'N') {
        foundDeamidation = pep->mods.find((int) (pep->stripped.length()) - 2);
        if (foundDeamidation == pep->mods.end()) {
            if (pep->nextAA == 'S' || pep->nextAA == 'T') {
                pep->mods[(int) pep->stripped.length() - 2] = "Deamidated";
            }
        } else if (foundDeamidation->second == "Deamidated") {
            if (pep->nextAA != 'S' && pep->nextAA != 'T') {
                pep->mods.erase(foundDeamidation);
            }
        }
    }

    // There's also the case where the N is the last AA on the C-terminus.
    // However, since we have no knowledge what AA is 2-residue down from this N, we assume that
    // this N is not part of a NXS/T motif.
    // In the grand scheme of things this probably shouldn't matter a whole lot, because this peptide
    // would be semitryptic and it's somewhat rare anyway.

    if (pep->stripped[pep->stripped.length() - 1] == 'N') {
        foundDeamidation = pep->mods.find((int) (pep->stripped.length()) - 1);
        if (foundDeamidation != pep->mods.end() && foundDeamidation->second == "Deamidated") {
            pep->mods.erase(foundDeamidation);
        }
    }

}

