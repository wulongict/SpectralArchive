#include "SpectraSTLibEntry.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

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

/* Class: SpectraSTLibEntry
 * 
 * The class that represents a library entry. 
 * 
 */


using namespace std;

extern bool g_verbose;
extern SpectraSTLog *g_log;

// constructor from arguments
SpectraSTLibEntry::SpectraSTLibEntry(Peptide *pep, string comments, string status, SpectraSTPeakList *peakList,
                                     string fragType) :
        m_name(""),
        m_pep(pep),
        m_charge(pep->charge),
        m_mw(pep->monoisotopicMH()),
        m_precursorMz(pep->monoisotopicMZ()),
        m_commentsStr(comments),
        m_comments(NULL),
        m_status(status),
        m_fullName(""),
        m_peakList(peakList),
        m_libId(0),
        m_libFileOffset(0),
        m_dotTimeProfile(NULL),
        m_ms1(NULL),
        m_fragType(fragType) {

    // the name is of the format AC[339]DEFGHIK/2
    m_name = pep->interactStyleWithCharge();

    // the full name is of the format K.AC[339]DEFGHIK.L/2 (ETD)
    m_fullName = pep->interactStyleFullWithCharge();

    if (!(m_fragType.empty())) {
        m_fullName += " (" + m_fragType + ")";
    }

    // append a comment with the average precursor m/z
    stringstream avemzss;
    avemzss << " AvePrecursorMz=";
    avemzss.precision(4);
    avemzss << fixed << showpoint << pep->averageMZ();
    m_commentsStr += avemzss.str();

    // if peakList is not instantiated, will do so here
    if (!m_peakList) {
        m_peakList = new SpectraSTPeakList(m_precursorMz, m_charge, 0, false, m_fragType);
    }

    // tell the peak list what peptide ion it is
    m_peakList->setPeptidePtr(pep);

}

// constructor from arguments
SpectraSTLibEntry::SpectraSTLibEntry(string name, double precursorMz, string comments, string status,
                                     SpectraSTPeakList *peakList, string fragType) :
        m_name(name),
        m_pep(NULL),
        m_charge(1),
        m_mw(precursorMz),
        m_precursorMz(precursorMz),
        m_commentsStr(comments),
        m_comments(NULL),
        m_status(status),
        m_fullName(name),
        m_peakList(peakList),
        m_libId(0),
        m_libFileOffset(0),
//  m_dotHistogram(50, 0.0),
        m_dotTimeProfile(NULL),
        m_ms1(NULL),
        m_fragType(fragType) {

    if (!(m_fragType.empty())) {
        m_fullName += " (" + m_fragType + ")";
    }

    string::size_type slashPos = m_name.rfind('/');
    if (slashPos != string::npos && slashPos < m_name.length() - 1) {
        m_charge = atoi(m_name.substr(slashPos + 1).c_str());
        m_mw = m_precursorMz * m_charge;
    }

    // if peakList is not instantiated, will do so here
    if (!m_peakList) {
        m_peakList = new SpectraSTPeakList(m_precursorMz, m_charge, 0, false, m_fragType);
    }

}


// constructor from library file
SpectraSTLibEntry::SpectraSTLibEntry(ifstream &libFin, bool binary, bool forSearch) :
        m_name(""),
        m_charge(0),
        m_mw(0.0),
        m_precursorMz(0.0),
        m_commentsStr(""),
        m_pep(NULL),
        m_comments(NULL),
        m_status(""),
        m_fullName(""),
        m_peakList(NULL),
        m_libId(0),
        m_libFileOffset(0),
        m_dotTimeProfile(NULL),
        m_ms1(NULL),
        m_fragType("") {

    // construct the object by reading from a library file
    if (binary) {
        readFromBinaryFile(libFin, forSearch);
    } else {
        readFromFile(libFin, forSearch);
    }
}

// readFromFile - reads from a text .splib file.
void SpectraSTLibEntry::readFromFile(ifstream &libFin, bool forSearch) {

    m_libFileOffset = libFin.tellg();

    string mods("");

    string line;
    while (nextLine(libFin, line, "Num")) { // read until the NumPeaks line
        string::size_type separatorPos = 0;
        string head = nextToken(line, 0, separatorPos, ":");
        string tail = nextToken(line, separatorPos + 1, separatorPos, "\r\n", " \t");

        if (head == "Name") {
            m_name = tail;
        } else if (head == "LibID") {
            m_libId = atoi(tail.c_str());
        } else if (head == "MW") {
            m_mw = atof(tail.c_str());
        } else if (head == "PrecursorMZ") {
            m_precursorMz = atof(tail.c_str());
        } else if (head == "Charge") { // deprecated - newer library doesn't have this - charge is specified in the name
            m_charge = atoi(tail.c_str());
        } else if (head ==
                   "Mods") { // deprecated - newer library doesn't have this - mod is specified in the peptide name with []
            mods = tail;
        } else if (head == "FullName") {
            m_fullName = tail;
        } else if (head == "Status") {
            m_status = tail;
        } else if (head == "Comment") {
            m_commentsStr = tail;
        }
    }

    if (line == "_EOF_") {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }

    if (m_name.empty()) {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }

    if (m_name[0] == '_') {
        // not a peptide
        m_charge = 1;
        m_pep = NULL;

        string fullNameMinusFrag(m_fullName);
        string::size_type spacePos = m_fullName.find(" (");
        if (spacePos != string::npos) {

            if (spacePos + 2 < m_fullName.length()) {
                string::size_type closeParenPos = m_fullName.find(')', spacePos + 2);
                if (closeParenPos != string::npos) {
                    m_fragType = m_fullName.substr(spacePos + 2, closeParenPos - spacePos - 2);
                }
            }

            fullNameMinusFrag = m_fullName.substr(0, spacePos);
        }

        string::size_type slashPos = m_name.rfind('/');
        if (slashPos != string::npos && slashPos < m_name.length() - 1) {
            m_charge = atoi(m_name.substr(slashPos + 1).c_str());
        }

    } else {

        if (m_fullName.empty()) {
            // if fullName is not provided (older version of .splib), find it in the Comments
            string mspFullName("");
            if (getOneComment("Fullname", mspFullName)) {
                m_pep = new Peptide(m_name, m_charge, mods);
                m_pep->prevAA = mspFullName[0];
                m_pep->nextAA = mspFullName[mspFullName.rfind('/') - 1];
            } else {
                m_pep = new Peptide(m_name, m_charge, mods);
            }
        } else {
            string fullNameMinusFrag(m_fullName);
            string::size_type spacePos = m_fullName.find(" (");
            if (spacePos != string::npos) {

                if (spacePos + 2 < m_fullName.length()) {
                    string::size_type closeParenPos = m_fullName.find(')', spacePos + 2);
                    if (closeParenPos != string::npos) {
                        m_fragType = m_fullName.substr(spacePos + 2, closeParenPos - spacePos - 2);
                    }
                }

                fullNameMinusFrag = m_fullName.substr(0, spacePos);
            }

            m_pep = new Peptide(fullNameMinusFrag, m_charge, mods);

            // if (!(m_pep->isGood())) {
            //   g_log->error("GENERAL", "Problematic peptide ID (Unrecognized modification?): " + m_fullName + " .");
            //   g_log->crash();
            // }

        }

        m_charge = m_pep->charge;
        m_name = m_pep->interactStyleWithCharge();

    }

    if (m_pep) {
        // if mw and precursor m/z not specified, calculate from peptide
        if (m_mw < 0.00001) {
            m_mw = m_pep->monoisotopicMH();
        }
        if (m_precursorMz < 0.00001) {
            m_precursorMz = m_pep->monoisotopicMZ();
        }
    } else {
        if (m_mw < 0.00001) {
            m_mw = m_precursorMz * m_charge;
        }
    }

    // line should start with "NumPeaks:" now

    string::size_type separatorPos = 0;
    string dummy = nextToken(line, 0, separatorPos, ":");
    int numPeaks = atoi(nextToken(line, separatorPos + 1, separatorPos, "\r\n", " \t").c_str());

    m_peakList = new SpectraSTPeakList(m_precursorMz, m_charge, numPeaks, true, m_fragType);

    if (m_pep) {
        m_peakList->setPeptidePtr(m_pep);
    }
    // read peak list
    string peakLine;
    while (nextLine(libFin, peakLine, "")) {

        string::size_type pos = 0;
        double mz = atof((nextToken(peakLine, pos, pos)).c_str());
        float intensity = atof((nextToken(peakLine, pos, pos)).c_str());
        string annotation = nextToken(peakLine, pos, pos);
        string info = nextToken(peakLine, pos, pos, "\r\n"); // read till end of line

        if (!forSearch) {
            m_peakList->insert(mz, intensity, annotation, info);
        } else {
            m_peakList->insertForSearch(mz, intensity, (annotation.empty() ? "" : annotation.substr(0, 1)));
        }
    }
    if (numPeaks != m_peakList->getNumPeaks()) {
        cerr << "SpectraSTLibEntry: Number of peaks listed does not match specified number. Ignore specified number."
             << endl;
    }

}

// readFromBinaryFile - reads an entry from a binary library file
// the binary format is:
// <libId (int)> <fullName (\n-terminated string)> <precursorM/Z (double)>
// <status (\n-terminated string)> <numPeaks (int)>
// numPeaks times: <m/z (double)> <intensity (double)> <annotation (\n-terminated string)> <info (\n-terminated string)>
// <comment (\n-terminated string)>

void SpectraSTLibEntry::readFromBinaryFile(ifstream &libFin, bool forSearch) {

    m_libFileOffset = libFin.tellg();

    string line("");

    // lib ID
    libFin.read((char *) (&m_libId), sizeof(int));

    // full name
    if (!nextLine(libFin, line)) {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }
    m_fullName = line;

    if (m_fullName.size() < 2) {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }

    if (m_fullName[0] == '_') {
        // not a peptide
        m_charge = 1;
        m_pep = NULL;
        m_name = m_fullName;

        string fullNameMinusFrag(m_fullName);
        string::size_type spacePos = m_fullName.find(" (");
        if (spacePos != string::npos) {

            if (spacePos + 2 < m_fullName.length()) {
                string::size_type closeParenPos = m_fullName.find(')', spacePos + 2);
                if (closeParenPos != string::npos) {
                    m_fragType = m_fullName.substr(spacePos + 2, closeParenPos - spacePos - 2);
                }
            }

            fullNameMinusFrag = m_fullName.substr(0, spacePos);
        }

        string::size_type slashPos = m_name.rfind('/');
        if (slashPos != string::npos && slashPos < m_name.length() - 1) {
            m_charge = atoi(m_name.substr(slashPos + 1).c_str());
        }

    } else {

        if (m_fullName[1] != '.') {
            g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
            g_log->crash();
        }

        string fullNameMinusFrag(m_fullName);
        string::size_type spacePos = m_fullName.find(" (");
        if (spacePos != string::npos) {

            if (spacePos + 2 < m_fullName.length()) {
                string::size_type closeParenPos = m_fullName.find(')', spacePos + 2);
                if (closeParenPos != string::npos) {
                    m_fragType = m_fullName.substr(spacePos + 2, closeParenPos - spacePos - 2);
                }
            }

            fullNameMinusFrag = m_fullName.substr(0, spacePos);
        }

        m_pep = new Peptide(fullNameMinusFrag, 0, "");

        //   if (!(m_pep->isGood())) {
        //     g_log->error("GENERAL", "Problematic peptide ID (Unrecognized modification?): " + m_fullName + " . Please check that your usermods file is in the current directory.");
        //     g_log->crash();
        //   }

        m_name = m_pep->interactStyleWithCharge();

        m_charge = m_pep->charge;

    }
    // precursor m/z
    libFin.read((char *) (&m_precursorMz), sizeof(double));

    m_mw = m_precursorMz * m_charge;

    // status
    if (!nextLine(libFin, line)) {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }
    m_status = line;

    // num peaks
    unsigned int numPeaks = 0;
    libFin.read((char *) (&numPeaks), sizeof(unsigned int));

    m_peakList = new SpectraSTPeakList(m_precursorMz, m_charge, numPeaks, true, m_fragType);

    if (m_pep) {
        m_peakList->setPeptidePtr(m_pep);
    }

    // peaks
    for (unsigned int i = 0; i < numPeaks; i++) {
        double mz = 0.0;
        double intensity = 0.0;
        libFin.read((char *) (&mz), sizeof(double));
        libFin.read((char *) (&intensity), sizeof(double));

        if (!nextLine(libFin, line)) {
            g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
            g_log->crash();
        }
        string annotation = line;

        if (!nextLine(libFin, line)) {
            g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
            g_log->crash();
        }
        string info = line;

        float floatIntensity = (float) intensity;

        if (!forSearch) {
            m_peakList->insert(mz, floatIntensity, annotation, info);
        } else {
            m_peakList->insertForSearch(mz, floatIntensity, (annotation.empty() ? "" : annotation.substr(0, 1)));
        }
    }

    // comments
    if (!nextLine(libFin, line)) {
        g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
        g_log->crash();
    }
    m_commentsStr = line;

}

// destructor
SpectraSTLibEntry::~SpectraSTLibEntry() {

    if (m_peakList) {
        delete m_peakList;
        m_peakList = NULL;
    }

    if (m_pep) {
        delete m_pep;
        m_pep = NULL;
    }

    if (m_comments) {
        delete m_comments;
        m_comments = NULL;
    }

    if (m_dotTimeProfile) {
        delete (m_dotTimeProfile);
        m_dotTimeProfile = NULL;
    }

    if (m_ms1) {
        delete (m_ms1);
        m_ms1 = NULL;
    }

}

// copy constructor
SpectraSTLibEntry::SpectraSTLibEntry(SpectraSTLibEntry &other) :
        m_pep(NULL),
        m_comments(NULL),
        m_peakList(NULL),
        m_ms1(NULL),
        m_dotTimeProfile(NULL) {

    (*this) = other;
}

// assignment operator
SpectraSTLibEntry &SpectraSTLibEntry::operator=(SpectraSTLibEntry &other) {

    this->m_libId = other.m_libId;
    this->m_libFileOffset = other.m_libFileOffset;
    this->m_name = other.m_name;
    this->m_mw = other.m_mw;
    this->m_commentsStr = other.m_commentsStr;
    this->m_precursorMz = other.m_precursorMz;
    this->m_fullName = other.m_fullName;
    this->m_charge = other.m_charge;
    this->m_status = other.m_status;
    this->m_fragType = other.m_fragType;

    // deep copy m_pep
    if (this->m_pep) {
        delete (this->m_pep);
        this->m_pep = NULL;
    }
    if (other.m_pep) {
        this->m_pep = new Peptide(*(other.m_pep));
    } else {
        this->m_pep = NULL;
    }

    // deep copy m_comments
    if (this->m_comments) {
        delete (this->m_comments);
        this->m_comments = NULL;
    }
    if (other.m_comments) {
        this->m_comments = new map<string, string>;
        for (map<string, string>::iterator c = other.m_comments->begin(); c != other.m_comments->end(); c++) {
            this->m_comments->insert(*c);
        }
    } else {
        this->m_comments = NULL;
    }

    // deep copy peakList
    if (this->m_peakList) {
        delete (this->m_peakList);
        this->m_peakList = NULL;
    }
    if (other.m_peakList) {
        this->m_peakList = new SpectraSTPeakList(*(other.m_peakList));

        // set peakList's m_pep
        this->m_peakList->setPeptidePtr(this->m_pep);

    } else {
        this->m_peakList = NULL;
    }


    // deep copy m_ms1
    if (this->m_ms1) {
        delete (this->m_ms1);
        this->m_ms1 = NULL;
    }
    if (other.m_ms1) {
        this->m_ms1 = new SpectraSTPeakList(*(other.m_ms1));
    } else {
        this->m_ms1 = NULL;
    }

    // deep copy m_dotTimeProfile
    if (this->m_dotTimeProfile) {
        delete (this->m_dotTimeProfile);
        this->m_dotTimeProfile = NULL;
    }
    if (other.m_dotTimeProfile) {
        this->m_dotTimeProfile = new vector<pair<int, float> >(*(other.m_dotTimeProfile));
    } else {
        this->m_dotTimeProfile = NULL;
    }

    return (*this);
}

// getSafeName - a name that can be used as part of a file name
string SpectraSTLibEntry::getSafeName() {

    if (m_pep) {

        stringstream ss;
        ss << m_pep->getSafeName() << '_' << m_libId;
        return (ss.str());
    } else {
        stringstream safeName;
        string forbidden("~`!@#$%^&*()+\\|:<>,.?/{}'\"");
        for (string::size_type pos = 0; pos < m_name.length(); pos++) {
            if (forbidden.find(m_name[pos]) != string::npos) {
                safeName << '_';
            } else {
                safeName << m_name[pos];
            }
        }

        return (safeName.str());
    }

}

void SpectraSTLibEntry::setPrecursor(double precursorMz, int charge) {

    m_precursorMz = precursorMz;
    if (charge != 0) m_charge = charge;

    m_mw = m_precursorMz * (m_charge == 0 ? 1 : m_charge);

    m_peakList->setParentMz(m_precursorMz);
}

void SpectraSTLibEntry::setFragType(string &fragType) {

    m_fragType = fragType;
    m_peakList->setFragType(fragType);

    string::size_type spacePos = m_fullName.find(" (", 0);

    if (spacePos != string::npos) {
        string pre = m_fullName.substr(0, spacePos);
        if (fragType.empty()) {
            m_fullName = pre;
        } else {
            m_fullName = pre + " (" + m_fragType + ")";
        }
    } else {
        if (!(fragType.empty())) {
            m_fullName += " (" + m_fragType + ")";
        }
    }
}

// writeToFile - writes the entry to a file in text format
void SpectraSTLibEntry::writeToFile(ofstream &libFout) {
    libFout << "Name: " << m_name << endl;
    libFout << "LibID: " << m_libId << endl;
    libFout.precision(4);
    libFout << "MW: " << fixed << m_mw << endl;
    libFout.precision(4);
    libFout << "PrecursorMZ: " << fixed << m_precursorMz << endl;
    libFout << "Status: " << m_status << endl;
    libFout << "FullName: " << m_fullName << endl;
    libFout << "Comment: " << getCommentsStr() << endl;

    if (m_peakList) {
        m_peakList->writeToFile(libFout);
    } else {
        libFout << "NumPeaks: 0" << endl;
    }
    libFout << endl;

}

// writeToBinaryFile - write the entry to file in binary format
void SpectraSTLibEntry::writeToBinaryFile(ofstream &libFout) {

    libFout.write((char *) (&m_libId), sizeof(int));
    libFout << m_fullName << endl;
    libFout.write((char *) (&m_precursorMz), sizeof(double));
    libFout << m_status << endl;
    m_peakList->writeToBinaryFile(libFout);
    libFout << getCommentsStr() << endl;

}

// writeDtaFile - write the entry as a .dta (i.e. no header) 
void SpectraSTLibEntry::writeDtaFile(string dtaFileName) {

    ofstream dtaFout;

    if (!myFileOpen(dtaFout, dtaFileName)) {
        g_log->error("CREATE", "Cannot open DTA file \"" + dtaFileName +
                               "\" for writing the library spectrum to dta format. Skipped.");
        return;
    }

    if (m_pep) {
        dtaFout.precision(4);
        dtaFout << fixed << m_pep->monoisotopicNeutralM() + 1.00783 << '\t';  // SEQUEST needs M+1H
        dtaFout.precision(4);
        dtaFout << fixed << m_pep->charge << endl;
    } else {
        dtaFout.precision(4);
        dtaFout << fixed << m_precursorMz * m_charge - m_charge + 1.00783 << '\t';  // SEQUEST needs M+1H
        dtaFout.precision(4);
        dtaFout << fixed << m_charge << endl;
    }

    if (m_peakList) {
        m_peakList->writeToDtaFile(dtaFout);
    }

    dtaFout.flush();
    dtaFout.close();

}

void SpectraSTLibEntry::writeMgfFile(ofstream &mgfFout) {

    mgfFout << "BEGIN IONS" << endl;
    mgfFout.precision(6);
    mgfFout << "PEPMASS=" << fixed << getPrecursorMz() << endl;
    mgfFout << "CHARGE=" << getCharge() << endl;
    mgfFout << "TITLE=" << getName() << endl;

    if (m_peakList) {
        m_peakList->writeToDtaFile(mgfFout);
    }

    mgfFout << "END IONS" << endl;

}

// writeMRM - write the entry as MRM transition table
void SpectraSTLibEntry::writeMRM(ofstream &mrmFout, string format) {

    map<string, pair<unsigned int, unsigned int> > samples;
    getSampleInfo(samples, "USED_ONLY");

    int totNumUsed = 0;
    int maxNumUsed = 0;
    string bestSamp("");
    for (map<string, pair<unsigned int, unsigned int> >::iterator si = samples.begin(); si != samples.end(); si++) {
        int curNumUsed = si->second.first;
        totNumUsed += curNumUsed;
        if (curNumUsed > maxNumUsed) {
            maxNumUsed = curNumUsed;
            bestSamp = si->first;
        }
    }

    double pI = 0.0;

    if (m_pep) {
        pI = m_pep->computePI();
    }

    double rt = 0.0;
    string rtstr("");
    string::size_type pos = 0;
    if (getOneComment("RetentionTime", rtstr)) {
        rt = atof(nextToken(rtstr, 0, pos, ",/ \t\r\n").c_str());
    }
    double irt = 0.0;
    string irtstr("");
    pos = 0;
    if (getOneComment("iRT", irtstr)) {
        irt = atof(nextToken(irtstr, 0, pos, ",/ \t\r\n").c_str());
    }
    string protein("");
    stringstream proteinss;
    pos = 0;
    if (getOneComment("Protein", protein)) {
        if (protein.find('/') != string::npos) {
            int proteinCount = atoi(nextToken(protein, 0, pos, "/\t\r\n").c_str());
            proteinss << proteinCount << '\t';
            bool isFirstProtein = true;
            while (pos < protein.length()) {
                if (!isFirstProtein) proteinss << ';';
                proteinss << nextToken(protein, pos, pos, "/\t\r\n", "/");
                isFirstProtein = false;
            }
        } else {
            proteinss << '1' << '\t' << protein;
        }
    } else {
        proteinss << '0';
    }

    // make sure we provide the monoisotopic mass for the precursor
    double monoPrecursorMz = m_precursorMz;
    string dummy;
    if (!getOneComment("AvePrecursorMz", dummy)) {
        monoPrecursorMz = m_pep->monoisotopicMZ();
    }

    double collisionEnergy = 0;
    if (m_charge <= 2) {
        collisionEnergy = 0.044 * monoPrecursorMz + 5.5;
    } else {
        collisionEnergy = 0.051 * monoPrecursorMz + 0.55;
    }
    stringstream cess;
    cess.precision(2);
    cess << fixed << showpoint << collisionEnergy;

    if (format == "DEFAULT" || format == "SHOWINFO") {
        stringstream press;
        press << bestSamp << '\t' << maxNumUsed << '/' << totNumUsed << '\t';
        press.precision(2);
        press << fixed << pI << '\t';
        press.precision(4);
        press << fixed << monoPrecursorMz << '\t';
        press.precision(2);
        press << fixed << rt;
        if (!(irtstr.empty())) {
            press.precision(2);
            press << fixed << '(' << irt << ')';
        }
        stringstream postss;
        postss << m_charge << '\t' << m_name << '\t' << cess.str() << '\t' << proteinss.str();
        m_peakList->writeMRM(mrmFout, press.str(), postss.str(), format);
    }

}

void SpectraSTLibEntry::writeInfo(ofstream &mrmFout) {

    double rt = 0.0;
    string rtstr("");
    string::size_type pos = 0;
    if (getOneComment("RetentionTime", rtstr)) {
        rt = atof(nextToken(rtstr, 0, pos, "/ \t\r\n").c_str());
        rt /= 60.0;
    }

    string protein("");
    stringstream proteinss;
    pos = 0;
    bool isDecoy = true;
    if (getOneComment("Protein", protein)) {
        if (protein.find('/') != string::npos) {
            int proteinCount = atoi(nextToken(protein, 0, pos, "/\t\r\n").c_str());
            proteinss << proteinCount << '\t';
            bool isFirstProtein = true;
            while (pos < protein.length()) {
                if (!isFirstProtein) proteinss << ';';
                string curProtein = nextToken(protein, pos, pos, "/\t\r\n", "/");
                if (curProtein.substr(0, 5) != "DECOY") isDecoy = false;
                proteinss << curProtein;
                isFirstProtein = false;
            }
        }
    }

    // make sure we provide the monoisotopic mass for the precursor
    double monoPrecursorMz = m_precursorMz;
    string dummy;
    if (!getOneComment("AvePrecursorMz", dummy)) {
        monoPrecursorMz = m_pep->monoisotopicMZ();
    }

    double massDiff = 0.0;
    if (getOneComment("MassDiff", dummy)) {
        massDiff = atof(dummy.c_str());
        monoPrecursorMz += (massDiff / (double) m_charge);
    }

    string rawSpectrum = getBestRawSpectrum();
    string::size_type lastDotPos = rawSpectrum.rfind('.');
    unsigned int scanNum = 9999999;
    if (lastDotPos != string::npos) scanNum = atoi(rawSpectrum.substr(lastDotPos + 1).c_str());
    string::size_type firstDotPos = rawSpectrum.find('.');
    string fileName("What");
    if (firstDotPos != string::npos) fileName = rawSpectrum.substr(0, firstDotPos);

    stringstream press;
    press << fileName << '\t' << scanNum - 1 << '\t';
    press.precision(4);
    press << fixed << monoPrecursorMz << '\t';
    press.precision(4);
    press << fixed << massDiff << '\t';
//  press << m_fragType << '\t';
    press << m_name << '\t';
    press << m_pep->mspMods() << '\t';
    press << proteinss.str() << '\t';
    press << (isDecoy ? "D" : "T") << '\t';
    press.precision(4);
    press << fixed << getProb() << '\t';

    mrmFout << press.str() << endl;

    stringstream postss;
    postss << m_charge << '\t' << m_name << '\t' << '\t' << proteinss.str();

//  mrmFout << press.str() << postss.str() << endl;

}


void SpectraSTLibEntry::writePAIdent(ofstream &fout, string baseName) {

    fout << "00000" << '\t'; // search ID

    fout << baseName << '.';
    fout.width(5);
    fout.fill('0');
    fout << right << m_libId + 1;
    fout << '.';
    fout.width(5);
    fout.fill('0');
    fout << right << m_libId + 1;
    fout << '.' << m_charge << '\t'; // spectrum reference

    fout << "PAp000000" << '\t'; // peptide accession

    fout << m_pep->stripped << '\t'; // peptide sequence

    fout << m_pep->prevAA << '\t'; // peptide prev AA

    fout << m_pep->interactStyle() << '\t'; // modified peptide

    fout << m_pep->nextAA << '\t'; // peptide next AA

    fout << m_charge << '\t'; // charge

    double prob = getProb();
    fout << prob << '\t'; // probability

    fout << "0.0000" << '\t'; // massdiff

    fout << getFirstProtein() << '\t'; // protein name

    fout << prob << '\t'; // NSP-adjusted probability

    fout << getNrepsUsed() << '\t'; // number of observations

    fout << "0.0000" << endl; // number of sibling peptides

}

string SpectraSTLibEntry::getFirstProtein(int &proteinCount) {

    string protein("");
    string::size_type pos = 0;
    if (getOneComment("Protein", protein)) {
        if (protein.find('/') != string::npos) {
            proteinCount = atoi(nextToken(protein, 0, pos, "/\t\r\n").c_str());
            return (nextToken(protein, pos, pos, "/\t\r\n", "/"));
        } else {
            // old format, no protein count
            proteinCount = 1;
            return (protein);
        }
    }
    proteinCount = 0;
    return (protein);
}

string SpectraSTLibEntry::getFirstProtein() {
    int dummy = 0;
    return (getFirstProtein(dummy));
}

void SpectraSTLibEntry::getAllProteins(vector<string> &proteins) {

    proteins.clear();

    string proteinStr("");
    string::size_type pos = 0;
    if (getOneComment("Protein", proteinStr)) {
        if (proteinStr.find('/') != string::npos) {
            unsigned int proteinCount = atoi(nextToken(proteinStr, 0, pos, "/\t\r\n").c_str());
            string protein("");
            while (!((protein = nextToken(proteinStr, pos, pos, "/\t\r\n", "/")).empty())) {
                proteins.push_back(protein);
            }
        } else {
            // old format, no protein count, single protein
            proteins.push_back(proteinStr);
        }
    }
}


string SpectraSTLibEntry::getRawSpectra() {

    string rawSpectra("");
    if (getOneComment("RawSpectrum", rawSpectra)) {
        return (rawSpectra);
    } else if (getOneComment("RawSpectra", rawSpectra)) {
        return (rawSpectra);
    }
    return ("");
}

string SpectraSTLibEntry::getBestRawSpectrum() {

    string rawSpectrum("");
    if (getOneComment("BestRawSpectrum", rawSpectrum)) {
        return (rawSpectrum);
    } else if (getOneComment("RawSpectrum", rawSpectrum)) {
        return (rawSpectrum);
    }
    return ("");
}


// parseCommentsStr - parses the comments string to create the m_comments object
// this is quite costly, so only do so when it needs to access and modify many comment fields,
// e.g. during consensus creation
void SpectraSTLibEntry::parseCommentsStr() {

    m_comments = new map<string, string>;

    string::size_type pos = 0;
    string field;

    while (!(field = nextToken(m_commentsStr, pos, pos, " \"\t\r\n")).empty()) {

        if (m_commentsStr[pos] == '\"') {
            // token ends with quote. this token is not complete,
            // should go all the way to the next quote
            // ignoring whitespace in between

            pos++; // into the quoted string
            string inQuote = nextToken(m_commentsStr, pos, pos, "\"");
            field += inQuote;
            pos++; // into the space following the end quote
        }


        string::size_type eqPos = 0;
        if ((eqPos = field.find('=', 0)) != string::npos) {
            // first equal sign separates the attribute and the value
            string attr = field.substr(0, eqPos);
            string value = field.substr(eqPos + 1);
            (*m_comments)[attr] = value;
        } else {
            // no equal sign, problem with token. just ignore
        }


        pos++;
    }

}

// getCommentsStr - returns the comments as a string; if m_comments is instantiated, it will
// prints that object as a string; if not, it will just return m_commentsStr.
string SpectraSTLibEntry::getCommentsStr() {

    if (!m_comments) {
        return (m_commentsStr);
    }

    stringstream ss;
    for (map<string, string>::iterator c = m_comments->begin(); c != m_comments->end(); c++) {

        if (c != m_comments->begin()) {
            ss << ' ';
        }

        ss << c->first << '=';

        if (c->second.find_first_of(" \t\r\n", 0) != string::npos) {
            // need quotes around it
            ss << '\"' << c->second << '\"';
        } else {
            ss << c->second;
        }
    }

    m_commentsStr = ss.str();
    return (ss.str());
}

// getOneComment - get one comment field. If m_comments is not instantiated, just
// do simple string parsing to find the field.
bool SpectraSTLibEntry::getOneComment(string attr, string &value) {

    if (!m_comments) {
        string::size_type found = m_commentsStr.find(attr + '=', 0);

        if (found == string::npos) {
            return (false);
        } else {

            string::size_type vs;
            if (found == 0) {
                vs = found + attr.length() + 1;
            } else {
                // make sure the attribute we're looking for isn't prefixed with something else
                found = m_commentsStr.find(' ' + attr + '=', 0);
                if (found == string::npos) {
                    return (false);
                }
                vs = found + attr.length() + 2;
            }

            if (vs >= m_commentsStr.length()) {
                return (false);
            }
            if (m_commentsStr[vs] == '\"') {
                // in quote, get rid of quotes
                value = nextToken(m_commentsStr, vs, vs, "\"\t\r\n", "\"");
            } else {
                value = nextToken(m_commentsStr, vs, vs, " \t\r\n");
            }
            return (true);
        }

    } else {
        // already parsed into m_comments. simply find it from the map object
        map<string, string>::iterator found = m_comments->find(attr);
        if (found == m_comments->end()) {
            return (false);
        } else {
            value = found->second;
            return (true);
        }
    }
}

// setOneComment - set one comment field. This will automatically parse the string into an object
// for easy insertion or modification. returns true if that comment is found.
bool SpectraSTLibEntry::setOneComment(string attr, string value) {

    if (!m_comments) {
        parseCommentsStr();
    }

    map<string, string>::iterator found = m_comments->find(attr);
    if (found == m_comments->end()) {
        (*m_comments)[attr] = value;
        return (false);
    } else {
        (*m_comments)[attr] = value;
        return (true);
    }

}

// deleteOneComment - delete one comment field. This will automatically parse the string into an object
// for easy deletion. returns true if that comment is found
bool SpectraSTLibEntry::deleteOneComment(string attr) {

    if (!m_comments) {
        parseCommentsStr();
    }

    map<string, string>::iterator found = m_comments->find(attr);
    if (found == m_comments->end()) {
        return (false);
    } else {
        m_comments->erase(found);
        return (true);
    }

}

// setPeakList - setting the peak list.
void SpectraSTLibEntry::setPeakList(SpectraSTPeakList *peakList) {

    if (m_peakList) {
        delete (m_peakList);
    }
    m_peakList = peakList;

}

void SpectraSTLibEntry::synchWithPep() {

    if (!m_pep) return;

    m_name = m_pep->interactStyleWithCharge();
    m_fullName = m_pep->interactStyleFullWithCharge();
    if (!(m_fragType.empty())) {
        m_fullName += " (" + m_fragType + ")";
    }
    m_mw = m_pep->monoisotopicMH();
    m_precursorMz = m_pep->monoisotopicMZ();
    m_charge = m_pep->charge;

    setOneComment("Mods", m_pep->mspMods());

    if (m_pep->NTT() == 2) {
        setOneComment("NTT", "2");
        setOneComment("Pep", "Tryptic");
    } else if (m_pep->NTT() == 1) {
        setOneComment("NTT", "1");
        setOneComment("Pep", "Semi-tryp");
    } else {
        setOneComment("NTT", "0");
        setOneComment("Pep", "Non-tryp");
    }

    stringstream nmcss;
    nmcss << m_pep->NMC();
    setOneComment("NMC", nmcss.str());

    stringstream naass;
    naass << m_pep->NAA();
    setOneComment("NAA", naass.str());


}

// getAveragePrecursorMz - get the isotopically averaged precursor m/z. If AvePrecursorMz is NOT found in the Comment, it
// means that the library is in the old format, that is, average mass is used for m_precursorMz.
double SpectraSTLibEntry::getAveragePrecursorMz() {

    double mz = m_precursorMz;
    string mzStr("");
    if (getOneComment("AvePrecursorMz", mzStr)) {
        mz = atof(mzStr.c_str());
    }
    return (mz);
}

// getProb - convenient getter for the probability in the comment
double SpectraSTLibEntry::getProb(double valueIfNotFound) {

    double prob = valueIfNotFound;

    string probStr("");
    if (getOneComment("Prob", probStr)) {
        prob = atof(probStr.c_str());
    } else {
        string tfRatioStr("");
        string::size_type tfRatioPos = 0;
        if (getOneComment("Tfratio", tfRatioStr)) {
            double tfRatio = atof(tfRatioStr.c_str());
            prob = tfRatio / (1.0 + tfRatio);
        }
    }

    return (prob);
}

// getNrepsUsed - convenient getter for the number of replicates used, i.e. the 'X' in 'Nrep=X/Y' in the comment.
unsigned int SpectraSTLibEntry::getNrepsUsed(unsigned int valueIfNotFound) {

    unsigned int numRepsUsed = valueIfNotFound;
    string nrepsStr("");
    string::size_type slashPos = 0;
    if (getOneComment("Nreps", nrepsStr)) {
        numRepsUsed = atoi(nextToken(nrepsStr, 0, slashPos, "/ \t\r\n").c_str());
    }

    return (numRepsUsed);
}

// getSeqInfo - parses out the Se= field in the comments, 
// and incorporates the result in the 'seqs' map. Note that if 'seqs' already contains
// some sequence search information, getSeqInfo will not delete it; instead, it merges
// that existing information with what it parses out from this comment.
// (This is useful if you want to combine multiple entries -- see SpectraSTSpLibImporter.cpp)
// the structure of the seqs map is: 
// char_denoting_search_engine => ptr of ( score_name => (A, B) ) where A is the sum of scores weighted
// by N (the number of replicates), B is the sum of ( <stdev>^2 * <mean>^2 ) weighted by N.
// To calculate the mean and stdev: mean = A / (sum of N), stdev = sqrt(B / (sum of N) - A^2). 
// This is done
// in SpectraSTLibEntry::setSeqInfo.

void SpectraSTLibEntry::getSeqInfo(map<char, map<string, pair<double, double> > *> &seqs) {

    string seqInfo("");
    if (!getOneComment("Se", seqInfo)) {
        // can't find Se field. nothing can be done.
        return;
    }

    // The structure of the Se= field is:
    // <num_engines>/^<char_denoting_search_engine1><num_reps_id_by_search_engine1>:<score1_name>=<score1_mean>/<score1_stdev>,
    // <score2_name>=<score2_mean>/<score2_stdev>,...^<char_denoting_search_engine2>...
    // If confused, see any .splib file created

    string::size_type hatPos = seqInfo.find('^', 0) + 1;
    string seq("");

    while (!((seq = nextToken(seqInfo, hatPos, hatPos, "^\r\n")).empty())) {

        // seq should contain the piece between two ^'s. this should correspond to
        // the info for one search engine

        char se = seq[0];
        string::size_type colonPos = 0;

        // N is the number of replicates identified by that engine
        double N = atof(nextToken(seq, 1, colonPos, ":\r\n").c_str());

        // try to find that search engine in the existing seqs
        map<char, map<string, pair<double, double> > *>::iterator foundse = seqs.find(se);
        map<string, pair<double, double> > *scs;

        if (foundse == seqs.end()) {
            // not found. create a new search engine (with no score for now)
            scs = new map<string, pair<double, double> >;
            seqs[se] = scs;
        } else {
            // found
            scs = (*foundse).second;
        }

        string::size_type commaPos = colonPos + 1;
        string sc("");

        // now parses the scores
        while (!((sc = nextToken(seq, commaPos, commaPos, ",\r\n")).empty())) {
            string::size_type pos = 0;
            string scName = nextToken(sc, 0, pos, "=\r\n");
            double value = atof(nextToken(sc, pos + 1, pos, "/\r\n").c_str());
            double dev = atof(nextToken(sc, pos + 1, pos, ",\r\n").c_str());

            // try to find that score in the existing scs map
            map<string, pair<double, double> >::iterator foundsc = (*scs).find(scName);

            if (foundsc == (*scs).end()) {

                // not found. create a new score
                pair<double, double> p;
                p.first = value * (double) N;
                p.second = (dev * dev + value * value) * (double) N;
                (*scs)[scName] = p;
            } else {
                // found. update that score
                (*foundsc).second.first += (value * (double) N);
                (*foundsc).second.second += ((dev * dev + value * value) * (double) N);
            }


            commaPos++;
        }

        // each sequence search engine will also carry a "nr" (num of replicates, N) score
        // to keep track of how many replicates are identified by this engine.
        map<string, pair<double, double> >::iterator foundnr = (*scs).find("nr");
        if (foundnr == (*scs).end()) {
            // not found. create

            pair<double, double> nr;
            nr.first = N;    // note that to make use of the structure, N is stored as a double rather than the more suitable int
            nr.second = 0.0; // never used. just sent to zero.
            (*scs)["nr"] = nr;
        } else {
            // found. update
            double oldN = (*foundnr).second.first;
            (*foundnr).second.first = oldN + N;
        }
        hatPos++;
    }
}

// setSeqInfo - This sets the Se= field in the comments given the information contained in the 'seqs' map.
// See SpectraSTLibEntry::getSeqInfo for the structure of the seqs map.
void SpectraSTLibEntry::setSeqInfo(map<char, map<string, pair<double, double> > *> &seqs) {

    stringstream ss;

    // number of search engines
    ss << seqs.size();

    map<char, map<string, pair<double, double> > *>::iterator se;
    for (se = seqs.begin(); se != seqs.end(); se++) {

        int N = 0;
        // N is the number of replicates, available in the "nr" field
        map<string, pair<double, double> >::iterator foundnr = se->second->find("nr");
        if (foundnr != se->second->end()) {
            N = (int) ((*foundnr).second.first + 0.5);
        }

        // search engine type, followed by num replicates IDed by this engine
        ss << '^' << (*se).first << N << ':';

        // now comes the scores
        map<string, pair<double, double> >::iterator sc;

        bool isFirst = true;
        for (sc = (*se).second->begin(); sc != (*se).second->end(); sc++) {
            if ((*sc).first == "nr") {
                // already printed N. don't print it again
                continue;
            }

            if (!isFirst) {
                ss << ',';
            }

            // calculate the mean and stdev from the values in the score map --
            // the formula are explained in SpectraSTLibEntry::getSeqInfo.
            double mean = (*sc).second.first / (double) N;
            double variance = ((*sc).second.second / (double) N - mean * mean);
            double stdev = 0.0;
            if (variance > 0.00001) stdev = sqrt(variance);
            ss.precision(4);
            if (mean > 1000.0 || mean < 0.001) {
                ss << (*sc).first << '=' << scientific << mean << '/' << scientific << stdev;
            } else {
                ss << (*sc).first << '=' << fixed << mean << '/' << fixed << stdev;
            }
            isFirst = false;
        }
    }

    setOneComment("Se", ss.str());

}

// getSampleInfo - parses out the Sample= field in the comments, 
// and incorporates the result in the 'samples' map. Note that if 'samples' already contains
// some sample information, getSampleInfo will not delete it; instead, it merges
// that existing information with what it parses out from this comment.
// (This is useful if you want to combine multiple entries -- see SpectraSTSpLibImporter.cpp)
// The structure of the 'samples' map is: 
// datasetName => (num_replicates_from_that_dataset_USED, num_replicates_from_that_dataset_TOTAL).

void SpectraSTLibEntry::getSampleInfo(map<string, pair<unsigned int, unsigned int> > &samples, string option) {

    string sampleInfo("");
    if (!getOneComment("Sample", sampleInfo)) {
        return;
    }

    // the structure of the Sample= comment is:
    // <num_different_datasets>/<dataset1>,<num_replicates_from_dataset1_USED>,<num_replicates_from_dataset1_TOTAL>/
    // <dataset2>,<num_replicates_from_dataset2_USED>,<num_replicates_from_dataset2_TOTAL>...

    string::size_type slashPos = sampleInfo.find('/', 0) + 1;

    string sa("");

    while (!((sa = nextToken(sampleInfo, slashPos, slashPos, "/\r\n")).empty())) {

        string::size_type commaPos = 0;

        string datasetName = nextToken(sa, commaPos, commaPos, ",\r\n");
        unsigned int N = atoi(nextToken(sa, commaPos + 1, commaPos, ",\r\n").c_str());
        if (option == "TOTAL_ONLY") N = 0;
        unsigned int T = atoi(nextToken(sa, commaPos + 1, commaPos, "\r\n").c_str());
        if (option == "USED_ONLY") T = 0;

        map<string, pair<unsigned int, unsigned int> >::iterator foundsa = samples.find(datasetName);
        if (foundsa == samples.end()) {
            pair<unsigned int, unsigned int> pr;
            pr.first = N;
            pr.second = T;
            samples[datasetName] = pr;
        } else {
            int oldN = (*foundsa).second.first;
            int oldT = (*foundsa).second.second;
            (*foundsa).second.first = oldN + N;
            (*foundsa).second.second = oldT + T;
        }

        slashPos++;
    }

}

// setSampleInfo - This sets the Sample= field in the comments given the information contained in the 'samples' map.
// See SpectraSTLibEntry::setSampleInfo for the structure of the 'samples' map.
void SpectraSTLibEntry::setSampleInfo(map<string, pair<unsigned int, unsigned int> > &samples) {

    stringstream ss;

    // number of different datasets
    ss << samples.size();

    map<string, pair<unsigned int, unsigned int> >::iterator sa;
    for (sa = samples.begin(); sa != samples.end(); sa++) {
        ss << '/' << (*sa).first << ',' << (*sa).second.first << ',' << (*sa).second.second;
    }

    setOneComment("Sample", ss.str());

}

// getInstrumentInfo - parses out the Instrument= field in the comments, 
// and incorporates the result in the 'instruments' map. Note that if 'instruments' already contains
// some instrument information, getInstrumentInfo will not delete it; instead, it merges
// that existing information with what it parses out from this comment.
// (This is useful if you want to combine multiple entries -- see SpectraSTSpLibImporter.cpp)
// The structure of the 'instruments' map is: 
// instrumentName => (num_replicates_from_that_instrument_USED, num_replicates_from_that_instrument_TOTAL).
void SpectraSTLibEntry::getInstrumentInfo(map<string, pair<unsigned int, unsigned int> > &instruments, string option) {

    string instrumentInfo("");
    if (!getOneComment("Inst", instrumentInfo)) {
        return;
    }

    string::size_type slashPos = instrumentInfo.find('/', 0) + 1;

    string inst("");

    while (!((inst = nextToken(instrumentInfo, slashPos, slashPos, "/\r\n")).empty())) {

        string::size_type commaPos = 0;

        string instName = nextToken(inst, commaPos, commaPos, ",\r\n");
        int N = atoi(nextToken(inst, commaPos + 1, commaPos, ",\r\n").c_str());
        if (option == "TOTAL_ONLY") N = 0;
        int T = atoi(nextToken(inst, commaPos + 1, commaPos, "\r\n").c_str());
        if (option == "USED_ONLY") T = 0;

        map<string, pair<unsigned int, unsigned int> >::iterator foundinst = instruments.find(instName);
        if (foundinst == instruments.end()) {
            pair<unsigned int, unsigned int> pr;
            pr.first = N;
            pr.second = T;
            instruments[instName] = pr;
        } else {
            unsigned int oldN = (*foundinst).second.first;
            unsigned int oldT = (*foundinst).second.second;
            (*foundinst).second.first = oldN + N;
            (*foundinst).second.second = oldT + T;
        }

        slashPos++;
    }
}

// setInstrumentInfo - This sets the Inst= field in the comments given the information contained in the 'samples' map.
// See SpectraSTLibEntry::setInstumentInfo for the structure of the 'instruments' map.
void SpectraSTLibEntry::setInstrumentInfo(map<string, pair<unsigned int, unsigned int> > &instruments) {

    stringstream ss;

    ss << instruments.size();

    map<string, pair<unsigned int, unsigned int> >::iterator inst;
    for (inst = instruments.begin(); inst != instruments.end(); inst++) {
        ss << '/' << (*inst).first << ',' << (*inst).second.first << ',' << (*inst).second.second;
    }

    setOneComment("Inst", ss.str());
}

// annotatePeaks - just calls m_peakList->annotate(), then add a comment of fracUnassigned
void SpectraSTLibEntry::annotatePeaks(bool redo, bool fixMz) {
    m_peakList->annotate(redo, fixMz);

    string unassignedStr(m_peakList->getFracUnassignedStr());
    if (!(unassignedStr.empty())) {
        setOneComment("FracUnassigned", m_peakList->getFracUnassignedStr());
    }
}

// getNTT - returns the number of tryptic termini
unsigned int SpectraSTLibEntry::getNTT() {
    if (m_pep) {
        return (m_pep->NTT());
    }
    string nttStr("");
    if (getOneComment("NTT", nttStr)) {
        return (atoi(nttStr.c_str()));
    }
    return (0);

}

// getNMC - returns the number of internal missed cleavage
unsigned int SpectraSTLibEntry::getNMC() {
    if (m_pep) {
        return (m_pep->NMC());
    }
    string nmcStr("");
    if (getOneComment("NMC", nmcStr)) {
        return (atoi(nmcStr.c_str()));
    }
    return (0);

}

// freePeakList - deletes the peak list (usually for memory saving)
void SpectraSTLibEntry::freePeakList() {

    if (m_peakList) {
        delete (m_peakList);
    }
    m_peakList = NULL;
}

// makeDecoy - makes a decoy entry. Takes the peptide ion ID, shuffles the sequence, then repositions the peaks in the peak list
// based on the shuffled sequence.
void SpectraSTLibEntry::makeDecoy(Peptide *decoyPep, unsigned int index) {

//  if (!m_pep) return;

    // shuffles the peptide ion ID sequence

//  Peptide* decoyPep = new Peptide(m_pep->removeModOfType("Label:13C(6)"), m_charge);
//  m_mw = decoyPep->averageMH();
//  m_precursorMz = decoyPep->averageMZ();
//  stringstream precursorMzss;
//  precursorMzss.precision(2);
//  precursorMzss << fixed << m_precursorMz;
//  setOneComment("Parent", precursorMzss.str());

//  Peptide* decoyPep = new Peptide(m_pep->shufflePeptideSequence(), m_charge);
//  decoyPep->prevAA = m_pep->prevAA;
//  decoyPep->nextAA = m_pep->nextAA;

    // modify the peak list
    m_peakList->annotate();
    m_peakList->setPeptidePtr(decoyPep);


    // deliberately not change the m_mw and m_precursorMz fields -- to maintain the precursor m/z distribution
    // even when the decoyPep has a different mass as the origPep. (This happens when the shuffle doesn't produce // sufficiently different sequences and two random AA have to be added.)

    // SHUFFLE
    m_peakList->repositionPeaks();
    // SHIFT
    //  m_peakList->shiftAllPeaks(20);

    // rename the entry, and add comments
    m_name = decoyPep->interactStyleWithCharge();
    m_fullName = decoyPep->interactStyleFullWithCharge();
    if (!(m_fragType.empty())) {
        m_fullName += " (" + m_fragType + ")";
    }

    setOneComment("OrigPeptide", m_pep->interactStyleFullWithCharge());
    setOneComment("Mods", decoyPep->mspMods());
    unsigned int ntt = decoyPep->NTT();
    stringstream nttss;
    nttss << ntt;
    setOneComment("NTT", nttss.str());
    if (ntt == 2) {
        setOneComment("Pep", "Tryptic");
    } else if (ntt == 1) {
        setOneComment("Pep", "Semi-tryp");
    } else {
        setOneComment("Pep", "Non-tryp");
    }

    stringstream nmcss;
    nmcss << decoyPep->NMC();
    setOneComment("NMC", nmcss.str());

    // setOneComment("Remark", "Artificial_Light");

    labelAsDecoy(index);


    if (m_pep) delete (m_pep);
    m_pep = decoyPep;


}

void SpectraSTLibEntry::labelAsDecoy(unsigned int index) {

    stringstream remarkss;
    remarkss << "DECOY_" << index;
    setOneComment("Remark", remarkss.str());
    setOneComment("Spec", "Decoy");

    string protein("");
    if (getOneComment("Protein", protein)) {
        string::size_type slashPos = protein.find('/', 0);
        if (slashPos != string::npos && slashPos < protein.length() - 1) {
            string nextPro("");
            stringstream dpss;
            unsigned int numPro = atoi(protein.substr(0, slashPos).c_str());
            dpss << numPro;
            while (!((nextPro = nextToken(protein, slashPos + 1, slashPos, "/\t\r\n")).empty())) {
                dpss << '/' << remarkss.str() << '_' << nextPro;
            }
            setOneComment("Protein", dpss.str());
        } else {
            setOneComment("Protein", remarkss.str() + "_" + protein);
        }
    } else {
        setOneComment("Protein", remarkss.str() + "_" + "Unknown");
    }

}

// makeSemiempiricalSpectrum - makes a semi-empirical entry. Takes the new peptide ID, then repositions the peaks 
// in the current peak list based on new ID. Note that ownership of newPep passes to SpectraSTLibEntry.
void SpectraSTLibEntry::makeSemiempiricalSpectrum(Peptide *newPep) {

    setOneComment("OrigPeptide", m_pep->interactStyleFullWithCharge());

    if (m_pep) delete (m_pep);
    m_pep = newPep;

    synchWithPep();

    m_peakList->annotate();

    m_peakList->setPeptidePtr(m_pep);

    m_peakList->repositionPeaks(false);

    setOneComment("Spec", "Semi-empirical");

    // TODO - may have to change the protein name...


}

// comparator used by SpectraSTSpLibImporter
bool SpectraSTLibEntry::sortEntryPtrsByPeakListWeightDesc(SpectraSTLibEntry *a, SpectraSTLibEntry *b) {

    return (a->getPeakList()->getWeight() > b->getPeakList()->getWeight());

}

// evaluatePhosphoSiteAssignment - experimental, try to assess how good the site assignments are by only looking at 
// those peaks that can distinguish the current assignment and the closest alternatives
void SpectraSTLibEntry::evaluatePhosphoSiteAssignment() {

    if (!m_pep || !m_peakList || !m_pep->isModsSet) {
        return;
    }

    vector<int> nophos; // positions of all non-phosphorylated S, T, Y
    vector<int> phos; // positions of all phosphorylated S, T, Y

    for (unsigned int pos = 0; pos < m_pep->NAA(); pos++) {
        char aa = m_pep->stripped[pos];
        if (aa == 'S' || aa == 'T' || aa == 'Y') {
            if (m_pep->mods.find(pos) != m_pep->mods.end() && m_pep->mods[pos] == "Phospho") {
                phos.push_back(pos);
            } else {
                nophos.push_back(pos);
            }
        }
    }


    if (phos.empty() || nophos.empty()) {
        if (!phos.empty()) setOneComment("PSite", "Good"); // no ambiguity
        return;
    }

    //  cerr << "PEPTIDE: " << m_pep->interactStyleWithCharge() << endl;

    stringstream sitess;
    sitess << phos.size();

    bool isGood = true;

    double phosMass = Peptide::getModMonoisotopicMass("Phospho");

    // find closest non-phosphorylated sites on each side; will only consider these as alternatives
    for (vector<int>::size_type p = 0; p < phos.size(); p++) {
        int left = -1;
        int right = -1;
        for (vector<int>::size_type s = 0; s < nophos.size(); s++) {
            if (nophos[s] < phos[p]) {
                left = nophos[s];
            } else if (nophos[s] > phos[p]) {
                right = nophos[s];
                break;
            }
        }

        // cerr << "  Site=" << phos[p] << " (Left=" << left << ", Right=" << right << ")" << endl;

        double leastConfidentRatio = 10000000.0;

        if (left != -1) {
            // potential site on the left

            // check b ions from left+1 to pos, y ions from len-pos to len-left-1
            double forEvidence = 0.001;
            double againstEvidence = 0.001;

            for (int i = left + 1; i <= phos[p]; i++) {
                for (unsigned int ch = 1; ch <= (unsigned int) m_pep->charge; ch++) {
                    double mz = m_pep->monoisotopicMZFragment('b', i, ch);
                    double foundInt = m_peakList->findPeak(mz, 0.5);
                    forEvidence += foundInt;
                    // cerr << "    LEFT(+): " << 'b' << i << "^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('y', m_pep->NAA() - i, ch);
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    forEvidence += foundInt;
                    // cerr << "    LEFT(+): " << 'y' << m_pep->NAA() - i << "^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('b', i, ch) + phosMass / (double) ch;
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    againstEvidence += foundInt;
                    // cerr << "    LEFT(-): " << 'b' << i << "+80^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('y', m_pep->NAA() - i, ch) - phosMass / (double) ch;
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    againstEvidence += foundInt;
                    // cerr << "    LEFT(-): " << 'y' << m_pep->NAA() - i << "-80^" << ch << " = " << foundInt << endl;

                }
            }

            double ratio = forEvidence / againstEvidence;
            // cerr << "    LEFT_RATIO = " << ratio << endl;
            if (ratio < leastConfidentRatio) leastConfidentRatio = ratio;
        }

        if (right != -1) {
            // potential site on the right

            // check b ions from left+1 to pos, y ions from len-pos to len-left-1
            double forEvidence = 0.001;
            double againstEvidence = 0.001;

            for (int i = phos[p] + 1; i <= right; i++) {
                for (unsigned int ch = 1; ch <= (unsigned int) m_pep->charge; ch++) {
                    double mz = m_pep->monoisotopicMZFragment('b', i, ch);
                    double foundInt = m_peakList->findPeak(mz, 0.5);
                    forEvidence += foundInt;
                    // cerr << "    RIGHT(+): " << 'b' << i << "^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('y', m_pep->NAA() - i, ch);
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    forEvidence += foundInt;
                    // cerr << "    RIGHT(+): " << 'y' << m_pep->NAA() - i << "^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('b', i, ch) - phosMass / (double) ch;
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    againstEvidence += foundInt;
                    // cerr << "    RIGHT(-): " << 'b' << i << "-80^" << ch << " = " << foundInt << endl;
                    mz = m_pep->monoisotopicMZFragment('y', m_pep->NAA() - i, ch) + phosMass / (double) ch;
                    foundInt = m_peakList->findPeak(mz, 0.5);
                    againstEvidence += foundInt;
                    // cerr << "    RIGHT(-): " << 'y' << m_pep->NAA() - i << "+80^" << ch << " = " << foundInt << endl;

                }
            }

            double ratio = forEvidence / againstEvidence;
            // cerr << "  RIGHT_RATIO = " << ratio << endl;
            if (ratio < leastConfidentRatio) leastConfidentRatio = ratio;
        }

        // print to Comment
        sitess << '/' << leastConfidentRatio;

        // if any of the assigned site is not overwhelmingly supported, then the whole phosphopeptide is labeled "Uncertain"
        if (leastConfidentRatio < 10.0) isGood = false;

    }

    setOneComment("PSiteConf", sitess.str());
    if (isGood) {
        setOneComment("PSite", "Good");
    } else {
        setOneComment("PSite", "Uncertain");
    }


}


int SpectraSTLibEntry::getMassDiffInt() {

    int massDiffInt = 0;
    string mdStr("");
    if (getOneComment("MassDiff", mdStr)) {
        double massdiff = atof(mdStr.c_str());
        if (massdiff >= 0) {
            massDiffInt = (int) (0.5 + massdiff);
        } else {
            massDiffInt = (int) (-0.5 + massdiff);
        }
    }
    return (massDiffInt);
}

void SpectraSTLibEntry::recordDot(int scanNum, double dot) {

    if (!m_dotTimeProfile) m_dotTimeProfile = new vector<pair<int, float> >;

    pair<int, float> p;
    p.first = scanNum;
    p.second = (float) dot;
    m_dotTimeProfile->push_back(p);
}

SpectraSTPeakList *SpectraSTLibEntry::generateMS1(SpectraSTSearchParams &params) {

    if (m_ms1) return (m_ms1); // won't do it again

    m_ms1 = new SpectraSTPeakList(0.0, 0, 4, true, "CID-QTOF");

    double precursorMass = m_precursorMz * (double) m_charge;
    double logPrecursorMass = log(precursorMass);

    double first = -0.34 * logPrecursorMass + 2.9;
    double second = -6.4e-8 * precursorMass * precursorMass + 2.6e-4 * precursorMass + 9.8e-2;
    double third = 0.12 * logPrecursorMass - 0.74;
    double fourth = 1.0 - first - second - third;

    m_ms1->insertForSearch(m_precursorMz, first, "");
    m_ms1->insertForSearch(m_precursorMz + 1.0 / (double) m_charge, second, "");
    m_ms1->insertForSearch(m_precursorMz + 2.0 / (double) m_charge, third, "");
    m_ms1->insertForSearch(m_precursorMz + 3.0 / (double) m_charge, fourth, "");

    m_ms1->prepareForSearch(params, true);

    return (m_ms1);
}


void SpectraSTLibEntry::printDotTimeProfile(ofstream &fout) {

    if (!m_dotTimeProfile) return;

    fout << m_libId << '\t' << m_name << '\t' << getFirstProtein() << '\t' << getPrecursorMz();

    float sum = 0.0;
    float sumSq = 0.0;
    int count = 0;

    for (vector<pair<int, float> >::iterator i = m_dotTimeProfile->begin(); i != m_dotTimeProfile->end(); i++) {
        //  if (*i > 0.01) {
        sum += i->second;
        sumSq += i->second * i->second;
        count++;
        //  }
    }

    float average = sum / (float) count;
    float stdev = sqrt(sumSq / (float) (count - 1) - average * average);
    bool good = false;

    ofstream sout;
    myFileOpen(sout, "tmp.spec");

    int scan = 0;
    float max = 0.0;
    float lastValue = 0.0;
    int direction = 1;

    vector<pair<int, float> > peaks;

    float peakPickThreshold = average * 5.0;

    sum = 0.0;
    sumSq = 0.0;
    count = 0;

    for (vector<pair<int, float> >::iterator i = m_dotTimeProfile->begin(); i != m_dotTimeProfile->end(); i++) {
        float value = i->second - average;
        if (value < 0.0) value = 0.0;

        if (value > 0.01) {
            count++;
            sum += value;
            sumSq += value * value;
        }

        if (value > peakPickThreshold) {
            if (lastValue < value) {
                direction = 1;
                lastValue = value;
            } else {
                if (direction > 0) {
                    // change direction, reach peak
                    pair<int, float> pr;
                    pr.first = i->first;
                    pr.second = lastValue;
                    peaks.push_back(pr);
                }
                direction = -1;
                lastValue = value;
            }
        } else {
            if (direction > 0 && lastValue > peakPickThreshold) {
                // change direction, reach peak
                pair<int, float> pr;
                pr.first = i->first;
                pr.second = lastValue;
                peaks.push_back(pr);
            }
            lastValue = 0.0;
        }

        if (max < value) max = value;
        sout.precision(3);
        sout << fixed << i->first << '\t';
        sout.precision(3);
        sout << fixed << value << endl;
    }

    sout.close();

    float top = peakPickThreshold;
    float sec = peakPickThreshold;
    float peakAverage = 0.0;
    int peakCount = 0;

    for (vector<pair<int, float> >::iterator pk = peaks.begin(); pk != peaks.end(); pk++) {

        peakCount++;
        peakAverage += pk->second;

        fout.precision(1);
        fout << '\t' << pk->first << ',' << fixed << pk->second;
        if (pk->second > top) {
            sec = top;
            top = pk->second;
        } else if (pk->second > sec) {
            sec = pk->second;
        }
    }

    peakCount--;
    peakAverage -= top;
    peakAverage /= (float) peakCount;

    float gap = (top - sec) / top;

    // sum -= max;
    // sumSq -= max * max;
    // count--;
    // average = sum / (float)count;
    // stdev = sqrt(sumSq / (float)(count - 1) - average * average);
    //float zscore = (max - average) / stdev;

    fout.precision(3);
    fout << gap;
    fout.precision(3);
    fout << '\t' << fixed << max / average << '\t' << fixed << average << '\t' << fixed << max << endl;
//  fout.precision(3);
//  fout << fixed << zscore << endl;

    if (max >= 500.0) {
//  if (zscore >= 10.0) {
        ofstream gout;
        myFileOpen(gout, "tmp.gp");
        gout << "set terminal png" << endl;
        gout << "set output \"" << getSafeName() << ".png\"" << endl;
        gout.precision(2);
        gout << "set label \"" << fixed << getPrecursorMz() << "\" at 5," << max * 0.4 << endl;
        gout.precision(1);
        gout << "set label \"" << fixed << max << '/' << fixed << average << "=" << fixed << max / average << "\" at 5,"
             << max * 0.6 << endl;
        // gout << "set label \"(" << fixed << max << '-' << fixed << average << ")/" << fixed << stdev << '=' << fixed << zscore << "\" at 5," << max * 0.6 << endl;
        gout.precision(2);
        gout << "set label \"" << fixed << gap << "\" at 5," << max * 0.8 << endl;
        gout << "plot \"tmp.spec\" using 1:2 with lines title \"" << getSafeName() << "  " << getFirstProtein() << "\""
             << endl;
        gout.close();
        myGnuplot("tmp.gp");

    }


}


/*
void SpectraSTLibEntry::recordDot(double dot) {
 
  int bin = (int)(dot * 50.0);
  
  if (bin < 0) bin = 0;
  if (bin > 49) bin = 49;
  
  m_dotHistogram[bin]++;
  
}

void SpectraSTLibEntry::printDotHistogram(ofstream& fout) {
 
  fout << m_libId << '\t' << m_peakList->getNumPeaks();
  
  for (vector<unsigned int>::iterator i = m_dotHistogram.begin(); i != m_dotHistogram.end(); i++) {
    fout << '\t' << (*i);
  }
    
  fout << endl;
}
*/
