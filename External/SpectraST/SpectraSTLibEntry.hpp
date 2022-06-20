#ifndef SPECTRASTLIBENTRY_HPP_
#define SPECTRASTLIBENTRY_HPP_

#include "SpectraSTPeakList.hpp"
#include "Peptide.hpp"
#include <iostream>
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

/* Class: SpectraSTLibEntry
 * 
 * The class that represents a library entry. 
 * 
 */

using namespace std;

class SpectraSTLibEntry {

public:

    // only two ways to construct an entry object:
    // 1. by giving all required information as arguments
    // 2. by giving a ifstream object pointing to the beginning of an entry in a library file
    SpectraSTLibEntry(Peptide *pep, string comments, string status, SpectraSTPeakList *peakList, string fragType = "");

    SpectraSTLibEntry(string name, double precursorMz, string comments, string status, SpectraSTPeakList *peakList,
                      string fragType = "");

    SpectraSTLibEntry(ifstream &libFin, bool binary = false, bool shortAnnotation = false);

    // copy constructor and assignment operator
    SpectraSTLibEntry(SpectraSTLibEntry &other);

    SpectraSTLibEntry &operator=(SpectraSTLibEntry &other);

    // destructor
    ~SpectraSTLibEntry();

    // Accessor methods
    unsigned int getLibId() { return m_libId; }

    fstream::off_type getLibFileOffset() { return m_libFileOffset; }

    string getName() { return m_name; }

    int getCharge() { return (m_pep ? m_pep->charge : m_charge); }

    string getMods() { return (m_pep ? m_pep->mspMods() : ""); }

    string getStatus() { return m_status; }

    string getFullName() { return (m_fullName); }

    string getSafeName(); // a name that can be used as part of a file name (i.e. no '.' or '/' or [], etc)
    double getMw() { return m_mw; }

    double getPrecursorMz() { return m_precursorMz; }

    double getAveragePrecursorMz();

    string getCommentsStr();

    Peptide *getPeptidePtr() { return m_pep; }

    SpectraSTPeakList *getPeakList() { return m_peakList; }

    double getProb(double valueIfNotFound = 0.0);

    unsigned int getNrepsUsed(unsigned int valueIfNotFound = 1);

    string getFragType() { return (m_fragType); }

    int getMassDiffInt();

    // Setters
    void setStatus(string status) { m_status = status; }

    void setLibId(unsigned int libId) { m_libId = libId; }

    void setLibFileOffset(fstream::off_type libFileOffset) { m_libFileOffset = libFileOffset; }

    void setPeakList(SpectraSTPeakList *peakList);

    void
    synchWithPep(); // sets m_name, m_fullName, m_mw, m_precursorMz and Mods= in Comments according to the m_pep object
    void setPrecursor(double precursorMz, int charge = 0); // charge = 0 means the charge is unknown
    void setFragType(string &fragType);

    // Predicates for checking the modification on cysteine (if any)
    bool isCleavableICAT() { return (m_pep ? m_pep->isCleavableICAT() : false); }

    bool isUncleavableICAT() { return (m_pep ? m_pep->isUncleavableICAT() : false); }

    bool isCAMCysteine() { return (m_pep ? m_pep->isCAMCysteine() : false); }

    bool hasUnmodifiedCysteine() { return (m_pep ? m_pep->hasUnmodifiedCysteine() : false); }

    unsigned int getNTT();

    unsigned int getNMC();

    // Accessing the comment
    bool getOneComment(string attr, string &value);

    bool setOneComment(string attr, string value);

    bool deleteOneComment(string attr);

    // Accessing specific fields in comment
    string getRawSpectra();

    string getBestRawSpectrum();

    void getAllProteins(vector<string> &proteins);

    string getFirstProtein();

    string getFirstProtein(int &proteinCount);

    // Accessing and setting the Se=, Sample= and Inst= fields in the comment
    void getSeqInfo(map<char, map<string, pair<double, double> > *> &seqs);

    void setSeqInfo(map<char, map<string, pair<double, double> > *> &seqs);

    void getSampleInfo(map<string, pair<unsigned int, unsigned int> > &samples, string option = "");

    void setSampleInfo(map<string, pair<unsigned int, unsigned int> > &samples);

    void getInstrumentInfo(map<string, pair<unsigned int, unsigned int> > &instruments, string option = "");

    void setInstrumentInfo(map<string, pair<unsigned int, unsigned int> > &instruments);

    // Annotate peaks of the peakList
    void annotatePeaks(bool redo = false, bool fixMz = false);

    // Evaluate phosphosite assignment
    void evaluatePhosphoSiteAssignment();

    // File I/O methods
    void writeToFile(ofstream &libFout);

    void writeToBinaryFile(ofstream &libFout);

    void writeDtaFile(string dtaFileName);

    void writeMRM(ofstream &mrmFout, string format);

    void writeInfo(ofstream &mrmFout);

    void writeMgfFile(ofstream &mgfFout);

    void writePAIdent(ofstream &fout, string baseName);

    // Deleting the peak list
    void freePeakList();

    // Making a decoy spectrum
    void makeDecoy(Peptide *decoyPep, unsigned int index = 0);

    void labelAsDecoy(unsigned int index = 0);

    // Making a semiempirical spectrum
    void makeSemiempiricalSpectrum(Peptide *newPep);

    // comparator for sorting by the weight of the peak list.
    static bool sortEntryPtrsByPeakListWeightDesc(SpectraSTLibEntry *a, SpectraSTLibEntry *b);

    // testing functions to record search results for plotting null model
    // void recordDot(double dot);
    // void printDotHistogram(ofstream& fout);
    void recordDot(int scanNum, double dot);

    vector<pair<int, float> > *getDotTimeProfile() { return (m_dotTimeProfile); }

    void clearDotTimeProfile() { if (m_dotTimeProfile) m_dotTimeProfile->clear(); }

    void printDotTimeProfile(ofstream &fout);

    SpectraSTPeakList *generateMS1(SpectraSTSearchParams &params);

private:

    // fields

    Peptide *m_pep; // defines the peptide ion.
    unsigned int m_libId; // just a zero-based unique identifier for each entry
    fstream::off_type m_libFileOffset; // file offset for random access
    string m_name; // m_name is an interact-style (with charge) string representing the peptide ion
    double m_mw;
    double m_precursorMz;
    string m_status;
    string m_fullName; // m_fullName is like m_name but with prevAA and nextAA
    string m_commentsStr;
    int m_charge;
    string m_fragType;

    // m_comments is an object representation of the Comment: field.
    // this is however not used (for efficiency sake) for simple accessing of a comment field
    // which will be done by string parsing of m_commentsStr instead.
    map<string, string> *m_comments;

    // m_peakList IS the property of SpectraSTLibEntry!
    SpectraSTPeakList *m_peakList;

    // vector<unsigned int> m_dotHistogram;
    vector<pair<int, float> > *m_dotTimeProfile;

    // m_ms1
    SpectraSTPeakList *m_ms1;

    // parse m_commentsStr to create m_comments
    void parseCommentsStr();

    // reading from files - these are private, so to read from files the constructor has to be called
    void readFromFile(ifstream &libFin, bool forSearch);

    void readFromBinaryFile(ifstream &libFin, bool forSearch);


};

#endif /*SPECTRASTLIBENTRY_HPP_*/
