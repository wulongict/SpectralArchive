//
// Created by wulong on 6/19/17.
//

#ifndef MYTOOL_XMLFILEPARSER_H
#define MYTOOL_XMLFILEPARSER_H

#include "../External/rapidxml-1.13/rapidxml.hpp"
#include <string>
#include <vector>

#include <zlib.h>
#include <map>
#include <numeric>
#include <memory>
#include <sstream>
#include "ICGtInfoUpdate.h"

using namespace std;

const double UNKOWN_PEPTIDEPROPHET_SCORE = 0;
// todo: put all of those values into one header or namespace.

class PeakList;
class SPsmAnnotation;
class modification;

using namespace rapidxml;

namespace XMLParser{

    void printAttributes(xml_node<> *node);
    void print(xml_node<> *node);
    void printChildren(xml_node<> *node);

    xml_node<> *breadth_first_search(string node_name, xml_node<> *root, bool verbose);
    vector<xml_node<>*> find_all_siblings(string node_name, xml_node<> * first_node);
    template<typename T> void first_attribute(xml_node<> * node, string attr_name, T & attr_val);
};
class SearchHit{
public:
    long getMemSize(){
        return sizeof(SearchHit)
        +  m_peptide.length()
        + m_preAA.length()
         +m_nextAA.length()
        + m_protein.length()
         +m_modified_peptide.length()
         +accumulate(m_alternative_proteins.begin(),m_alternative_proteins.end(),0L,
                     [](long b, const string& a){
             return b + a.length();
         })+(sizeof(int)+sizeof(double))*m_pos_mass.size();

    }
    typedef float real;

    real m_xcorr;
    real m_deltacn;
    real m_deltacnstar;
    real m_spscore;
    real m_sprank;
    real m_expect;
    real m_hyperscore;
    real m_nextscore;
    real m_massdiff;
    real m_peptideprophet_prob;
    real m_iprophet_prob;
    real m_nterm_mass;
    real m_cterm_mass;
    int m_hit_rank;
    int m_num_matched_ion;
    int m_tot_num_ions;
    int m_num_missed_cleavages;
    string m_peptide;
    string m_preAA;
    string m_nextAA;
    string m_protein;
    string m_modified_peptide;
    vector<string> m_alternative_proteins;
    vector<string> m_localizations;
    map<int, double> m_pos_mass;
//    double m_fvalue; todo: in the future I can add this one
public:

    void print();
    explicit SearchHit(xml_node<>* node);
    SearchHit(int hit_rank, string peptide, string preAA, string nextAA, string protein, int num_matched_ion,
                  int tot_num_ions, int num_missed_cleavages, string modified_peptide, double xcorr,
                  double deltacn, double deltacnstar, double spscore, double sprank, double expect,
                  double massdiff, double peptideprophet_prob, double iprophet_prob, double nterm_mass,
                  double cterm_mass, map<int, double> pos_mass, vector<string> &altprots);
    SearchHit(const SearchHit & other);
    string getFirstAlterProteinNotDecoy();
    bool isDecoy(bool useAlternativeProt);
    string getProtein_UseAlterProteinIfItsNotDecoy(bool useAlternativeProt=true);

};

class PSMInfo{
public:
    long getMemSize(){
        return sizeof(PSMInfo)
             +spectrum.length()+m_basename.length()
               +accumulate(searchhits.begin(),searchhits.end(),0L,
                           [](long b, const shared_ptr<SearchHit>& a){
                               return b + a->getMemSize();
                           });

    }
    double retention_time;
    double precursorNeutMass;
    int start_scan;
    int end_scan;
    int index; // index in pepxml
    int charge;
    string spectrum;
    string m_basename;
    vector<shared_ptr<SearchHit> > searchhits;
public:
    void print();
    explicit PSMInfo(xml_node<>* spectrum_query_node);
    PSMInfo(const PSMInfo & other);
    PSMInfo();
    PSMInfo(string basenameStr, int scan, int chg=-1):PSMInfo(){
        start_scan = scan;
        end_scan = scan;
        charge = chg;
        m_basename = basenameStr;
    }
    bool isDecoy(bool useAlternativeProt, int hitrank);
    void printSearchHit(int k);
    string getProtein_UseAlterProteinIfItsNotDecoy(bool useAlternativeProt, int hitrank);

    string getFirstAlterProteinNotDecoy();
    double getParentMz();
    string toTsvString(){
        // the header line
        string header="spectrum\tstart_scan\tindex\tprecursormass\tcharge\tRT\tpeptide\tmodified_peptide\tpProb\tiProb\tmassDiff\tprotein\tlocalization";
        ostringstream  oss;
        oss << spectrum << "\t"
         << start_scan << "\t"
         << index << "\t"
         << precursorNeutMass << "\t"
         << charge << "\t"
         << retention_time << "\t";
        if(not searchhits.empty()){
            oss << searchhits[0]->m_peptide << "\t"
            << searchhits[0]->m_modified_peptide << "\t"
            << searchhits[0]->m_peptideprophet_prob << "\t"
            << searchhits[0]->m_iprophet_prob << "\t"
            << searchhits[0]->m_massdiff << "\t"
            << searchhits[0]->m_protein << "\t";
            if(not searchhits[0]->m_localizations.empty()){
                for(const auto & each: searchhits[0]->m_localizations){
                    oss << each << "\t";
                }

            }
        }
        return oss.str();
    }
};

class Spectrum{
public:
    int m_ms_level;
    double m_RT_in_seconds;
    PeakList *pkl;
    double m_precursormz;
    int m_precursorcharge;
    int spectrum_scan;
    string spectrum_name;
public:
    Spectrum(const Spectrum &other);
    explicit Spectrum(xml_node<>* spec_node);
    ~Spectrum();
    void parse_ms_level(const xml_node<> *spec_node);
    void parse_retention_time(const xml_node<> *spec_node);
    void parse_peak_list(const xml_node<> *spec_node, const int peaknum);
    void parse_binary_list_param(const xml_node<> *x, bool &zlib, string &valuename, int &bits);
    void parse_precursorinfo( xml_node<> *spec_node);
};

class ICPepXMLParser{
public:
    ICPepXMLParser() =default;
    virtual ~ICPepXMLParser() = default;

    // get PSM information.
    virtual bool getPsmInfo(PSMInfo &psminfo)=0;
    virtual int getPsmNum() = 0;
    virtual void getPSMInfobyindex(int i , PSMInfo &psminfo) = 0;
    virtual pair<int, int> getIndexRange(string filename) = 0;
    virtual vector<string> getAllSoruceFiles()= 0 ;

};

class CometPepXMLParser: public ICGtInfoUpdate, public ICPepXMLParser
{
    string m_filename;
//    char * m_buf;
    xml_node<>* m_currentNode;
//    xml_document<> doc;
    vector<shared_ptr<PSMInfo>> psm;
    //vector<PSMInfo> psm;
    map<int,vector<int>> m_scan2idxvec;
    multimap<int, int> m_scan2psminfoidx;
    map<string, int> m_spectrumName2psminfoidx;
    bool m_use_scan2idxvec;

    void export_psm_info(vector<shared_ptr<PSMInfo>> & psm,xml_document<> &doc);
public:
    map<string, pair<int, int>> m_filename2indeRange;
    vector<string> m_allSourceFiles;
    bool getPsmInfo(PSMInfo &psminfo ) override;
    int getPsmNum()override;
    pair<int, int> getIndexRange(string filename) override;
    vector<string> getAllSoruceFiles() override;
    bool getPSMInfobyScan(int scan, PSMInfo &psminfo);
    bool getPSMInfobyScanCharge(int scan, int charge, PSMInfo &psminfo);
    void getPSMInfobyindex(int i , PSMInfo &psminfo);

    explicit CometPepXMLParser(string filename);
    ~CometPepXMLParser() ;

    void exportToTXT(string filename);
    bool updateGtInfo(SPsmAnnotation &gtinfo) override;
    void getscoreandpep_mod( int scan, string &peptide, double &searchscore,
                            modification & modinfo, double &parentMZ, int  & charge);
};

// Support PeptideProphet or iProphet results
class PeptideProphetParser:  public ICGtInfoUpdate, public ICPepXMLParser
{
    string m_filename;
    char * m_buf;
    xml_node<> * m_currentNode;
    xml_document<> doc;
    vector<PSMInfo> m_psm;
    multimap<int, int> m_scan2psminfoidx;
    map<string, int> m_spectrumName2psminfoidx;
    double m_threshold;
    bool m_use_iProb;
    void export_psm_info(vector<PSMInfo> & psm);
    map<string, pair<int, int>> m_filename2indeRange;

public:
    vector<string> m_allSourceFiles;
    explicit PeptideProphetParser(string &filename);
    void initialize();
    bool isPSMSignificant(int i);
    bool isPSMSignificant(PSMInfo &psminfo );
    PeptideProphetParser(const PeptideProphetParser & other);
    ~PeptideProphetParser() ;
    string getInputfile() const;
public:
    int getPsmNum();
    vector<string> getAllSoruceFiles() override;
    pair<int, int> getIndexRange(string filename) override;
    bool getPsmInfo(PSMInfo & psminfo) override;

    void getPSMInfobyindex(int i , PSMInfo &psminfo);
    bool getPSMInfobyScan(int scan, PSMInfo & psminfo);  // TO BE DELETED
    bool getPSMInfobyScanFileName(string filename, int scan, PSMInfo &psminfo);
    bool getAllPSMsWithScanFileName(string filename, int scan, vector<PSMInfo> &psms);

    bool updateGtInfo(SPsmAnnotation &gtinfo) override;
    void filter_with_FDR(double fdr_threshold, vector<PSMInfo> &newpsm);
    double getThresholdForFDR(double fdr_threshold, bool use_iProb);
private:
    // CALL different charge states
    bool getPSMInfo(int scan, PSMInfo &psminfo, string &basename);
    bool getPSMInfobySpectrumName(string spectrumName, PSMInfo &psminfo);
    // CALL using spectrumName
    bool getPSMInfo(int scan, int chg, PSMInfo &psminfo, string basename);
    bool updateGtInfo(string name, int scan, int chg, SPsmAnnotation &gtinfo);
    bool updateGtInfoOnSpectrumName(string spectrumName, SPsmAnnotation &gtinfo );
    void filterWithThreshold(double min_probability_score, vector<PSMInfo> &newpsm);
    string generateSpectrumName(int scan, int chg, const string &basename) const;

};

ICPepXMLParser *createPepXML(string filename);

// todo write a mzML reader here
class mzMLReader
{
    string m_filename;
    char *m_buf;
    xml_node<> * m_currentNode;
    xml_document<> doc;
    vector<Spectrum > Spectra;
    void get_all_spectra();
public:
    explicit mzMLReader(string filename);
    ~mzMLReader();
    vector<Spectrum>* getSpectra(){
            return & Spectra;
    }
};

double get_peptidePropeht_prob(SearchHit &sh);
double get_iprophet_prob(SearchHit &sh);



void extract_pep_prob(vector<double> & tProbs, vector<double>& dProbs, PeptideProphetParser &ppp, double (*getscore)(SearchHit &hit), bool useAlternativeProt=true);



#endif //MYTOOL_XMLFILEPARSER_H
