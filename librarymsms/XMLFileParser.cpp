//
// Created by wulong on 6/19/17.
//

#include <zlib.h>
#include <iostream>
#include <string>
#include <queue>
#include <cstring>
#include <spdlog/spdlog.h>
#include <fstream>
#include <set>

#include "XMLFileParser.h"
#include "Util.h"
#include "CPSMAnnotation.h"
#include "PeakList.h"
#include "base64github.h"
#include "../External/rapidxml-1.13/rapidxml.hpp"
using namespace rapidxml;
using namespace std;
using namespace XMLParser;


void printXYZ(int length, char *buffer) {
    for (int i = 0; i < length; ++i) {
        cout << "buffer: " << buffer[i];
    }
}

xml_node<> *XMLParser::breadth_first_search(string node_name, xml_node<> *root, bool verbose)  {
    xml_node<> *result = nullptr;
    queue<xml_node<> *> queue_nodes;
    queue_nodes.push(root);

    while (not queue_nodes.empty())  {
        xml_node<> *node = queue_nodes.front();
        if (strcmp(node->name(), node_name.c_str()) == 0) {
            result = node;
            break;
        }
        for (xml_node<> *child = node->first_node(); child; child = child->next_sibling()) {
            queue_nodes.push(child);
        }
        queue_nodes.pop();
    }
    return result;
}

vector<xml_node<> *> XMLParser::find_all_siblings(string node_name, xml_node<> *first_node) {
    vector<xml_node<> *> siblings;
    for (; first_node; first_node = first_node->next_sibling(node_name.c_str())) {
        siblings.push_back(first_node);
    }
    return siblings;
}

void XMLParser::printAttributes(xml_node<> *node) {
    cout << "[Info] Attribute -------------" << endl;
    for (xml_attribute<> *attr = node->first_attribute();
         attr; attr = attr->next_attribute()) {
        cout << " " << attr->name() << " = " << attr->value() << "\n";
    }
    cout << "[Info] End ------- " << endl;
}

void XMLParser::print(xml_node<> *node) {
    cout << "Node name: " << node->name() << endl;
    printAttributes(node);
}

void XMLParser::printChildren(xml_node<> *node) {
    int childrenNum = 0;
    cout << "[Info] Children --------------" << endl;
    for (xml_node<> *newnode = node; newnode; newnode = newnode->next_sibling()) {
        childrenNum++;
        cout << "Child Num: " << childrenNum << endl;
        print(newnode);
        cout << "[Info] End -------" << endl;
        if (childrenNum > 5) {
            return;;
        }
    }
}

template<typename T>
void XMLParser::first_attribute(xml_node<> *node, string attr_name, T &attr_val) {
    auto attr = node->first_attribute(attr_name.c_str());
    if(attr){
        string attr_str = attr->value();
        attr_val = stringTo<T>(attr_str);
    }
}

void CometPepXMLParser::export_psm_info(vector<shared_ptr<PSMInfo>> &psm,xml_document<> &doc) {
    m_currentNode = breadth_first_search("msms_run_summary", &doc, false);
    vector<xml_node<> *> sourcefileNodes = find_all_siblings("msms_run_summary", m_currentNode);

    for (auto &x: sourcefileNodes) {
        string newfile = x->first_attribute("base_name")->value();
        string extension =  x->first_attribute("raw_data")->value();
        // sometimes, the raw_data is mzXML, NOT .mzXML, the dot is missing. So this is an ad hoc fix.
        if(extension.at(0)!='.') extension = '.' + extension;

        // if newfile does not end with extension, append it.
        if(newfile.length() < extension.length() || newfile.substr(newfile.length() - extension.length()) != extension)
        {
            newfile += extension;
        }

        m_allSourceFiles.push_back(newfile);
        spdlog::get("A")->info("Source filename: {}", newfile);

        string keywords = "spectrum_query";
        m_currentNode = breadth_first_search(keywords, &doc, false);
        vector<xml_node<> *> psm_nodes = find_all_siblings("spectrum_query", m_currentNode);

        // Progress ps(psm_nodes.size(), "Parsing spectrum in CometPepXML file");
        for (auto spectrum_query_node: psm_nodes) {
            // ps.increase();
            psm.push_back(make_shared<PSMInfo>(spectrum_query_node));
        }

        m_filename2indeRange[newfile]=make_pair(psm.size()-psm_nodes.size(), psm.size());
    }
    for(int i = 0; i < psm.size(); i ++)
    {
        int scan = (*psm[i]).start_scan;
        if(m_scan2idxvec.count(scan)==0)    {
            m_scan2idxvec[scan]={i};
        }else{
            m_scan2idxvec[scan].push_back(i);
        }
    }
    m_use_scan2idxvec = true;
    // cout << "Mapping of scan to idx created" << endl;
}

CometPepXMLParser::CometPepXMLParser(string filename) {
    m_use_scan2idxvec = false;
    m_currentNode = nullptr;
    spdlog::get("A")->info("Loading file: {}", filename);
    // cout << "[Info] Loading file " << filename << endl;
    long len = 0;
    File::getfilesize(filename, len);

    // cout << "File size :" << len * 1.0 / 1024 / 1024 << " MB" << endl;

    FILE *pfile = fopen(filename.c_str(), "r");
    char *  m_buf = new char[len + 1];
    if (nullptr == m_buf ) {
        cout << "Fail to get memory of size " << len << endl;
        throw invalid_argument("can not get enought memory");
    }
    long long nsize = fread(m_buf, sizeof(char), len, pfile);
    m_buf[len] = '\0';
    fclose(pfile);

    {
        xml_document<> doc;
        doc.parse<0>(m_buf);
        // cout << "[Info] extract PSMs " << endl;
        export_psm_info(psm,doc);
        // cout << "[Info] " << psm.size() << " PSMs extracted from input file" << endl;
        doc.clear();

        for(int i = 0 ; i < psm.size(); i ++)   {
            // what if one scan occurs twice?
//        if(m_scan2psminfoidx.count(psm[i].start_scan)==1){
//            int pre_idx = m_scan2psminfoidx.find(psm[i].start_scan)->second;
//            cout << "duplicate scan found (old): " << pre_idx << "\t" << psm[pre_idx].start_scan << "\t" << psm[pre_idx].charge << endl;
//            cout << "duplicate scan found (new): " << i << "\t" << psm[i].start_scan << "\t" << psm[i].charge << endl;
//        }
            m_scan2psminfoidx.insert(make_pair(psm[i]->start_scan, i));
//        if(m_scan2psminfoidx.count(psm[i].start_scan)==2){
//            cout << "which is used?: idx used: " << m_scan2psminfoidx.find(psm[i].start_scan)->second   << endl;
//        }
            m_spectrumName2psminfoidx.insert(make_pair(psm[i]->spectrum, i));
        }

    }

    delete[] m_buf;

    long memSize = accumulate(psm.begin(), psm.end(), 0L,[](long a, shared_ptr<PSMInfo> b){return a + b->getMemSize();});
    // cout << "[PSM]  Total size of memory taken " << memSize/1024/1024 << " MB" << endl;
}

CometPepXMLParser::~CometPepXMLParser() {
    //cout << "[Info] calling " << __FUNCTION__  << endl;
};

bool CometPepXMLParser::getPsmInfo(PSMInfo &psminfo) {
    int charge = psminfo.charge;
    int scan = psminfo.start_scan;
    return getPSMInfobyScanCharge(scan, charge, psminfo);

}
bool CometPepXMLParser::getPSMInfobyScanCharge(int scan, int charge, PSMInfo &psminfo) {
    bool ret = false;
    if(m_use_scan2idxvec){
        if(m_scan2idxvec.count(scan)==0){
            cout << "Error: scan=" << scan << " out of range" << endl;
        } else{
            // looking for proper charge state.
            vector<int> idxlist = m_scan2idxvec[scan];
            for(auto i : idxlist){
                if(psm[i]->charge == charge){
                    psminfo=*psm[i];
                    ret = true;
                    break;
                }
            }
        }

    }else{
        cout << "Error: please create map of scan to index first" << endl;
    }
    return ret;
}


// Attention: this function assumes that scan numbers are sorted ascendingly
//
bool CometPepXMLParser::getPSMInfobyScan(int scan, PSMInfo &psminfo) {
    // WARNING: DO NOT USE THIS FUNCTION
    //cout << "[Warning] DO NOT USE THIS FUNCTION: DEPRECATED" << endl;
    bool found = false;
    int low = 0, high = psm.size();
    while(low < high){
        // check mid
        int mid = low + (high-low)/2;
        PSMInfo &currentPSM = *psm[mid];
        if(currentPSM.end_scan>scan){
            high = mid;
        } else if(currentPSM.end_scan <scan)
        {
            low = mid+1;
        } else{
            found = true;
            psminfo = currentPSM;
            break;
        }

    }
    return found;
}

void CometPepXMLParser::exportToTXT(string filename) {
    ofstream fout(filename);
    fout << "#spectrum\tscan\tprecursor\tcharge\trt\thit_rank\tpeptide\txcorr\texpect\tspscore\thyperscore\tnextscore" << endl;
    for(int i = 0; i < psm.size(); i ++)    {
        PSMInfo &currentPSM = *psm[i];
        fout << currentPSM.spectrum << "\t" << currentPSM.start_scan << "\t" << currentPSM.precursorNeutMass
             << "\t" << currentPSM.charge << "\t" << currentPSM.retention_time ;
        for(int j = 0; j < currentPSM.searchhits.size(); j ++) {
            SearchHit &sh = *currentPSM.searchhits[j];
            fout << "\t" <<"hit_" << j+1 << "\t" << sh.m_peptide << "\t" <<  sh.m_xcorr << "\t" << sh.m_expect
                 << "\t" << sh.m_spscore << "\t" << sh.m_hyperscore<< "\t" << sh.m_nextscore;
        }
        fout << endl;
    }
    fout.close();
}

bool CometPepXMLParser::updateGtInfo(SPsmAnnotation &gtinfo) {
    bool ret = false;
    PSMInfo psminfo;

    if(getPSMInfobyScan(gtinfo.ms2_scan,psminfo))    {
        if(psminfo.charge!=0){
            gtinfo.charge = psminfo.charge;
            gtinfo.precursorMz = psminfo.precursorNeutMass / psminfo.charge + proton_mass;
        }

        if (psminfo.searchhits.size()>0)        {
            SearchHit &sh0 = *psminfo.searchhits[0];
            gtinfo.peptideseq = sh0.m_peptide;
            gtinfo.score = sh0.m_expect;
            gtinfo.nterm_mass = sh0.m_nterm_mass;
            gtinfo.cterm_mass = sh0.m_cterm_mass;
            gtinfo.modificationstr = CPosMass(sh0.m_pos_mass).toString();
            ret = true;
        }
    }else{
        gtinfo.peptideseq = "UNKNOWN";
        gtinfo.score = -1;
        gtinfo.nterm_mass = -1;
        gtinfo.cterm_mass = -1;


    }
    return ret;
}

void CometPepXMLParser::getscoreandpep_mod(int scan, string &peptide, double &searchscore, modification &modinfo,
                                           double &parentMZ, int &charge) {
    PSMInfo psminfo;
    if(getPSMInfobyScan(scan,psminfo))    {
        charge = psminfo.charge;
        parentMZ = psminfo.precursorNeutMass / psminfo.charge + proton_mass;
        if(psminfo.searchhits.size()>0)        {
            SearchHit &sh0 = *psminfo.searchhits[0];
            peptide = sh0.m_peptide;
            searchscore = sh0.m_expect;
            modinfo.nterm_mass = sh0.m_nterm_mass;
            modinfo.cterm_mass = sh0.m_cterm_mass;
            modinfo.pos_mass = sh0.m_pos_mass;
        }
    }else{
        peptide = "UNKNOWN";
        searchscore = -1;
        modinfo.nterm_mass = -1;
        modinfo.cterm_mass = -1;
        parentMZ = -1;
        charge = -1;
    }
}

int CometPepXMLParser::getPsmNum() {return psm.size();}

void CometPepXMLParser::getPSMInfobyindex(int i, PSMInfo &psminfo) {
    psminfo= *psm[i];
}

pair<int, int> CometPepXMLParser::getIndexRange(string filename) {
    return m_filename2indeRange[filename];
}

// get source file name from the spectrum name.
vector<string> CometPepXMLParser::getAllSoruceFiles() {
        std::set<string> allSourceFiles;
    // add all source files to set 
    // for each item in m_allSourceFiles, find its file name with CFile
    for(auto &x: m_allSourceFiles) {
        
        string basename = File::CFile(x).basename;
        // if(extension.length()>0) {
        //     filename = x.substr(0, x.length() - extension.length() - 1);
        // }
        allSourceFiles.insert(basename);
    }
    // allSourceFiles.insert(m_allSourceFiles.begin(), m_allSourceFiles.end());

    
    // find all source filenames from the spectrum name.
    for(auto &x: psm) {
        vector<string> items;
        split_string(x->m_basename, items, '.');
        string name = "";
        if(items.size() == 4){
            name = items[0];
            // m_allSoruceFiles.push_back(items[0]);
        }else{
            cout << "Found '.' in filename: " << x->m_basename << endl;
            for(auto &y: items) {
                if(name.length()>0) name += ".";
                name += y;
            }
        }
        allSourceFiles.insert(name);
        
        m_allSourceFiles = vector<string>(allSourceFiles.begin(), allSourceFiles.end());
        
    }
    // print all the source file name 
        // for(auto &x: m_allSourceFiles) {
        //     cout << "Source file: " << x << endl;
        // }
    return m_allSourceFiles;
}


void PSMInfo::print() {
    cout << "start_scan = " << start_scan << endl
         << "end_scan = " << end_scan << endl
         << "precursor_neutral_mass = " << precursorNeutMass << endl
         << "index = " << index << endl
         << "retention_time = " << retention_time << endl
         << "charge = " << charge << endl
         << "spectrum = " << spectrum << endl;
    for (auto each : searchhits) each->print();

}
/*
 *   <spectrum_query spectrum="decoy_human_hcd_selected_new.00001.00001.2"
 *   spectrumNativeID="AAAAAAAAAAAAAAAGAGAGAK/2"
 *   start_scan="1"
 *   end_scan="1"
 *   precursor_neutral_mass="1595.837979"
 *   assumed_charge="2"
 *   index="1">
 */

// todo: to be updated!
PSMInfo::PSMInfo(xml_node<> *spectrum_query_node):PSMInfo(){
    if(spectrum_query_node== nullptr)
    {
        cout << "empty spectrum " << spectrum_query_node  << endl;
        throw runtime_error("empty spectrum query node!");
    }    
    first_attribute(spectrum_query_node, "spectrum", spectrum);

    spectrum = spectrum_query_node->first_attribute("spectrum")->value();
    m_basename = spectrum; // make the basename the same as spectrum.
    first_attribute(spectrum_query_node,"start_scan", start_scan);
    start_scan = stoi(spectrum_query_node->first_attribute("start_scan")->value(),nullptr);
    first_attribute(spectrum_query_node,"end_scan", end_scan);
    end_scan = stoi(spectrum_query_node->first_attribute("end_scan")->value(), nullptr);
    first_attribute(spectrum_query_node,"assumed_charge", charge);
    charge = stoi(spectrum_query_node->first_attribute("assumed_charge")->value(),nullptr);
    first_attribute(spectrum_query_node,"index", index);
    index = stoi(spectrum_query_node->first_attribute("index")->value(), nullptr);
    first_attribute(spectrum_query_node,"precursor_neutral_mass", precursorNeutMass);
    precursorNeutMass = strtod(spectrum_query_node->first_attribute("precursor_neutral_mass")->value(), nullptr);
    // Remember always test before you use it
    first_attribute(spectrum_query_node,"retention_time_sec", retention_time);
    retention_time=0;
    auto attr=spectrum_query_node->first_attribute("retention_time_sec");
    if(attr){
        retention_time = strtod(attr->value(),nullptr);
    }

    xml_node<> *first_hit = breadth_first_search("search_hit", spectrum_query_node, false);
    if(first_hit!= nullptr)   {
        vector<xml_node<> *> siblings = find_all_siblings("search_hit", first_hit);
        for (auto & sibling : siblings) {
            searchhits.push_back(make_shared<SearchHit>(sibling));
        }
    }
}

PSMInfo::PSMInfo(const PSMInfo &other) {
    start_scan = other.start_scan;
    end_scan = other.end_scan;
    precursorNeutMass = other.precursorNeutMass;
    index = other.index;
    retention_time = other.retention_time;
    charge = other.charge;
    spectrum = other.spectrum;
    searchhits = other.searchhits;
    m_basename = other.m_basename;
}

PSMInfo::PSMInfo() {
    start_scan = -1;
    end_scan = -1;
    precursorNeutMass = -1;
    index = -1;
    retention_time = -1;
    charge = -1;
    spectrum = "";
    m_basename="";
}


// changed!
bool PSMInfo::isDecoy(bool useAlternativeProt, int hitrank) {
    if(searchhits.empty()) {
        cout << "Error: empty hits!" << endl;
        return false;
    }
    return searchhits.at(hitrank)->isDecoy(useAlternativeProt);

}

void PSMInfo::printSearchHit(int k) {
    searchhits.at(k)->print();
}

string PSMInfo::getProtein_UseAlterProteinIfItsNotDecoy(bool useAlternativeProt, int hitrank) {
    return searchhits.at(hitrank)->getProtein_UseAlterProteinIfItsNotDecoy(useAlternativeProt);
}

string PSMInfo::getFirstAlterProteinNotDecoy() {
    string prot = "";
    if (not searchhits.empty()) {
        prot=searchhits.at(0)->getFirstAlterProteinNotDecoy();
    }
    return prot;
}

double PSMInfo::getParentMz() {
    double protonmass = 1.007276;
    if(charge==0) {
        throw logic_error("zero charge found, divided by zero on precursor mz calculation!");
    }
    return precursorNeutMass/charge + protonmass;

}

void SearchHit::print() {
    cout << "m_hit_rank = " << m_hit_rank << endl
         << "m_peptide = " << m_peptide << endl
         << "m_preAA = " << m_preAA << endl
         << "m_nextAA = " << m_nextAA << endl
         << "m_protein = " << m_protein << endl;
    for(auto &x: m_alternative_proteins) cout << "alternative Protein: " << x << endl;
    cout << "m_num_matched_ion = " << m_num_matched_ion << endl
         << "m_tot_num_ions = " << m_tot_num_ions << endl
         << "m_num_missed_cleavages = " << m_num_missed_cleavages << endl
         << "m_modified_peptide = " << m_modified_peptide << endl
         << "m_nterm_mass = " << m_nterm_mass << endl
         << "m_cterm_mass = " << m_cterm_mass << endl
         << "m_xcorr = " << m_xcorr << endl
         << "m_deltacn = " << m_deltacn << endl
         << "m_deltacnstar = " << m_deltacnstar << endl
         << "m_spscore = " << m_spscore << endl
         << "m_sprank = " << m_sprank << endl
         << "m_expect = " << m_expect << endl
         << "m_hyperscore = " << m_hyperscore << endl
         << "m_nextscore = " << m_nextscore << endl
         << "m_massdiff = " << m_massdiff << endl
         << "m_peptideprophet_prob = " << m_peptideprophet_prob << endl
         << "m_iprophet_prob= " << m_iprophet_prob << endl << endl;
}

/*
 * <search_hit hit_rank="1" peptide="GGLESASGRRGAPEPR" peptide_prev_aa="R"
 * peptide_next_aa="L" protein="DECOY_sp|Q15773|MLF2_HUMAN" num_tot_proteins="3" num_matched_ions="6"
 tot_num_ions="60" calc_neutral_pep_mass="1595.812815" massdiff="0.025168"
 num_tol_term="2" num_missed_cleavages="2" num_matched_peptides="6087">
    <alternative_protein protein="DECOY_tr|Q5U0N1|Q5U0N1_HUMAN"/>
    <alternative_protein protein="DECOY_tr|A8K1F4|A8K1F4_HUMAN"/>
    <search_score name="xcorr" value="1.166"/>
    <search_score name="deltacn" value="0.021"/>
    <search_score name="deltacnstar" value="0.000"/>
    <search_score name="spscore" value="45.2"/>
    <search_score name="sprank" value="5"/>
    <search_score name="expect" value="6.59E+00"/>
   </search_hit>
 *
 * */


SearchHit::SearchHit(xml_node<> *node)  {
    // node attr list:
    // hit_rank
    m_hit_rank = atoi(node->first_attribute("hit_rank")->value());

/*
<search_hit hit_rank="1" peptide="HNLGHGHK" protein="P01042" num_tot_proteins="1" num_matched_ions="0" calc_neutral_pep_mass="1356.77805" massdiff="0.0005770" 
protein_descr="Kininogen-1 [OS=Homo sapiens]" protein_mw="71.912" calc_pI="6.81">
          <modification_info modified_peptide="[K].hNLGHGHk.[H]">
            <mod_aminoacid_mass position="8" mass="357.25789496" />
          </modification_info>
          <modification_info modified_peptide="[K].hNLGHGHk.[H]" mod_nterm_mass="229.162932" />
          <search_score name="XCorr" value="4.49736452102661" />
          <search_score name="Percolator q-Value" value="0.0003473" />
          <search_score name="Percolator PEP" value="1.58E-05" />
          <search_score name="Percolator SVMScore" value="1.796" />
          <analysis_result analysis="percolator">
            <percolator_result probability="0.9999842" PEP="1.58E-05" q-Value="0.0003473" />
          </analysis_result>
        </search_hit>
*/


    m_num_matched_ion = -1; // unknown
    if (node->first_attribute("num_matched_ions")) {
        m_num_matched_ion = atoi(node->first_attribute("num_matched_ions")->value());
    }
    
    m_num_missed_cleavages = -1;
    if(node->first_attribute("num_missed_cleavages")){
        m_num_missed_cleavages = atoi(node->first_attribute("num_missed_cleavages")->value());
    }

    m_tot_num_ions = -1;
    if (node->first_attribute("tot_num_ions")) {
        m_tot_num_ions = atoi(node->first_attribute("tot_num_ions")->value());
    }
    m_peptide="";
    if (node->first_attribute("peptide")) {
        m_peptide = node->first_attribute("peptide")->value();
        // m_peptide = node->first_attribute("peptide")->value();
    }
    
    
    m_preAA = "";
    if (node->first_attribute("peptide_prev_aa")) {
        m_preAA = node->first_attribute("peptide_prev_aa")->value();
    }
    m_nextAA = "";
    if (node->first_attribute("peptide_next_aa")) {
        m_nextAA = node->first_attribute("peptide_next_aa")->value();
    }
    // m_preAA = node->first_attribute("peptide_prev_aa")->value();
    // m_nextAA = node->first_attribute("peptide_next_aa")->value();
    m_protein = "";
    if (node->first_attribute("protein")) {
        m_protein = node->first_attribute("protein")->value();
    }
    // m_protein = node->first_attribute("protein")->value();

    m_massdiff = 0;
    if (node->first_attribute("massdiff")) {
        m_massdiff = atof(node->first_attribute("massdiff")->value());
    }
    // m_massdiff = atof(node->first_attribute("massdiff")->value());
    m_modified_peptide = m_peptide;
    m_peptideprophet_prob = UNKOWN_PEPTIDEPROPHET_SCORE;
    m_iprophet_prob = UNKOWN_PEPTIDEPROPHET_SCORE;

    m_nterm_mass = 0;
    m_cterm_mass = 0;

//    cout << "where" << endl;
//    <ptm_result localization="14" best_score_with_ptm="17.478317" score_without_ptm="15.668178" second_best_score_with_ptm="15.668178" ptm_mass="41.030273"/>
    xml_node<> *ptm_results = node->first_node("ptm_result");
    if(ptm_results){
//        cout << "gt ptm" << endl;
        // how many siblings?
        vector<xml_node<>* > ptmResAll = find_all_siblings("ptm_result", ptm_results);
        for(auto &eachPtmRes : ptmResAll)  {
            string localization = eachPtmRes->first_attribute("localization")->value();
            m_localizations.push_back(localization);
        }

    }


    xml_node<> *alternative_protein = node->first_node("alternative_protein");
    if(alternative_protein)
    {
        /*
         * <search_hit hit_rank="1" peptide="EIEEKLR" peptide_prev_aa="K" peptide_next_aa="K"
         *      protein="DECOY_sp|Q5VT25|MRCKA_HUMAN" num_tot_proteins="6" num_matched_ions="10"
         *      tot_num_ions="12" calc_neutral_pep_mass="915.502546" massdiff="-0.025705" num_tol_term="2"
         *      num_missed_cleavages="1" num_matched_peptides="19147" protein_descr="Decoy sequence">
         *
         * <alternative_protein protein="sp|Q5PSV4|BRM1L_HUMAN" protein_descr="Breast cancer ... PE=1 SV=2"
         *      num_tol_term="1" peptide_prev_aa="S" peptide_next_aa="R"/>
         *
         * <alternative_protein protein="sp|Q08379|GOGA2_HUMAN" protein_descr="Golgin subfamily ... GN=GOLGA2 PE=1 SV=3"
         *      num_tol_term="1" peptide_prev_aa="S" peptide_next_aa="V"/>
         * ...
         *
         * */

        vector<xml_node<>* > alterProteins = find_all_siblings("alternative_protein", alternative_protein);
        for(auto &eachprot : alterProteins)  {
            string alterProt = eachprot->first_attribute("protein")->value();
            m_alternative_proteins.push_back(alterProt);
        }
    }

    xml_node<> *modification = node->first_node("modification_info");
    if (modification) {
        /*
         * This is the modification in X!Tandem
         *  <modification_info mod_nterm_mass="29.0391">
               <mod_aminoacid_mass position="4" mass="156.1263" />
               <mod_aminoacid_mass position="7" mass="156.1263" />
            </modification_info>

            This is the modification from Comet:
            <modification_info modified_peptide="QNRM[147]MQNYRGLLSPSDSR">
                <mod_aminoacid_mass position="4" mass="147.035385"/>
            </modification_info>

            MSFragger:
            <modification_info modified_peptide="M[147]M[147]EEEKDCK">
            <mod_aminoacid_mass position="1" mass="147.035385"/>
            <mod_aminoacid_mass position="2" mass="147.035385"/>
            <mod_aminoacid_mass position="8" mass="160.030649"/>
            </modification_info>

         * */
        // attention, the modified_peptide here is not correct! we need nterm modification, and C[160]
        // <modification_info modified_peptide="M[147]M[147]EEEKDCK">
        // <mod_aminoacid_mass position="1" mass="147.035385"/>
        // <mod_aminoacid_mass position="2" mass="147.035385"/>
        // <mod_aminoacid_mass position="8" mass="160.030649"/>
        // </modification_info>
        // 


        xml_attribute<char> *attr = modification->first_attribute("modified_peptide");
        bool modified_pep_updated = false;
        if (attr) {
            modified_pep_updated = true;
            m_modified_peptide = attr->value();
        }
        attr = modification->first_attribute("mod_nterm_mass");
        if (attr) m_nterm_mass = atof(attr->value());
        attr = modification->first_attribute("mod_cterm_mass");
        if (attr) m_cterm_mass = atof(attr->value());

        


        // get all the mod_aminoacid_mass info
        xml_node<> *mod_aminoacid_mass_one = modification->first_node("mod_aminoacid_mass");// then find all sibings;
        vector<xml_node<> *> all_mod_aminoacid_mass = find_all_siblings("mod_aminoacid_mass", mod_aminoacid_mass_one);

        for (auto eachmod: all_mod_aminoacid_mass) {
            int pos = atoi(eachmod->first_attribute("position")->value());
            double mass = strtod(eachmod->first_attribute("mass")->value(),nullptr);
            m_pos_mass[pos] = mass;
        }

        // refresh peptide.
        if(not modified_pep_updated){
            // Note: without modified_pep_updated variable, the program crash!
            // what is the reason? C[160] modified or not?
            string tmp_mod_string = "";
            if(m_nterm_mass>1e-6){
                // non zero nterm mass
                tmp_mod_string+='n';
                int imass=(int)round(m_nterm_mass);
                tmp_mod_string += "["+to_string(imass)+"]";

            }
            for(int k=0; k < m_peptide.length(); k ++){
                tmp_mod_string += m_peptide.at(k);
                if(m_pos_mass.find(k+1)!=m_pos_mass.end()){
                    double mass = m_pos_mass.at(k+1);
                    // we need amino acid mass table
                    double amino_acid_mass = 0; // no need.
                    int imass=(int)round(mass+amino_acid_mass);
                    tmp_mod_string += "["+to_string(imass)+"]";
                }
            }
            m_modified_peptide = tmp_mod_string;

        }


    }


/*
 * <search_hit hit_rank="1" peptide="GGLESASGRRGAPEPR" peptide_prev_aa="R"
 * peptide_next_aa="L" protein="DECOY_sp|Q15773|MLF2_HUMAN" num_tot_proteins="3" num_matched_ions="6"
 tot_num_ions="60" calc_neutral_pep_mass="1595.812815" massdiff="0.025168"
 num_tol_term="2" num_missed_cleavages="2" num_matched_peptides="6087">
    <alternative_protein protein="DECOY_tr|Q5U0N1|Q5U0N1_HUMAN"/>
    <alternative_protein protein="DECOY_tr|A8K1F4|A8K1F4_HUMAN"/>
    <search_score name="xcorr" value="1.166"/>
    <search_score name="deltacn" value="0.021"/>
    <search_score name="deltacnstar" value="0.000"/>
    <search_score name="spscore" value="45.2"/>
    <search_score name="sprank" value="5"/>
    <search_score name="expect" value="6.59E+00"/>
   </search_hit>
 *
 * */

    xml_node<> *score = node->first_node("search_score");
    vector<xml_node<> *> siblings = find_all_siblings("search_score", score);

    m_xcorr = 0;
    m_deltacn = 0;
    m_deltacnstar = 0;
    m_sprank = 0;
    m_spscore = 0;
    m_expect = 0;
    m_hyperscore = 0;
    m_nextscore = 0;
    for (auto each_score: siblings) {
        if (strcmp("xcorr", each_score->first_attribute("name")->value()) == 0) {
            m_xcorr = strtod(each_score->first_attribute("value")->value(),nullptr);
        }
        if (strcmp("deltacn", each_score->first_attribute("name")->value()) == 0) {
            m_deltacn = strtod(each_score->first_attribute("value")->value(), nullptr);
        }
        if (strcmp("deltacnstar", each_score->first_attribute("name")->value()) == 0) {
            m_deltacnstar = strtod(each_score->first_attribute("value")->value(),nullptr);
        }
        if (strcmp("spscore", each_score->first_attribute("name")->value()) == 0) {
            m_spscore = strtod(each_score->first_attribute("value")->value(), nullptr);
        }
        if (strcmp("sprank", each_score->first_attribute("name")->value()) == 0) {
            m_sprank = strtod(each_score->first_attribute("value")->value(), nullptr);
        }
        if (strcmp("expect", each_score->first_attribute("name")->value()) == 0) {
            m_expect = strtod(each_score->first_attribute("value")->value(), nullptr);
        }
        if (strcmp("hyperscore", each_score->first_attribute("name")->value()) == 0) {
            m_hyperscore = strtod(each_score->first_attribute("value")->value(),nullptr);
        }
        if (strcmp("nextscore", each_score->first_attribute("name")->value()) == 0) {
            m_nextscore = strtod(each_score->first_attribute("value")->value(), nullptr);
        }
    }
    xml_node<> *peptideprophet = nullptr, *iprophet = nullptr, *analysis_results = nullptr;
    analysis_results = node->first_node("analysis_result");
    vector<xml_node<> *> ar_siblings = find_all_siblings("analysis_result", analysis_results);
    for (auto each_ar : ar_siblings) {
        if (strcmp(each_ar->first_attribute("analysis")->value(), "peptideprophet") == 0) {
            peptideprophet = each_ar->first_node("peptideprophet_result");
            if (peptideprophet != nullptr) {
                string tmpvalue = peptideprophet->first_attribute("probability")->value();
                m_peptideprophet_prob = atof(tmpvalue.c_str());
            }
        }

        if (strcmp(each_ar->first_attribute("analysis")->value(), "interprophet") == 0) {
            iprophet = each_ar->first_node("interprophet_result");
            if (iprophet != nullptr) {
                string tmpvalue = iprophet->first_attribute("probability")->value();
                m_iprophet_prob = atof(tmpvalue.c_str());
            }
        }

    }

}

SearchHit::SearchHit(int hit_rank, string peptide, string preAA, string nextAA, string protein, int num_matched_ion,
                     int tot_num_ions, int num_missed_cleavages, string modified_peptide, double xcorr,
                     double deltacn, double deltacnstar, double spscore, double sprank, double expect,
                     double massdiff, double peptideprophet_prob, double iprophet_prob, double nterm_mass,
                     double cterm_mass, map<int, double> pos_mass, vector<string> &altProts) {
    cout << "I guess it is never used so this will never be printed " << endl;
    m_hit_rank = hit_rank;
    m_peptide = peptide;
    m_preAA = preAA;
    m_nextAA = nextAA;
    m_protein = protein;
    m_alternative_proteins = altProts;
    m_num_matched_ion = num_matched_ion;
    m_tot_num_ions = tot_num_ions;
    m_num_missed_cleavages = num_missed_cleavages;
    m_modified_peptide = modified_peptide;
    m_xcorr = xcorr;
    m_deltacn = deltacn;
    m_deltacnstar = deltacnstar;
    m_spscore = spscore;
    m_sprank = sprank;
    m_expect = expect;
    m_massdiff = massdiff;
    m_peptideprophet_prob = peptideprophet_prob;
    m_iprophet_prob = iprophet_prob;

    m_cterm_mass = cterm_mass;
    m_nterm_mass = nterm_mass;
    m_pos_mass = pos_mass;
}

SearchHit::SearchHit(const SearchHit &other) {
    m_hit_rank = other.m_hit_rank;
    m_peptide = other.m_peptide;
    m_preAA = other.m_preAA;
    m_nextAA = other.m_nextAA;
    m_protein = other.m_protein;
    m_alternative_proteins = other.m_alternative_proteins;
    m_num_matched_ion = other.m_num_matched_ion;
    m_tot_num_ions = other.m_tot_num_ions;
    m_num_missed_cleavages = other.m_num_missed_cleavages;
    m_modified_peptide = other.m_modified_peptide;
    m_xcorr = other.m_xcorr;
    m_deltacn = other.m_deltacn;
    m_deltacnstar = other.m_deltacnstar;
    m_spscore = other.m_spscore;
    m_sprank = other.m_sprank;
    m_expect = other.m_expect;
    m_hyperscore = other.m_hyperscore;
    m_nextscore = other.m_nextscore;

    m_massdiff = other.m_massdiff;
    m_peptideprophet_prob = other.m_peptideprophet_prob;
    m_iprophet_prob = other.m_iprophet_prob;
    m_cterm_mass = other.m_cterm_mass;
    m_nterm_mass = other.m_nterm_mass;
    m_pos_mass = other.m_pos_mass;

}

string SearchHit::getFirstAlterProteinNotDecoy() {
    string prot = "";

    for (auto &x : m_alternative_proteins) {
        if (x.find("DECOY") != string::npos)  {  continue; }
        else {  prot = x; break; }
    }

    return prot;
}

bool SearchHit::isDecoy(bool useAlternativeProt) {
    if(useAlternativeProt){
        return m_protein.find("DECOY")!=string::npos and getFirstAlterProteinNotDecoy()=="";

    } else{
        return m_protein.find("DECOY")!=string::npos;
    }

}

string SearchHit::getProtein_UseAlterProteinIfItsNotDecoy(bool useAlternativeProt) {
    // prot--> not decoy , return ! YES prot
    // prot & alt_prot --> isDecoy, Can not improve! return prot

    if(m_protein.find("DECOY")==string::npos){
        return m_protein;
    }else if ( isDecoy(useAlternativeProt)){
        return m_protein;
    } else {
        // alt_prot --> NOT decoy, return alt_prot
        return getFirstAlterProteinNotDecoy();
    }


}

PeptideProphetParser::PeptideProphetParser(string &filename) {
    m_filename = filename;
    initialize();
}

PeptideProphetParser::~PeptideProphetParser() {
    // cout << "[Info] calling " << __FUNCTION__  << endl;
    delete[] m_buf;
}



void PeptideProphetParser::export_psm_info(vector<PSMInfo> &psm) {
    m_currentNode = breadth_first_search("msms_run_summary", &doc, false);
    vector<xml_node<> *> sourcefileNodes = find_all_siblings("msms_run_summary", m_currentNode);


    for (auto &x: sourcefileNodes) {
        string newfile = x->first_attribute("base_name")->value();
        newfile += x->first_attribute("raw_data")->value();
        m_allSourceFiles.push_back(newfile);
        spdlog::get("A")->info("Get source fileanme: {}", newfile);

        m_currentNode = breadth_first_search("spectrum_query", x, false);
        if (m_currentNode == nullptr) {
            cout << "NO spectrum_query node found!!" << endl;
            continue;
        }
        vector<xml_node<> *> spectrum_queries = find_all_siblings("spectrum_query", m_currentNode);

        // Progress ps(spectrum_queries.size(), "Parsing PSMs");
        for (auto &y: spectrum_queries) {
            // ps.increase();
            psm.emplace_back(y);
        }
        m_filename2indeRange[newfile]=make_pair(psm.size()-spectrum_queries.size(), psm.size());
    }
    cout << "[Info] " << psm.size() << " PSMs extracted from input file" << endl;
}

void PeptideProphetParser::filterWithThreshold(double min_probability_score, vector<PSMInfo> &newpsm) {
    Progress ps(m_psm.size());
    for (int i = 0; i < m_psm.size(); i++) {
        ps.increase();
        bool ret = false;
        if(not m_psm[i].searchhits.empty())       {
            if(m_use_iProb){
                if(m_psm[i].searchhits[0]->m_iprophet_prob >= min_probability_score) ret = true;
            }
            else
            {
                if(m_psm[i].searchhits[0]->m_peptideprophet_prob >= min_probability_score) ret = true;
            }
        }
        if(ret)  newpsm.push_back(m_psm[i]);
    }


}

double PeptideProphetParser::getThresholdForFDR(double fdr_threshold, bool use_iProb) {
    cout << "[Info] Exporting PeptideProphet results..." << endl;
    string summary_tag = "peptideprophet_summary";
    if(use_iProb) summary_tag = "interprophet_summary";
    xml_node<> *any_prophet_summary = breadth_first_search(summary_tag, &doc, false);
    // input file
    cout << "[Info] Input file section: " << endl;
    m_currentNode = breadth_first_search("roc_error_data", any_prophet_summary, false);
    vector<xml_node<> *> roc_error_data_notes = find_all_siblings("roc_error_data", m_currentNode);

    cout << "[Info] " << roc_error_data_notes.size() << " roc error sections found!" << endl;
    for (auto x: roc_error_data_notes) {
        if (strcmp(x->first_attribute("charge")->value(), "all") == 0) {
            m_currentNode = x;
            cout << "[Info] found roc error section with charge=all" << endl;
            break;
        } else {
            cout << x->first_attribute("charge")->value() << endl;
        }
    }
    // get all the information inside currentNode
    m_currentNode = m_currentNode->first_node("error_point");
    vector<xml_node<> *> error_points = find_all_siblings("error_point", m_currentNode);

    double min_probability_score = 0;
    int num_corr = 0, num_incorr = 0;
    for (auto x: error_points) {
        if (fabs(std::strtod(x->first_attribute("error")->value(), nullptr) - fdr_threshold) < 1e-6) {
            min_probability_score = std::strtod(x->first_attribute("min_prob")->value(), nullptr);
            num_corr = atoi(x->first_attribute("num_corr")->value());
            num_incorr = atoi(x->first_attribute("num_incorr")->value());
            cout << "[Info] find FDR<=" << setprecision(4) << fdr_threshold << " when probability >= "
                 << min_probability_score << "\tnum_corr=" << num_corr << "\t num_incorr=" << num_incorr << endl;
            break;
        }
    }
    return min_probability_score;
}

void PeptideProphetParser::filter_with_FDR(double fdr_threshold, vector<PSMInfo> &newpsm) {
    double min_probability_score = this->getThresholdForFDR(fdr_threshold, m_use_iProb);
    this->filterWithThreshold(min_probability_score, newpsm);
    cout << "[Info] Filter with FDR < "<< fdr_threshold <<", PSM number is " << newpsm.size() << endl;
}

bool PeptideProphetParser::isPSMSignificant(PSMInfo &psminfo) {
    if(psminfo.searchhits.empty()) {
        return false;
    }
    else if(m_use_iProb){
        return psminfo.searchhits[0]->m_iprophet_prob >= m_threshold;
    }
    else{
        return psminfo.searchhits[0]->m_peptideprophet_prob >= m_threshold;
    }
}

void PeptideProphetParser::initialize() {
    cout << "[Info] Loading file " << m_filename << endl;
    long len = 0;
    File::getfilesize(m_filename, len);
    FILE *pfile = fopen(m_filename.c_str(), "r");

    m_buf = new char[len+1];
    if (m_buf == nullptr)  {
        cout << "Fail to get pointer" << endl;
        exit(0);
    }
    long nsize = fread(m_buf, sizeof(char),len,pfile);
    cout << "[Info] Loaded bytes: " << nsize << " allowed size: " << len+1 << endl;
    m_buf[len] = '\0';
    fclose(pfile);
    doc.parse<0>(m_buf);    // 0 means default parse flags
    cout << "[Info] start to export: " << endl;
    export_psm_info(m_psm);
    cout << "[Info] Finish export" << endl;
    for(int i = 0 ; i < m_psm.size(); i ++)   {
        // what if one scan exist twice?
//        if(m_scan2psminfoidx.count(psm[i].start_scan)==1){
//            int pre_idx = m_scan2psminfoidx.find(psm[i].start_scan)->second;
//            cout << "duplicate scan found (old): " << pre_idx << "\t" << psm[pre_idx].start_scan << "\t" << psm[pre_idx].charge << endl;
//            cout << "duplicate scan found (new): " << i << "\t" << psm[i].start_scan << "\t" << psm[i].charge << endl;
//        }
        m_scan2psminfoidx.insert(make_pair(m_psm[i].start_scan, i));
//        if(m_scan2psminfoidx.count(psm[i].start_scan)==2){
//            cout << "which is used?: idx used: " << m_scan2psminfoidx.find(psm[i].start_scan)->second   << endl;
//        }
        m_spectrumName2psminfoidx.insert(make_pair(m_psm[i].spectrum, i));
    }

    m_use_iProb = false;
    if(m_filename.find("ipro.pep.xml")!=string::npos)    {
        m_use_iProb = true;
    }
    m_threshold = getThresholdForFDR(0.01, m_use_iProb);
    cout << "Find FDR < 1% threshold (default) for " << (m_use_iProb? "iProphet: ": "PeptidePropeht: ") << m_threshold << endl;
}

bool PeptideProphetParser::isPSMSignificant(int i) {
    PSMInfo &psminfo = m_psm[i];
    return isPSMSignificant(psminfo);
}

PeptideProphetParser::PeptideProphetParser(const PeptideProphetParser &other) {
    m_filename = other.getInputfile();
    initialize();
}

string PeptideProphetParser::getInputfile() const {
    return m_filename;
}

void PeptideProphetParser::getPSMInfobyindex(int i, PSMInfo &psminfo) {
    if(i<0 or i >= m_psm.size()) {
        cout << "index out of range" << endl;
    }
    psminfo = m_psm[i];
}

// todo: is it possible that one scan corresponds to multiple PSM...
// Warning: DO NOT USE THIS FUNCTION, ONE SPECTRA MIGHT CORRESPODING TO TWO PEPTIDES/PSMs.
bool PeptideProphetParser::getPSMInfobyScan(int scan, PSMInfo &psminfo) {
    const int cnt = m_scan2psminfoidx.count(scan);
    if (cnt == 0) {
        return false;
    } else if(cnt > 1){
        cout << "[Warning] One scan corresponds to " << cnt << " PSMs, first one will be used: scan " << scan << endl;
//        auto v = m_scan2psminfoidx.equal_range(scan);
//        for(auto x=v.first; x!=v.second; x++){
//            int k = x->second;
//            cout << "print -------------PSM " << endl;
//            m_psm[k].print();
//        }
    }
    psminfo = m_psm[m_scan2psminfoidx.find(scan)->second];
    return true;
}

// input: SpectrumName == SCAN, CHARGE, FILENAME
bool PeptideProphetParser::getPSMInfobySpectrumName(string spectrumName, PSMInfo &psminfo) {
    if (m_spectrumName2psminfoidx.count(spectrumName) == 0) return false;
    psminfo = m_psm[m_spectrumName2psminfoidx[spectrumName]];return true;
}

// input: SCAN, FILENAME
bool PeptideProphetParser::getPSMInfobyScanFileName(string filename, int scan, PSMInfo &psminfo) {
    string name = File::CFile(filename).basename;
    return getPSMInfo(scan, psminfo, name);
}

// input: SCAN, FILENAME
// go through charge states 1+ to 7+, once found a valid PSM, stop and return true
// otherwise, return false
bool PeptideProphetParser::getPSMInfo(int scan, PSMInfo &psminfo, string &basename) {
    bool found = false;
    for(int chg = 1; chg < 7; chg ++)    {
        found = getPSMInfo(scan, chg, psminfo, basename);
        if(found) break;
    }
    return found;
}


// The API: user should put charge, scan, and other information into the psminfo object.
// If found, return true, and the object psminfo will be refreshed;
// If not found, return false, and psminfo is untouched.
// input: index, charge, scan, filename, basename
// methd:
bool PeptideProphetParser::getPsmInfo(PSMInfo &psminfo) {
    int scan = psminfo.start_scan, chg = psminfo.charge;
    string name = psminfo.m_basename;
    bool found = false;
    if(chg == -1){
        found = getPSMInfo(scan, psminfo, name);
    }else{
        found = getPSMInfo(scan, chg, psminfo, name);
    }
    return found;
}

// PRIVATE
bool PeptideProphetParser::getPSMInfo(int scan, int chg, PSMInfo &psminfo, string basename) {
    string spectrumName = generateSpectrumName(scan, chg, basename);
    return getPSMInfobySpectrumName(spectrumName,psminfo);
}

bool PeptideProphetParser::getAllPSMsWithScanFileName(string filename, int scan, vector<PSMInfo> &psms) {
    bool found =false;
    int charge_start = 1, charge_end = 7;
    string name = File::CFile(filename).basename;
    for(int chg = charge_start; chg <= charge_end; chg ++)    {
        PSMInfo psminfo;
        if(getPSMInfo(scan, chg, psminfo, name)) {
            psms.push_back(psminfo);
            found = true;
        }
    }
    if(psms.size()>1) {
        cout << "Error: Unexpected number of PSMs "<< psms.size() << "  matched to one scan: " << scan << " of filename: " << filename << endl;
    }
    return found;
}

string PeptideProphetParser::generateSpectrumName(int scan, int chg, const string &basename) const {
    string scanstr = to_string(scan);
    int digit_len = 5; // this is a magic number.
    if (digit_len > scanstr.length()) {
        scanstr = string(digit_len - scanstr.length(), '0') + scanstr;
    }
    string spectrumName = basename + "." + scanstr + "." + scanstr + "." + to_string(chg); // how about the zeros
    return spectrumName;
}

bool PeptideProphetParser::updateGtInfo(SPsmAnnotation &gtinfo) {
    bool updated = false;
    string onlyname = File::CFile(gtinfo.mzxml_filename).basename;

    const int MAX_CHARGE = 7;
    for (int chg = 1; chg <= MAX_CHARGE; chg++) {
        if (updateGtInfo(onlyname,gtinfo.ms2_scan, chg,gtinfo)) {
            updated = true;
            break;
        }
    }
    return updated;

}

bool PeptideProphetParser::updateGtInfo(string name, int scan, int chg, SPsmAnnotation &gtinfo) {
    int digit_len = 5;
    string scanstr = to_string(scan);

    if (digit_len > scanstr.length()) {
        scanstr = string(digit_len - scanstr.length(), '0') + scanstr;
    }

    string spectrumName = name + "." + scanstr + "." + scanstr + "." + to_string(chg); // how about the zeros
    return  updateGtInfoOnSpectrumName(spectrumName, gtinfo);
}

bool PeptideProphetParser::updateGtInfoOnSpectrumName(string spectrumName, SPsmAnnotation &gtinfo) {
    bool ret = false;
    PSMInfo psminfo;
    if(getPSMInfobySpectrumName(spectrumName,psminfo))    {
        gtinfo.precursorMz = psminfo.precursorNeutMass;  // replacing gtinfo.precursor mass to psminfo precursor mass...
        gtinfo.charge = psminfo.charge;
        gtinfo.significance = isPSMSignificant(psminfo) ? 1 : 0;
        if(psminfo.searchhits.size()>0)   {
            SearchHit & sh0 = *psminfo.searchhits[0];
            gtinfo.peptideseq = sh0.m_peptide;
            gtinfo.score = sh0.m_expect;
            gtinfo.nterm_mass = sh0.m_nterm_mass;
            gtinfo.cterm_mass = sh0.m_cterm_mass;
            gtinfo.modificationstr = CPosMass(sh0.m_pos_mass).toString();

            gtinfo.pProb =sh0.m_peptideprophet_prob;
            gtinfo.iProb = sh0.m_iprophet_prob;
            gtinfo.protein = psminfo.getProtein_UseAlterProteinIfItsNotDecoy(true, 0);
            if(gtinfo.protein.find("DECOY") != string::npos)   {
                gtinfo.isDecoy = 1;
            } else{
                gtinfo.isDecoy = 0;
            }
            ret = true;
        }
    }
    else{
        // xinteract output file always fill scan number with '0' to get fixed width 5.
//            cout << spectrumName << " does not exist in file" << endl;
        gtinfo.peptideseq = "UNKNOWN";
        gtinfo.score = -1;
        gtinfo.nterm_mass = -1;// todo: should be zero
        gtinfo.cterm_mass = -1; // should be 0 todo:
        gtinfo.precursorMz = -1;
        gtinfo.charge = 0;
        gtinfo.pProb = 0;
        gtinfo.iProb = 0;

    }
    return ret;
}

int PeptideProphetParser::getPsmNum() {
    return m_psm.size();
}

vector<string> PeptideProphetParser::getAllSoruceFiles() {
    return m_allSourceFiles;
}

pair<int, int> PeptideProphetParser::getIndexRange(string filename) {
    return m_filename2indeRange[filename];
}


mzMLReader::mzMLReader(string filename){
    m_filename = filename; // keep an eye on m_filename in other parsers.
    cout << "[Info] Loading file " << filename << endl;
    FILE *pfile = fopen(filename.c_str(), "r");
    if (pfile == nullptr) {
        cout << "[Info] Fail to open file " << filename << endl;
        cout << "[Info] Please make sure the file does exist." << endl;
        cout << "[Info] Program will EXIT! " << endl;
        exit(-1);
    }
    fseek(pfile, 0L, SEEK_END);
    long len = ftell(pfile);
    rewind(pfile);
    m_buf = new char[len + 1];
    if (m_buf == nullptr) {
        cout << "Fail to get pointer" << endl;
        exit(0);
    }
    long nsize = fread(m_buf, sizeof(char), len, pfile);
    m_buf[len] = '\0';
    fclose(pfile);
    doc.parse<0>(m_buf);
    get_all_spectra();
}

mzMLReader::~mzMLReader() {
    cout << "Release buffer" << endl;
    delete[] m_buf;
}

// todo: to be continued from here; this is the core of the function
void mzMLReader::get_all_spectra() {
    m_currentNode = breadth_first_search("spectrum", &doc, false);
    vector<xml_node<> *> spectra = find_all_siblings("spectrum", m_currentNode);

    cout << "Number of spectra found : " << spectra.size() << endl;
    Progress ps(spectra.size(), "Loading spectra");
    for (auto x: spectra) {
        ps.increase();
        Spectra.emplace_back(x);
    }
}

Spectrum::Spectrum(const Spectrum &other) {
    pkl = new PeakList(*other.pkl);
    m_ms_level = other.m_ms_level;
    m_RT_in_seconds = other.m_RT_in_seconds;
    m_precursorcharge = other.m_precursorcharge;
    m_precursormz = other.m_precursormz;
    spectrum_name = other.spectrum_name;
    spectrum_scan = other.spectrum_scan;
}

// this is something helpful for parsing mzML
//
//MIN_REQ = [
//#
//#!NOTE!   exact names will be extracted of current OBO File, comments are just an orientation
//#         pymzml comes with a little script (queryOBO.py) to query the obo file
//#
//#         $ ./example_scripts/queryOBO.py "scan time"
//#         MS:1000016
//#         scan time
//#         "The time taken for an acquisition by scanning analyzers." [PSI:MS]
//#         Is a: MS:1000503 ! scan attribute
//#
//('MS:1000016',['value']             ), #"scan time"
//# -> Could also be ['value','unitName'] to retrieve a
//# tuple of time and unit by calling spectrum['scan time']
//('MS:1000040',['value']             ), #"m/z"
//('MS:1000041',['value']             ), #"charge state"
//('MS:1000127',['name']              ), #"centroid spectrum"
//('MS:1000128',['name']              ), #"profile spectrum"
//('MS:1000133',['name']              ), #"collision-induced dissociation"
//('MS:1000285',['value']             ), #"total ion current"
//('MS:1000422',['name']              ), #"high-energy collision-induced dissociation"
//('MS:1000511',['value']             ), #"ms level"
//('MS:1000512',['value']             ), #"filter string"
//('MS:1000514',['name']              ), #"m/z array"
//('MS:1000515',['name']              ), #"intensity array"
//('MS:1000521',['name']              ), #"32-bit float"
//('MS:1000523',['name']              ), #"64-bit float"
//('MS:1000744',['value']             ), # legacy precursor mz value ...
//]

// this mzML reader should be tested carefully with different data sets.

Spectrum::Spectrum(xml_node<> *spec_node) {
    m_precursormz = -1;
    m_precursorcharge = -1;
    m_ms_level = -1;
    m_RT_in_seconds = -1;
    pkl = new PeakList();
    int peaknum = atoi(spec_node->first_attribute("defaultArrayLength")->value());

    spectrum_name = spec_node->first_attribute("id")->value();
    int pos = spectrum_name.find("scan="); // to be fixed
    if (string::npos == pos) {
//        cout << "Warning: could not find scan number, Making scan with index + 1!" << endl;
        int index = atoi(spec_node->first_attribute("index")->value());
        spectrum_scan = index + 1;
    } else {
        string scanstr = spectrum_name.substr(pos + 5);
        spectrum_scan = stoi(scanstr.c_str(), nullptr);
    }
    parse_ms_level(spec_node);
    parse_retention_time(spec_node);
    parse_peak_list(spec_node, peaknum);
    parse_precursorinfo(spec_node);
}

void Spectrum::parse_precursorinfo( xml_node<> *spec_node) {
    if (m_ms_level == 2) {
        xml_node<> *precursorinfo = XMLParser::breadth_first_search("selectedIon", spec_node, false);
        if (precursorinfo == nullptr) {
            cout << "No precursor infomation for this ms2" << endl;
        } else {
            xml_node<> *precursorParam = XMLParser::breadth_first_search("cvParam", precursorinfo, false);
            vector<xml_node<> *> allParam = find_all_siblings("cvParam", precursorParam);
            for (auto x: allParam) {
                string accession = x->first_attribute("accession")->value();
                if (strcmp(accession.c_str(), "MS:1000744") == 0)  {
                    m_precursormz = atof(x->first_attribute("value")->value());
                } else if (strcmp(accession.c_str(), "MS:1000041") == 0)    {
                    m_precursorcharge = atoi(x->first_attribute("value")->value());
                }
            }
        }
    }
}
//#include "Util.h"
// bugfixed: strncpy to memcpy. as strncpy stop with null charactor
void Spectrum::parse_peak_list(const xml_node<> *spec_node, const int peaknum) {
    xml_node<> *peaknode = spec_node->first_node("binaryDataArrayList");

    xml_node<> *firstarray = peaknode->first_node("binaryDataArray");
    vector<xml_node<> *> allarray = find_all_siblings("binaryDataArray", firstarray);
    for (auto x: allarray) {

        bool zlib = false;
        string valuename = "mz";
        int bits = 32;
        parse_binary_list_param(x, zlib, valuename, bits);

        xml_node<> *binaryinfo = x->first_node("binary");

        string out;
        macaron::Base64::Decode(binaryinfo->value(), out);

//        cout << "parsing mzML: bits: " << bits << " mz: " << valuename << " zlib: " << zlib << " " << peaknum<< endl;
        if (bits == 32) {
            float *values = new float[peaknum];
            unsigned long length = peaknum * sizeof(float);
            if (zlib) {
                uncompress((unsigned char *) values, &length, out.c_str(), out.size());
            } else {
                memcpy((unsigned char *) values, out.c_str(),out.size());
            }
            vector<double> v(values, values + peaknum);
            if (valuename == "mz") {
                pkl->setM_mzList(v);
            } else {
                pkl->setM_intensityList(v);
            }
            delete[] values;

        } else if (bits == 64) {

//            cout << "64 bits version never tested! Use it at your own risk. " << sizeof(double)*8 << "bits peaknum=" << peaknum << endl;
            assert(peaknum!=0);
            double *values =nullptr;
            try{
                values = new double[peaknum];

            } catch (exception &ex) {
                cout << "caught error: " << ex.what() << endl;
            }

            if(nullptr==values) throw runtime_error("Fail to get vector");
            long length = peaknum * sizeof(double);
            if (zlib) {
                uncompress((unsigned char *) values, &length, out.c_str(), out.size());
            } else {
                memcpy((unsigned char *) values, out.c_str(),out.size());
            }

            vector<double> v(peaknum, 0);
            v.assign(values, values + peaknum);
            if (valuename == "mz") {
                pkl->setM_mzList(v);
            } else {
                pkl->setM_intensityList(v);
            }
            delete[] values;
        }
    }
}

void Spectrum::parse_binary_list_param(const xml_node<> *x, bool &zlib, string &valuename, int &bits) {
    zlib = false;
    valuename = "mz";
    xml_node<> *oneParam = x->first_node("cvParam");
    vector<xml_node<> *> allparams = find_all_siblings("cvParam", oneParam);
    for (auto y: allparams) {
        string accession = y->first_attribute("accession")->value();
        if (strcmp(accession.c_str(), "MS:1000521") == 0) {
            bits = 32;
        } else if (strcmp(accession.c_str(), "MS:1000523") == 0) {
            bits = 64;
        } else if (strcmp(accession.c_str(), "MS:1000574") == 0) {
            zlib = true;
        } else if (strcmp(accession.c_str(), "MS:1000514") == 0) {
            valuename = "mz";
        } else if (strcmp(accession.c_str(), "MS:1000515") == 0) {
            valuename = "intensity";
        }
    }
}

void Spectrum::parse_retention_time(const xml_node<> *spec_node) {
    xml_node<> *scaninfo = spec_node->first_node("scanList")->first_node("scan");
    if (scaninfo == nullptr) {
        cout << "no scan info found" << endl;
    } else {
        xml_node<> *scanparam = scaninfo->first_node("cvParam");
        vector<xml_node<> *> allscanparam = find_all_siblings("cvParam", scanparam);
        double retention_time = 0;
        for (auto eachscanparam: allscanparam) {
            string accession = eachscanparam->first_attribute("accession")->value();
            if (strcmp(accession.c_str(), "MS:1000016") == 0) {
                retention_time = atof(eachscanparam->first_attribute("value")->value());
                string unit_RT = eachscanparam->first_attribute("unitAccession")->value();
                if (strcmp(unit_RT.c_str(), "UO:0000031") == 0) {
                    retention_time *= 60; // minutes to seconds
                }
            }
        }
        m_RT_in_seconds = retention_time;
        pkl->setRTinSeconds(m_RT_in_seconds);
    }
}

void Spectrum::parse_ms_level(const xml_node<> *spec_node) {
    xml_node<> *specparams = spec_node->first_node("cvParam");
    vector<xml_node<> *> allspecparams = find_all_siblings("cvParam", specparams);
    for (auto eachparam: allspecparams) {
        string accession = eachparam->first_attribute("accession")->value();
        if (strcmp(accession.c_str(), "MS:1000511") == 0) {
            m_ms_level = atoi(eachparam->first_attribute("value")->value());
        } else {
        }
    }
}

Spectrum::~Spectrum() {delete pkl;}

double get_peptidePropeht_prob(SearchHit &sh) {
    return sh.m_peptideprophet_prob;
}

double get_iprophet_prob(SearchHit &sh) {
    return sh.m_iprophet_prob;
}

// and this is so wrong!!!??
void extract_pep_prob(vector<double> &tProbs, vector<double> &dProbs, PeptideProphetParser &ppp,
                      double (*getscore)(SearchHit &), bool useAlternativeProt) {
    int num = ppp.getPsmNum();
    cout << "total number of PSMs to be extracted: " << num << endl;
    for (int i = 0; i < num; i++) {
        PSMInfo psm;
        ppp.getPSMInfobyindex(i, psm);
        if (psm.searchhits.size() == 0) {
            cout << "empty search hits, skip this one" << endl;
            psm.print();
            continue;
        }

        double prob = getscore(*psm.searchhits[0]);

        if (psm.isDecoy(useAlternativeProt, 0)) {
            dProbs.push_back(prob);
        } else {
            tProbs.push_back(prob);
        }
    }
    cout << "Done data extraction targetScore: " << tProbs.size() << " decoy score " << dProbs.size() << endl;
}

// the pointer should be released by the caller.
ICPepXMLParser *createPepXML(string filename) {
    if(filename.find(".pep.xml")==string::npos and filename.find(".pepXML")==string::npos){
        cout << "file with wrong format. " << filename << endl;
        return nullptr;
    } else{
        if(filename.find("interact-")==string::npos){
            // not by xinteract
            cout << "COMET PepXML detected" << endl;
            return new  CometPepXMLParser(filename);
        } else{
            cout << "XINTERACT PepXML detected" << endl;
            return new PeptideProphetParser(filename);
        }
    }
}

// the two function below are the same.
// todo: merge into one
//void export_peptideprophetinfo(string filename, PeptideProphetParser &ppp, vector<tuple<double, double>> &fdr_counts) {
//    vector<double> tProbs, dProbs;
//    extract_pep_prob(tProbs, dProbs, ppp, get_peptidePropeht_prob);
//
//    vector<tuple<double, double>> FDR_CorrectNum;
//    get_FDR_CorrectNum(tProbs, dProbs, FDR_CorrectNum);
//    // export to file
//    CFdrNumCorr fdr_num_corr;
//    fdr_num_corr.set(FDR_CorrectNum);
//    string outputfile = filename + "_FDR_CorrectNum_pepprob.txt";
//    fdr_num_corr.saveAs(outputfile);
//
//    fdr_counts = FDR_CorrectNum;
//
//}
//
//void export_iprophetinfo(string filename, PeptideProphetParser &ppp, vector<tuple<double, double>> &fdr_counts) {
//    vector<double> tProbs, dProbs;
//    extract_pep_prob(tProbs, dProbs, ppp, get_iprophet_prob);
//
//    vector<tuple<double, double>> FDR_CorrectNum;
//    get_FDR_CorrectNum(tProbs, dProbs, FDR_CorrectNum);
//
//    string outputfile = filename + "_FDR_CorrectNum_iprob.txt";
//    ofstream fout(outputfile.c_str(), ios::out);
//    for (auto x: FDR_CorrectNum) {
//        fout << get<0>(x) << "\t" << get<1>(x) << endl;
//    }
//    fout.close();
//    cout << "Data exported to file : " << outputfile << endl;
//    fdr_counts = FDR_CorrectNum;
//
//}
