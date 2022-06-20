#include "SpectraSTPepXMLSearchOutput.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"

#ifndef STANDALONE_LINUX
#include "Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#endif

#include <time.h>

#ifdef _MSC_VER // MSVC
#include <direct.h>
#else

#include <unistd.h>

#endif

#include <sstream>

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

/* Class: SpectraSTPepXMLSearchOutput
 * 
 * Subclass of SpectraSTSearchOutput that prints to a pepXML file.
 * 
 * (STILL NEED A FEW THINGS TO WORK PROPERLY -- especially recording where the
 * spectra files are for the visualization CGI's in the TPP INTERACT interface)
 * 
 */


// constructor
SpectraSTPepXMLSearchOutput::SpectraSTPepXMLSearchOutput(string outFullFileName, string inputFileExt,
                                                         SpectraSTSearchParams &searchParams, string searchFilePath) :
        SpectraSTSearchOutput(outFullFileName, inputFileExt, searchParams),
        m_outFullFileName(outFullFileName),
        m_searchFilePath(searchFilePath),
        m_count(0) {

    // for pepXML format, NO_MATCH has to be excluded so as not to mess up downstream processing
    m_searchParams.hitListExcludeNoMatch = true;

#ifndef STANDALONE_LINUX
    m_proteolyticEnzyme = NULL;
#endif

}

// destructor
SpectraSTPepXMLSearchOutput::~SpectraSTPepXMLSearchOutput() {

#ifndef STANDALONE_LINUX
    if (m_proteolyticEnzyme) delete (m_proteolyticEnzyme);
#endif

}

// printHeader - prints all the information before the spectrum_query elements
void SpectraSTPepXMLSearchOutput::printHeader() {

    if (!m_fout) openFile();

    FileName fn;
    parseFileName(m_outFullFileName, fn);


    // almost the same as Sequest2XML here. uses the same #define's.
    (*m_fout) << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    (*m_fout) << "<?xml-stylesheet type=\"text/xsl\" href=\"" << PEPXML_STD_XSL << "pepXML_std.xsl" << "\"?>" << endl;

    (*m_fout) << "<msms_pipeline_analysis date=\"" << SpectraSTPepXMLSearchOutput::getDateTime(USE_LOCAL_TIME) << "\" ";
    (*m_fout) << "xmlns=\"" << PEPXML_NAMESPACE << "\" ";
    (*m_fout) << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ";
    (*m_fout) << "xsi:schemaLocation=\"" << PEPXML_NAMESPACE << " " << PEPXML_STD_XSL << PEPXML_SCHEMA << "\" ";
    (*m_fout) << "summary_xml=\"" << fn.name + fn.ext << "\">" << endl;

    // need to get full path -- not relative -- here? why?
    // const int dirSize = 1024;
    // char dirBuf[dirSize];
    // char *ret=getcwd(dirBuf, dirSize);

    string baseName(fn.path);
    if (!(m_searchFilePath.empty())) baseName = m_searchFilePath;

    // baseName += "/";

    // special case: if new .pep.xml extension is used, then fn.name will contain the .pep piece. Take that off
    //  if (fn.name.length() > 4 && fn.name.compare(fn.name.length() - 4, 4, ".pep") == 0) {
    //  baseName += fn.name.substr(0, fn.name.length() - 4);
    //} else {

    baseName += fn.name;

    //}

    (*m_fout) << "<msms_run_summary base_name=\"" << baseName << "\" ";

    // if instrument info is available (which is only the case for .mzXML input), also print them
    if (m_instrInfo && m_instrInfo->m_instrumentStructPtr) {
        (*m_fout) << "msManufacturer=\"" << m_instrInfo->m_instrumentStructPtr->manufacturer << "\" ";
        (*m_fout) << "msModel=\"" << m_instrInfo->m_instrumentStructPtr->model << "\" ";
        (*m_fout) << "msIonization=\"" << m_instrInfo->m_instrumentStructPtr->ionisation << "\" ";
        (*m_fout) << "msMassAnalyzer=\"" << m_instrInfo->m_instrumentStructPtr->analyzer << "\" ";
        (*m_fout) << "msDetector=\"" << m_instrInfo->m_instrumentStructPtr->detector << "\" ";
    }

    // HENRY'S Q: what's the difference between raw_data_type and raw_data? Just put .msp or .mzXML for both now
    // MH: Someone please answer this question. The data types seem redundant, although I see files that list both
    //     .raw and .mzXML. Why both? And how do you know which label to put with which attribute?
    (*m_fout) << "raw_data_type=\"" << m_searchFileExt << "\" ";
    (*m_fout) << "raw_data=\"" << m_searchFileExt << "\">" << endl;

    if (m_searchParams.enzymeForPepXMLOutput.empty()) {
        // default to trypsin if -s_ENZ option is not set.
        // hard-coded sample_enzyme element for trypsin
        (*m_fout) << "<sample_enzyme name=\"trypsin\">" << endl;
        (*m_fout) << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
        (*m_fout) << "</sample_enzyme>" << endl;

    } else {

#ifndef STANDALONE_LINUX
        ProteolyticEnzymeFactory* enzFactory = new ProteolyticEnzymeFactory();
        ProteolyticEnzyme* enz = enzFactory->getProteolyticEnzyme(m_searchParams.enzymeForPepXMLOutput.c_str());

        if (enz) {
          enz->writePepXMLTags(*m_fout);
          m_proteolyticEnzyme = enz;
        }
        delete (enzFactory);
#else
        (*m_fout) << "<sample_enzyme name=\"trypsin\">" << endl;
        (*m_fout) << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
        (*m_fout) << "</sample_enzyme>" << endl;
#endif

    }

    // search_summary element (not quite applicable to SpectraST), but have to be hard-coded
    (*m_fout) << "<search_summary base_name=\"" << baseName << "\" ";
    (*m_fout) << "search_engine=\"" << "SpectraST" << "\" ";
    (*m_fout) << "precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" ";
    (*m_fout) << "out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">" << endl;


    if (!m_searchParams.databaseFile.empty()) {
        (*m_fout) << "<search_database local_path=\"" << m_searchParams.databaseFile << "\" type=\""
                  << m_searchParams.databaseType << "\"/>" << endl;
    }

    // SpectraSTSearchParams can be entered here...
    m_searchParams.printPepXMLSearchParams(*m_fout);

    /*
    string icatType("");
    if (m_searchParams.expectedCysteineMod == "ICAT_cl") {
      icatType = "cl";
    } else if (m_searchParams.expectedCysteineMod == "ICAT_uc") {
      icatType = "uc";
    }

    (*m_fout) << "<parameter name=\"icat_type\" value=\"" << icatType << "\"/>" << endl;
    */

    // etc

    (*m_fout) << "</search_summary>" << endl;


}

// printStartQuery - prints the query information (before the search results)
void SpectraSTPepXMLSearchOutput::printStartQuery(string query, double precursorMz, int assumedCharge,
                                                  double retentionTime) {

    // parse the query name to find the startScan, endScan etc.
    // the query name is of the format: <baseName>.<startScan>.<endScan>.<charge>
    // note that in a typical SpectraST search, the charge is not specified, so it will be set to zero
    string::size_type dotpos;
    dotpos = query.rfind('.');
    dotpos = query.rfind('.', dotpos - 1);
    dotpos = query.rfind('.', dotpos - 1);

    string startScan("0");
    string endScan("0");

    if (dotpos != string::npos) {
        startScan = nextToken(query, dotpos + 1, dotpos, ".");
        endScan = nextToken(query, dotpos + 1, dotpos, ".");

        // now dotpos points to the dot before the "charge", which should be zero. Replace this zero with the
        // assumed charge
        stringstream ss;
        ss << query.substr(0, dotpos + 1) << assumedCharge;
        query = ss.str();

    }

    (*m_fout) << "<spectrum_query spectrum=\"" << query << "\" ";
    (*m_fout) << "start_scan=\"" << startScan << "\" end_scan=\"" << endScan << "\" ";

    // calculate precursor_neutral_mass from precursorMz and assumedCharge
    double protonMass = (*Peptide::AAMonoisotopicMassTable)['+'];
    if (m_searchParams.precursorMzUseAverage) {
        protonMass = (*Peptide::AAAverageMassTable)['+'];
    }

    double precursorNeutralMass = precursorMz * assumedCharge - assumedCharge * protonMass;

    (*m_fout) << "precursor_neutral_mass=\"" << precursorNeutralMass << "\" ";
    (*m_fout) << "assumed_charge=\"" << assumedCharge << "\" ";
    (*m_fout) << "index=\"" << ++m_count << "\"";

    if (retentionTime >= 0.0) {
        m_fout->precision(2);
        (*m_fout) << " retention_time_sec=\"" << fixed << retentionTime << "\"";
    }

    (*m_fout) << ">" << endl;

    // HOPEFULLY WE'LL HAVE A NEW PEPXML SCHEMA FOR SPECTRAST SOON

    (*m_fout) << "<search_result>" << endl;


}

// printEndQuery - prints whatever is needed to finish up a spectrum_query
void SpectraSTPepXMLSearchOutput::printEndQuery(string query) {

    (*m_fout) << "</search_result>" << endl;
    (*m_fout) << "</spectrum_query>" << endl;


}

// printHit - prints a search_hit
void SpectraSTPepXMLSearchOutput::printHit(string query, unsigned int hitRank, SpectraSTLibEntry *entry,
                                           SpectraSTSimScores &simScores) {

    if (!entry) {
        // for pepXML format, don't display NO_MATCH (since it'll really screw up interact)
        return;
    }

    Peptide *p = entry->getPeptidePtr();

    if (!p) {
        // non-peptide

        double thMass = entry->getPrecursorMz() * entry->getCharge() - entry->getCharge();
        double massDiff = simScores.precursorMzDiff * entry->getCharge();

        (*m_fout) << "<search_hit hit_rank=\"" << hitRank << "\" peptide=\"" << "XNNPEPTIDEX" << "\" ";
        (*m_fout) << "peptide_prev_aa=\"" << "X" << "\" peptide_next_aa=\"" << "X" << "\" ";
        (*m_fout) << "protein=\"" << entry->getSafeName() << "\" num_tot_proteins=\"" << 1 << "\" ";
        (*m_fout) << "calc_neutral_pep_mass=\"";
        (*m_fout).precision(4);
        (*m_fout) << fixed << thMass;
        (*m_fout) << "\" massdiff=\"";
        (*m_fout).precision(4);
        (*m_fout) << fixed << showpoint << showpos << massDiff;
        (*m_fout) << noshowpos << "\" num_tol_term=\"2\" ";
        (*m_fout) << "num_missed_cleavages=\"0\">" << endl;


    } else {

        // parse out the Protein and precursorMH fields from the comments
        // note that for now, the library only stores one protein per entry (no alternative proteins for shared peptides)
        string protein("");
        entry->getOneComment("Protein", protein); // will get "" if Protein is not found

        vector<string> proteins;

        string::size_type slashPos = protein.find('/', 0);

        if (slashPos == string::npos || slashPos >= protein.length() - 1) {
            // single protein, old format
            proteins.push_back(protein);

        } else {

            unsigned int numProteins = atoi(protein.substr(0, slashPos).c_str());
            // parse out all proteins
            string nextProtein("");
            while (!((nextProtein = nextToken(protein, slashPos + 1, slashPos, "/\t\r\n")).empty())) {
                proteins.push_back(nextProtein);
            }

            if (numProteins != (unsigned int) (proteins.size())) {
                // something wrong... but don't complain
                // cerr << "Protein count doesn't match number of listed proteins!" << endl;
            }
        }

        // calculate the massdiff
        double thMass = 0.0;
        if (m_searchParams.precursorMzUseAverage) {
            thMass = p->averageNeutralM();
        } else {
            thMass = p->monoisotopicNeutralM();
        }
        double massDiff = simScores.precursorMzDiff * entry->getCharge();

        (*m_fout) << "<search_hit hit_rank=\"" << hitRank << "\" peptide=\"" << p->stripped << "\" ";
        (*m_fout) << "peptide_prev_aa=\"" << p->prevAA << "\" peptide_next_aa=\"" << p->nextAA << "\" ";
        (*m_fout) << "protein=\"" << proteins[0] << "\" num_tot_proteins=\"" << proteins.size() << "\" ";
        (*m_fout) << "calc_neutral_pep_mass=\"";
        (*m_fout).precision(4);
        (*m_fout) << fixed << thMass;
        (*m_fout) << "\" massdiff=\"";
        (*m_fout).precision(4);
        (*m_fout) << fixed << showpoint << showpos << massDiff;


        unsigned int ntt = 2;
        unsigned int nmc = 0;

#ifndef STANDALONE_LINUX
        if (m_proteolyticEnzyme) {
          ntt = m_proteolyticEnzyme->getNumTolTerm(p->prevAA, p->stripped.c_str(), p->nextAA);
          nmc = m_proteolyticEnzyme->getNumMissedCleavages(p->stripped.c_str(), NULL);	// assuming no cleavage site will become a non-site if modified
        } else {
          ntt = p->NTT();
          nmc = p->NMC();
        }
#else
        ntt = p->NTT();
        nmc = p->NMC();
#endif

        (*m_fout) << noshowpos << "\" num_tol_term=\"" << ntt;
        (*m_fout) << "\" num_missed_cleavages=\"" << nmc << "\">" << endl;

        if (proteins.size() > 1) {
            for (vector<string>::size_type i = 1; i < proteins.size(); i++) {
                (*m_fout) << "<alternative_protein protein=\"" << proteins[i] << "\"/>" << endl;
                // no protein desc, num_tol_term or prev/next AA, schema doesn't require them. but will this break something?
            }
        }

        if (!p->mods.empty() || !p->nTermMod.empty() || !p->cTermMod.empty()) {
            // peptide is modified

            (*m_fout) << "<modification_info ";
            if (!p->nTermMod.empty()) {
                if (m_searchParams.precursorMzUseAverage) {
                    (*m_fout) << "mod_nterm_mass=\"" << Peptide::getAAPlusModAverageMass('n', p->nTermMod) << "\" ";
                } else {
                    (*m_fout) << "mod_nterm_mass=\"" << Peptide::getAAPlusModMonoisotopicMass('n', p->nTermMod)
                              << "\" ";
                }
            }
            if (!p->cTermMod.empty()) {
                if (m_searchParams.precursorMzUseAverage) {
                    (*m_fout) << "mod_cterm_mass=\"" << Peptide::getAAPlusModAverageMass('c', p->cTermMod) << "\" ";
                } else {
                    (*m_fout) << "mod_cterm_mass=\"" << Peptide::getAAPlusModMonoisotopicMass('c', p->cTermMod)
                              << "\" ";
                }
            }

            (*m_fout) << "modified_peptide=\"" << p->interactStyle() << "\">" << endl;

            map<int, string>::iterator i;
            for (i = p->mods.begin(); i != p->mods.end(); i++) {
                // plus 1 because we want the mod position to be one-based; in the Peptide class it is stored as zero-based
                (*m_fout) << "<mod_aminoacid_mass position=\"" << (*i).first + 1;

                if (m_searchParams.precursorMzUseAverage) {
                    (*m_fout) << "\" mass=\"" << Peptide::getAAPlusModAverageMass(p->stripped[(*i).first], (*i).second)
                              << "\"/>" << endl;
                } else {
                    (*m_fout) << "\" mass=\""
                              << Peptide::getAAPlusModMonoisotopicMass(p->stripped[(*i).first], (*i).second) << "\"/>"
                              << endl;
                }
            }
            (*m_fout) << "</modification_info>" << endl;
        }

    }

    // call simScores's output method (specialized for pepXML) to print the sim scores.
    simScores.printPepXML((*m_fout));


    (*m_fout) << "<search_score name=\"charge\" value=\"" << entry->getCharge() << "\"/>" << endl;
    (*m_fout) << "<search_score name=\"lib_file_offset\" value=\"" << entry->getLibFileOffset() << "\"/>" << endl;

    double libProb = entry->getProb(1.0);

    (*m_fout).precision(4);
    (*m_fout) << "<search_score name=\"lib_probability\" value=\"" << fixed << libProb << "\"/>" << endl;

    (*m_fout) << "<search_score name=\"lib_status\" value=\"" << entry->getStatus() << "\"/>" << endl;

    unsigned int numUsed = entry->getNrepsUsed();

    (*m_fout) << "<search_score name=\"lib_num_replicates\" value=\"" << numUsed << "\"/>" << endl;

    string remarkStr("");
    if (entry->getOneComment("Remark", remarkStr)) {
        (*m_fout) << "<search_score name=\"lib_remark\" value=\"" << remarkStr << "\"/>" << endl;
    } else {
        (*m_fout) << "<search_score name=\"lib_remark\" value=\"" << "_NONE_" << "\"/>" << endl;
    }

    (*m_fout) << "</search_hit>" << endl;


}

// printFooter - prints whatever is needed at the end of the pepXML file
void SpectraSTPepXMLSearchOutput::printFooter() {

    (*m_fout) << "</msms_run_summary>" << endl;
    (*m_fout) << "</msms_pipeline_analysis>" << endl;

}

// getDateTime - generates the date/time string.
string SpectraSTPepXMLSearchOutput::getDateTime(bool local) {
    // copied from TPP: Parser::getDateTime

    time_t now;
    time(&now);
    struct tm *tmstruct = NULL;
    if (local) {
        tmstruct = localtime(&now);
    } else {
        tmstruct = gmtime(&now);
    }
    stringstream ss;
    ss.width(4);
    ss << (tmstruct->tm_year + 1900) << '-';
    ss.width(2);
    ss.fill('0');
    ss << (tmstruct->tm_mon + 1) << '-';
    ss.width(2);
    ss.fill('0');
    ss << (tmstruct->tm_mday) << 'T';
    ss.width(2);
    ss.fill('0');
    ss << (tmstruct->tm_hour) << ':';
    ss.width(2);
    ss.fill('0');
    ss << (tmstruct->tm_min) << ':';
    ss.width(2);
    ss.fill('0');
    ss << (tmstruct->tm_sec);

    return (ss.str());
}



