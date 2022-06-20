//
// Created by wulong on 4/12/17.
//

#include <regex>
#include <cmath>
#include <utility>
#include "BasePSMFeature.h"
#include "../../../librarymsms/Util.h"
#include "CDebugMode.h"
#include "peak_annotation.h"

#define STANDALONE_LINUX
#include "../../../External/SpectraST/FileUtils.hpp"
#include "../../../External/SpectraST/SpectraSTPeakList.hpp"


void parse_annotation(string annotation, map<char, vector<int> > &ion_position_map) {
    ion_position_map.insert(std::pair<char, vector<int>>('b', vector<int>()));
    ion_position_map.insert(std::pair<char, vector<int>>('y', vector<int>()));
    trim_space_only(annotation);
    vector<string> ptokens;
    split_string(annotation, ptokens,',');
    for (auto & ptoken : ptokens) {
        unsigned long end = 0;
        string ion_name = nextToken(ptoken, 0, end, "/-^i");
        char iontype = ion_name[0];
        if (iontype != 'b' && iontype != 'y') {
            continue;
        }

        ion_name = ion_name.substr(1);
        int pos = stoi(ion_name);
        (ion_position_map)[iontype].push_back(pos);
    }
}


PSMFeature::PSMFeature(SpectraSTPeakList *pSPL, Peptide *pPep, int PSMID, string sourcefilename) {

//    double highestIntensity = 10000.0, noiseIntensityRatio = 0.01; // peak intensity less than 1% are removed
//    double error_tolerance = 0.02;
    m_sourcefilename = std::move(sourcefilename);
    m_pSPL = pSPL;
    m_pPep = pPep;
    m_PSMID = PSMID;
    m_initialized=false;

    initialize();
}

double PSMFeature::get_Matched_Y_Intensity() const {
    return m_Matched_Y_Intensity;
}

void PSMFeature::initialize() {
    if (m_initialized) {
        return;
    }

    // visit Spectra
    int peaknum = get_Peak_Num();
    m_is_B_ion = vector<bool>(peaknum, false);
    m_is_Y_ion = vector<bool>(peaknum, false);

    m_Explained_Intensity = 0;
    m_Total_Intensity = 0;
    m_Intensity_Unassigned_Peaks = 0;

    m_Matched_B_Intensity = 0;
    m_Matched_Y_Intensity = 0;

    m_Num_B_ion = 0;
    m_Num_Y_ion = 0;

    m_Num_B_or_Y_ion = 0;

    m_Num_Unassigned_Peaks = 0;


    std::regex b_with_NL_regex("b[0-9]+/[0-9.]+|b[0-9]+-[0-9.]+/[0-9.]+");
    std::regex y_with_NL_regex("y[0-9]+/[0-9.]+|y[0-9]+-[0-9.]+/[0-9.]+");

    std::regex b_noNL_regex("b[0-9]+/[0-9.]+");
    std::regex y_noNL_regex("y[0-9]+/[0-9.]+");

    std::regex b_noNL_regex_chg2("b[0-9]+/[0-9.]+|b[0-9]+[\^][2]+/[0-9.]+");
    std::regex y_noNL_regex_chg2("y[0-9]+/[0-9.]+|y[0-9]+[\^][2]+/[0-9.]+");
    int charge = m_pSPL->getParentCharge();

    regex b_regex, y_regex;
    // b- y- ions with NL in charge state up to precursor charge state
    b_regex = std::regex("b[0-9]+[-0-9\\.]*[\\^1-" + to_string(charge) + "]*\\/[-0-9\\.]+");
    y_regex = std::regex("y[0-9]+[-0-9\\.]*[\\^1-" + to_string(charge) + "]*\\/[-0-9\\.]+");


    double highestIntensity = 10000.0, noiseIntensityRatio = 0.01; // peak intensity less than 1% are removed

    double error_tolerance=0.5;
    if (isHighMassAcc()){
        error_tolerance = 0.02;
    }
    //double error_tolerance = 0.02;
    for (int i = 0; i < peaknum; i++) {
        Peak peak;

        m_pSPL->getPeak(i, peak);
        m_Total_Intensity += peak.intensity;
    }
    int parentChg = m_pSPL->getParentCharge();

    for (int i = 0; i < peaknum; ++i) {
        Peak peak;
        m_pSPL->getPeak(i, peak);
        bool skip_peak = false;
        if (peak.intensity < highestIntensity * noiseIntensityRatio) { skip_peak = true; }
//        (m_Total_Intensity) += peak.intensity;
        if (CDebugMode::callDebug()->getMdebug()) {
            cout << "Skip: " << skip_peak << "\t" << peak.mz << "\t" << peak.intensity << "\t" << peak.annotation
                 << endl;

        }
        if (skip_peak) continue;

        vector<peak_annotation> pka;
        parse_annotation_as_structX(peak.annotation, pka);

        for (const auto& pa : pka) {
            if (pa.isotopic) continue;
            if (fabs(pa.masserror) > error_tolerance) continue;
            if (pa.ion_NL != "0") continue;
            if (pa.charge > parentChg and parentChg > 1) continue;

            if (pa.ion_base_type == "b") {
                m_is_B_ion[i] = true;
                m_Mass_Error_of_Matched_b_Ions.push_back(pa.masserror);
            }
            else if (pa.ion_base_type == "y") {
                m_is_Y_ion[i] = true;
                m_Mass_Error_of_Matched_y_Ions.push_back(pa.masserror);
            }
        }

        if (m_is_B_ion[i]) {
            m_Num_B_ion += 1;
            if (CDebugMode::callDebug()->getMdebug()) {
                cout << "B-ion counts " << (m_Num_B_ion) << endl;
            }

            m_Matched_B_Intensity += peak.intensity;
        }

        if (m_is_Y_ion[i]) {
            m_Num_Y_ion += 1;
            if (CDebugMode::callDebug()->getMdebug()) {
                cout << "Y-ion counts " << m_Num_Y_ion << endl;
            }
            m_Matched_Y_Intensity += peak.intensity;
        }


        if (m_is_B_ion[i] || m_is_Y_ion[i]) {
            m_Num_B_or_Y_ion += 1;
            m_Explained_Intensity += peak.intensity;

        }

        if (peak.annotation == "?") {
            m_Intensity_Unassigned_Peaks += peak.intensity;
            m_Num_Unassigned_Peaks += 1;
        }
    }
    m_Mass_Error_of_Matched_Peaks.insert(m_Mass_Error_of_Matched_Peaks.end(), m_Mass_Error_of_Matched_b_Ions.begin(), m_Mass_Error_of_Matched_b_Ions.end());
    m_Mass_Error_of_Matched_Peaks.insert(m_Mass_Error_of_Matched_Peaks.end(), m_Mass_Error_of_Matched_y_Ions.begin(), m_Mass_Error_of_Matched_y_Ions.end());

    m_Std_Mass_Error_of_Matched_Peaks_Normalized = statistic::calcstd(m_Mass_Error_of_Matched_Peaks)/error_tolerance;
    m_Abs_Mass_Error_of_Matched_Peak_Normalized = statistic::calcabsavg(m_Mass_Error_of_Matched_Peaks)/error_tolerance;

    //=============OK-----------------------------

    string pepseq = m_pPep->full();  // X.ABCDEFGHIJK.X/z

    pepseq = pepseq.substr(2, pepseq.length() - 4); // remove X. and .X/z, get ABCDEFGHIJK
    int peplen = pepseq.length();
    m_Exaplained_Cleavages.assign(peplen,0);
    m_Exaplained_b_Cleavages.assign(peplen,0);
    m_Exaplained_y_Cleavages.assign(peplen,0);

    m_Num_Explained_Cleavages = 0.0;
    m_Num_Explained_b_Cleavages = 0.0;
    m_Num_Explained_y_Cleavages = 0.0;

    m_Exaplained_chg1_Cleavages.assign(peplen,0);
    m_Exaplained_b_chg1_Cleavages.assign(peplen,0);
    m_Exaplained_y_chg1_Cleavages.assign(peplen,0);

    m_Num_Explained_chg1_Cleavages = 0.0;
    m_Num_Explained_b_chg1_Cleavages = 0.0;
    m_Num_Explained_y_chg1_Cleavages = 0.0;


    // for each peak
    for (int m = 0; m < m_pSPL->getNumPeaks(); ++m) {

        Peak peak;
        m_pSPL->getPeak(m, peak);
        if (peak.intensity < noiseIntensityRatio * highestIntensity) continue;

        vector<peak_annotation> pka;
        parse_annotation_as_structX(peak.annotation, pka);

        for (const auto& pa : pka) {
            if (pa.isotopic) continue;
            if (fabs(pa.masserror) > error_tolerance) continue;
            if (pa.ion_NL != "0") continue;
            if (pa.charge > parentChg and parentChg > 1) continue;

            if ("b" == pa.ion_base_type) // not isotopic
            {
                m_Exaplained_Cleavages[pa.pos] += peak.intensity;
                m_Exaplained_b_Cleavages[pa.pos] += peak.intensity;
                if (CDebugMode::callDebug()->getMdebug()) {
                    cout << "Look peak ! " << peak.mz << "\t" << peak.intensity << endl;
                    pa.print();
                    cout << "B-ion cleavage found " << pa.pos << endl;
                }
                if (pa.charge == 1) {
                    m_Exaplained_chg1_Cleavages[pa.pos] += peak.intensity;
                    m_Exaplained_b_chg1_Cleavages[pa.pos] += peak.intensity;
                }
            } else if ("y" == pa.ion_base_type) {
                m_Exaplained_Cleavages[pepseq.length() - pa.pos] += peak.intensity;
                m_Exaplained_y_Cleavages[pepseq.length() - pa.pos] += peak.intensity;
                if (CDebugMode::callDebug()->getMdebug()) {
                    cout << "Look peak ! " << peak.mz << "\t" << peak.intensity << endl;
                    pa.print();
                    cout << "Y-ion cleavage found " << pepseq.length() - pa.pos << endl;

                }
                if (pa.charge == 1) {
                    m_Exaplained_chg1_Cleavages[pepseq.length() - pa.pos] += peak.intensity;
                    m_Exaplained_y_chg1_Cleavages[pepseq.length() - pa.pos] += peak.intensity;
                }
            }

        }


    }


    m_Longest_Ion_Series = 0.0;
    m_Longest_b_Ion_Series = 0.0;
    m_Longest_y_ion_Series = 0.0;

    m_Longest_chg1_Ion_Series = 0.0;
    m_Longest_b_chg1_Ion_Series = 0.0;
    m_Longest_y_chg1_ion_Series = 0.0;

//        *m_Longest_Ion_Series = 0;
    get_longest_ion_series(m_Exaplained_Cleavages, m_Num_Explained_Cleavages, m_Longest_Ion_Series);
    get_longest_ion_series(m_Exaplained_b_Cleavages, m_Num_Explained_b_Cleavages, m_Longest_b_Ion_Series);
    get_longest_ion_series(m_Exaplained_y_Cleavages, m_Num_Explained_y_Cleavages, m_Longest_y_ion_Series);

    get_longest_ion_series(m_Exaplained_chg1_Cleavages, m_Num_Explained_chg1_Cleavages, m_Longest_chg1_Ion_Series);
    get_longest_ion_series(m_Exaplained_b_chg1_Cleavages, m_Num_Explained_b_chg1_Cleavages,
                           m_Longest_b_chg1_Ion_Series);
    get_longest_ion_series(m_Exaplained_y_chg1_Cleavages, m_Num_Explained_y_chg1_Cleavages,
                           m_Longest_y_chg1_ion_Series);

    //--------------------------------------OK------------------------------------------------

    // for AAs
    try {
        string legalAAs = "ACDEFGHIKLMNOPQRSTVWY";
//    m_AAcounts = map<char, int>;
        for (char &legalAA : legalAAs) {
            (m_AAcounts)[legalAA] = 0;
//            cout << "(*m_AAcounts)[" << legalAAs[l]<< "] = " << (*m_AAcounts)[legalAAs[l]] << endl;
        }

        for (char &l : pepseq) {
            if (legalAAs.find(l) != string::npos) {
                (m_AAcounts)[l] += 1;
            }

        }
    }
    catch (exception &e) {
        cout << "caught error: " << e.what() << endl ;

    }
    m_initialized=true;
}

// for input vector, which is filled with intensities: ExplainedCleavages
// calculate the number of explained cleavages and longest ion series.
void PSMFeature::get_longest_ion_series(vector<double> &ExplainedCleavages, double &Num_Explained_Cleavages,
                                        double &Longest_Ion_Series) {
    double eps = 1e-9;
    double current_length_Ion_Sereies = 0;
    for (int k = 1; k < ExplainedCleavages.size(); ++k) {
        if (ExplainedCleavages[k] > eps) {
            Num_Explained_Cleavages += 1;
            current_length_Ion_Sereies += 1;
            Longest_Ion_Series = Longest_Ion_Series > current_length_Ion_Sereies ? Longest_Ion_Series
                                                                                     : current_length_Ion_Sereies;
        } else {
            current_length_Ion_Sereies = 0;
        }
    }
}

double PSMFeature::get_Matched_B_Intensity() const {
//    if (m_Matched_B_Intensity == nullptr) {
//        initialize();
//    }
    return m_Matched_B_Intensity;
}

double PSMFeature::get_Explained_Intensity() const {
//    if (m_Explained_Intensity == nullptr) {
//        initialize();
//    }
    return m_Explained_Intensity;
}


double PSMFeature::get_Longest_b_Ion_Series() const {

//    if (m_Num_B_ion == nullptr) {
//        initialize();
//    }
    return m_Longest_b_Ion_Series;
}

double PSMFeature::get_Longest_y_Ion_Series() const {

//    if (m_Num_B_ion == nullptr) {
//        initialize();
//    }
    return m_Longest_y_ion_Series;
}


double PSMFeature::get_Longest_Ion_Series() const {

//    if (m_Num_B_ion == nullptr) {
//        initialize();
//    }
    return m_Longest_Ion_Series;
}

double PSMFeature::get_Num_B_Y_Ions() const {
//    if (m_Num_B_or_Y_ion == nullptr) {
//        initialize();
//    }
    return m_Num_B_or_Y_ion;
}

double PSMFeature::getM_Longest_chg1_Ion_Series() const {
    return m_Longest_chg1_Ion_Series;
}

double PSMFeature::getM_Longest_b_chg1_Ion_Series() const {
    return m_Longest_b_chg1_Ion_Series;
}

double PSMFeature::getM_Longest_y_chg1_ion_Series() const {
    return m_Longest_y_chg1_ion_Series;
}

double PSMFeature::getM_Num_Explained_b_chg1_Cleavages() const {
    return m_Num_Explained_b_chg1_Cleavages;
}

double PSMFeature::getM_Num_Explained_y_chg1_Cleavages() const {
    return m_Num_Explained_y_chg1_Cleavages;
}

double PSMFeature::getM_Num_Explained_chg1_Cleavages() const {
    return m_Num_Explained_chg1_Cleavages;
}

double PSMFeature::get_Num_Explained_y_Cleavages() const {
//        if (m_Num_Explained_y_Cleavages == NULL) {
//            initialize();
//        }
    return m_Num_Explained_y_Cleavages;
}

double PSMFeature::get_Num_Explained_b_Cleavages() const {
//        if (m_Num_Explained_b_Cleavages == NULL) {
//            initialize();
//        }
    return m_Num_Explained_b_Cleavages;
}

double PSMFeature::get_Num_Explained_Cleavages() const {
//        if (m_Num_Explained_Cleavages == NULL) {
//            initialize();
//        }
    return m_Num_Explained_Cleavages;
}

double PSMFeature::get_Num_B_Ions() const {
//        if (m_Num_B_ion == NULL) {
//            initialize();
//        }
    return m_Num_B_ion;
}

double PSMFeature::get_Peak_Num() {

    m_Num_Peaks = m_pSPL->getNumPeaks();


    return m_Num_Peaks;
}

double PSMFeature::get_Pep_Len() {


    string pepseq = m_pPep->full();
//        cout << pepseq << endl;
    pepseq = pepseq.substr(2, pepseq.length() - 4);
    m_Pep_Len = (double) pepseq.length();

    return m_Pep_Len;
}

double PSMFeature::get_Num_AAs(char AA) {
//        if (m_AAcounts == NULL) {
//            initialize();
//        }
    return double((m_AAcounts)[AA]);
}

double PSMFeature::get_total_intensity() const {
//        if (m_Total_Intensity == NULL) {
//            initialize();
//        }
    return m_Total_Intensity;

}

double PSMFeature::get_Num_Y_Ions() const {
//        if (m_Num_Y_ion == NULL) {
//            initialize();
//        }
    return m_Num_Y_ion;
}

double PSMFeature::get_Unassigned_Peak_Intensity() const {
//        if (m_Intensity_Unassigned_Peaks == NULL) {
//            initialize();
//        }
    return m_Intensity_Unassigned_Peaks;
}

Peptide *PSMFeature::getPeptidePtr() {
    return m_pPep;
}

SpectraSTPeakList *PSMFeature::getSpecPKLPtr() {
    return m_pSPL;
}

int PSMFeature::getPSMID() const {
    return m_PSMID;
}

PSMFeature::~PSMFeature() = default;

string PSMFeature::getsourcefilename() {return m_sourcefilename;}

int PSMFeature::getCharge() {
    return m_pSPL->getParentCharge();
}

bool PSMFeature::isHighMassAcc() {
    bool ret = false;
    string m_fragType = m_pSPL->getFragType();
    if(m_fragType == "CID-QTOF" || m_fragType == "HCD" || m_fragType == "ETD-HR"){
        ret= true;
    }
    return ret;
}
