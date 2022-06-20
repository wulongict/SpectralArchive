//
// Created by wulong on 4/12/17.
//
// These features will be used for evaluation of the library entry

#ifndef MYTOOL_PSM_FEATURE_H
#define MYTOOL_PSM_FEATURE_H

#include <string>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <regex>

using namespace std;
class SpectraSTPeakList;
class Peptide;



// this function is not safe, because we are using pointers,
void parse_annotation(string annotation, map<char, vector<int>> &ion_position_map);

class PSMFeature {

    SpectraSTPeakList *m_pSPL;
    Peptide *m_pPep;
    int m_PSMID;

    map<char, int> m_AAcounts;
    double m_Explained_Intensity;
    vector<bool> m_is_B_ion;
    vector<bool> m_is_Y_ion;
    double m_Longest_Ion_Series;
    double m_Longest_b_Ion_Series;
    double m_Longest_y_ion_Series;

    double m_Longest_chg1_Ion_Series;
    double m_Longest_b_chg1_Ion_Series;
    double m_Longest_y_chg1_ion_Series;

    double m_Num_Peaks;
    double m_Pep_Len;
    double m_Num_B_ion;
    double m_Num_Y_ion;
    double m_Num_B_or_Y_ion;
    double m_Total_Intensity;
    double m_Matched_B_Intensity;
    double m_Matched_Y_Intensity;
    vector<double> m_Mass_Error_of_Matched_Peaks;
    vector<double> m_Mass_Error_of_Matched_b_Ions;
    vector<double> m_Mass_Error_of_Matched_y_Ions;
    vector<double> m_Std_Mass_Error_of_Matched_Peaks;
    vector<double> m_Abs_Mass_Error_of_Matched_Peaks;
    double m_Delta_Parent_Mass;

    double m_Std_Mass_Error_of_Matched_Peaks_Normalized;
    double m_Abs_Mass_Error_of_Matched_Peak_Normalized;

    vector<double> m_Exaplained_Cleavages;
    vector<double> m_Exaplained_b_Cleavages;
    vector<double> m_Exaplained_y_Cleavages;

    vector<double> m_Exaplained_chg1_Cleavages;
    vector<double> m_Exaplained_b_chg1_Cleavages;
    vector<double> m_Exaplained_y_chg1_Cleavages;

    double m_Num_Explained_b_chg1_Cleavages;
    double m_Num_Explained_y_chg1_Cleavages;

private:
    double m_Num_Explained_chg1_Cleavages;
    double m_Num_Complementry_b_y_Cleavages;
    double m_Num_Complementry_b_y_chg1_Cleavages;
    double m_Num_Explained_b_Cleavages;
    double m_Num_Explained_y_Cleavages;
    double m_Num_Explained_Cleavages;
    double m_Num_Unassigned_Peaks;
    double m_Intensity_Unassigned_Peaks;
    string m_sourcefilename;


public:
    bool isHighMassAcc();
    double get_MassError_Std_Normalized() const{ return m_Std_Mass_Error_of_Matched_Peaks_Normalized;}
    double get_MassError_AbsAvg_Normalized()const{return m_Abs_Mass_Error_of_Matched_Peak_Normalized;}
    double getM_Longest_chg1_Ion_Series() const;
    double getM_Longest_b_chg1_Ion_Series() const;
    double getM_Longest_y_chg1_ion_Series() const;
    double getM_Num_Explained_b_chg1_Cleavages() const;
    double getM_Num_Explained_y_chg1_Cleavages() const;
    double getM_Num_Explained_chg1_Cleavages() const;
    // get one spectra/peaklist and one peptide.
    // annotate and extract the features
    PSMFeature(SpectraSTPeakList *pSPL, Peptide *pPep, int PSMID, string sourcefilename);
    int getCharge();
    string getsourcefilename();

    ~PSMFeature();
    int getPSMID() const;
    SpectraSTPeakList * getSpecPKLPtr();
    Peptide * getPeptidePtr();
    double get_Matched_Y_Intensity() const;
    double get_Matched_B_Intensity() const;
    double get_Explained_Intensity() const;
    double get_Longest_Ion_Series() const;
    double get_Longest_b_Ion_Series() const;
    double get_Longest_y_Ion_Series() const;
    double get_Num_B_Y_Ions() const;
    double get_total_intensity() const;
    double get_Num_AAs(char AA);
    double get_Pep_Len();
    double get_Peak_Num();
    double get_Num_B_Ions() const;
    double get_Num_Explained_Cleavages() const;
    double get_Num_Explained_b_Cleavages() const;
    double get_Num_Explained_y_Cleavages() const;
    double get_Num_Y_Ions() const;
    double get_Unassigned_Peak_Intensity() const;

    // this function finishes all the feature calculation
    void initialize();
    void get_longest_ion_series(vector<double> &ExplainedCleavages, double &Num_Explained_Cleavages,
                                double &Longest_Ion_Series);

private:
    bool m_initialized;
};






#endif //MYTOOL_PSM_FEATURE_H
