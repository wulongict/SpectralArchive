//
// Created by wulong on 3/1/18.
//

#include "ConcretePSMFeatures.h"
#include "CPeakPairsImporter.h"
#include "FragmentationScore.h"
#include "classifier.h"
#include "BasePSMFeature.h"
#include "MGFReader.h"
#include "CPeakPair.h"


void Print_Ion_position_map(map<char, vector<int>> *ion_position_map) {
    for (int i = 0; i < (*ion_position_map)['b'].size(); ++i) {
        cout << "B " << (*ion_position_map)['b'][i] << endl;
    }

    for (int i = 0; i < (*ion_position_map)['y'].size(); ++i) {
        cout << "Y " << (*ion_position_map)['y'][i] << endl;
    }
}

FragmentationPatternScore::FragmentationPatternScore(bool ghost, double minInt, string model_file, string scoretype,
                                                     string binaryPath) {
//        p_mgf = mgf;
    m_ghost = ghost;
    m_minInt = minInt;
    cout << "[Info] Using minIntFC " << m_minInt << endl;
    m_feature_name = "FragmentationPatternScore";
    m_model_file = model_file;
    cout << "[Info] Using fragmentation model file: " << model_file << endl;
    m_scoretype = scoretype;
    m_binaryPath = binaryPath;
    m_lr = make_shared<classifier>(m_model_file,m_binaryPath);
}

double FragmentationPatternScore::calculate_feature(PSMFeature *pPSM) {
    CPeakPairsImporter PPFR;
    vector<CPeakPair> PPF;

    if (pPSM->getPeptidePtr() != nullptr) {
        get_feature_for_one_spectrum(m_ghost, PPF,pPSM->getPSMID() ,pPSM->getPeptidePtr() ,pPSM->getSpecPKLPtr());
    }
    PPFR.swap_as(PPF);
    PPFR.setSourceFileName(pPSM->getsourcefilename());

    CPeakPairsImporter PPFRminInt;
    PPFR.getSubsetOnMinIntenFoldChange(PPFRminInt, m_minInt);

    double finalscore = 0;
    if (PPFR.size() != 0) {
//        classifier lr(m_model_file, m_binaryPath);
//        lr.loadModel(m_model_file);
        FragmentationScore fragScore(PPFRminInt, m_scoretype, *m_lr);
        finalscore = fragScore.getscore(pPSM->getPSMID());
    }
    return finalscore;
}

FragmentationPatternScore::~FragmentationPatternScore() {
}

Longest_b_ion_series_Prop::~Longest_b_ion_series_Prop() {
}

double Longest_b_ion_series_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_b_Ion_Series() / (pPSM->get_Pep_Len() - 1);
}

Longest_b_ion_series_Prop::Longest_b_ion_series_Prop() {
    m_feature_name = "Longest_b_ion_series_Prop";
    eps = 1e-9;
}

Longest_y_ion_sereis_Prop::~Longest_y_ion_sereis_Prop() {
}

double Longest_y_ion_sereis_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_y_Ion_Series() / (pPSM->get_Pep_Len() - 1);

}

Longest_y_ion_sereis_Prop::Longest_y_ion_sereis_Prop() {
    m_feature_name = "Longest_y_ion_series_Prop";
    eps = 1e-9;
}

Longest_ion_series_Prop::~Longest_ion_series_Prop() {
}

double Longest_ion_series_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_Ion_Series() / (pPSM->get_Pep_Len() - 1);

}

Longest_ion_series_Prop::Longest_ion_series_Prop() {
    m_feature_name = "Longest_ion_series_Prop";
    eps = 1e-9;
}

Longest_y_cs1_ion_series_Prop::~Longest_y_cs1_ion_series_Prop() {
}

double Longest_y_cs1_ion_series_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_y_chg1_ion_Series()/(pPSM->get_Pep_Len()-1);
}

Longest_y_cs1_ion_series_Prop::Longest_y_cs1_ion_series_Prop() {
    m_feature_name = "Longest_y_cs1_ion_series_Prop";
    eps = 1e-9;
}

Longest_y_cs1_ion_series::~Longest_y_cs1_ion_series() {
}

double Longest_y_cs1_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_y_chg1_ion_Series();
}

Longest_y_cs1_ion_series::Longest_y_cs1_ion_series() {
    m_feature_name = "Longest_y_cs1_ion_series";
    eps = 1e-9;
}

Longest_y_ion_series::~Longest_y_ion_series() {
}

double Longest_y_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_y_Ion_Series();
}

Longest_y_ion_series::Longest_y_ion_series() {
    m_feature_name = "Longest_y_ion_series";
    eps = 1e-9;
}

Longest_b_cs1_ion_series_Prop::~Longest_b_cs1_ion_series_Prop() {
}

double Longest_b_cs1_ion_series_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_b_chg1_Ion_Series()/(pPSM->get_Pep_Len()-1);
}

Longest_b_cs1_ion_series_Prop::Longest_b_cs1_ion_series_Prop() {
    m_feature_name = "Longest_b_cs1_ion_series_Prop";
    eps = 1e-9;
}

Longest_b_cs1_ion_series::~Longest_b_cs1_ion_series() {
}

double Longest_b_cs1_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_b_chg1_Ion_Series();
}

Longest_b_cs1_ion_series::Longest_b_cs1_ion_series() {
    m_feature_name = "Longest_b_cs1_ion_series";
    eps = 1e-9;
}

Longest_b_ion_series::~Longest_b_ion_series() {
}

double Longest_b_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_b_Ion_Series();
}

Longest_b_ion_series::Longest_b_ion_series() {
    m_feature_name = "Longest_b_ion_series";
    eps = 1e-9;
}

Longest_cs1_ion_series_Prop::~Longest_cs1_ion_series_Prop() {
}

double Longest_cs1_ion_series_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_chg1_Ion_Series()/(pPSM->get_Pep_Len()-1);
}

Longest_cs1_ion_series_Prop::Longest_cs1_ion_series_Prop() {
    m_feature_name = "Longest_cs1_ion_series_Prop";
    eps = 1e-9;
}

Longest_cs1_ion_series::~Longest_cs1_ion_series() {
}

double Longest_cs1_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Longest_chg1_Ion_Series();
}

Longest_cs1_ion_series::Longest_cs1_ion_series() {
    m_feature_name = "Longest_cs1_ion_series";
    eps = 1e-9;
}

Longest_ion_series::~Longest_ion_series() {
}

double Longest_ion_series::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Longest_Ion_Series();
}

Longest_ion_series::Longest_ion_series() {
    m_feature_name = "Longest_ion_series";
    eps = 1e-9;
}

Explained_y_Cleavage_Prop::~Explained_y_Cleavage_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_y_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_y_Cleavages() / (pPSM->get_Pep_Len() - 1);
}

Explained_y_Cleavage_Prop::Explained_y_Cleavage_Prop() {
    m_feature_name = "Explained_y_Cleavage_Prop";
    eps = 1e-9;
}

Explained_b_Cleavage_Prop::~Explained_b_Cleavage_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_b_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_b_Cleavages() / (pPSM->get_Pep_Len() - 1);
}

Explained_b_Cleavage_Prop::Explained_b_Cleavage_Prop() {
    m_feature_name = "Explained_b_Cleavage_Prop";
    eps = 1e-9;
}

Explained_Cleavage_Prop::~Explained_Cleavage_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_Cleavages() / (pPSM->get_Pep_Len() - 1);
}

Explained_Cleavage_Prop::Explained_Cleavage_Prop() {
    m_feature_name = "Explained_Cleavage_Prop";
    eps = 1e-9;
}

Explained_y_cs1_Cleavage_Prop::~Explained_y_cs1_Cleavage_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_y_cs1_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_y_chg1_Cleavages()/(pPSM->get_Pep_Len()-1);
}

Explained_y_cs1_Cleavage_Prop::Explained_y_cs1_Cleavage_Prop() {
    m_feature_name = "Explained_y_cs1_Cleavage_Prop";
    eps = 1e-9;
}

Explained_y_cs1_Cleavage::~Explained_y_cs1_Cleavage() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_y_cs1_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_y_chg1_Cleavages();
}

Explained_y_cs1_Cleavage::Explained_y_cs1_Cleavage() {
    m_feature_name = "Explained_y_cs1_Cleavage";
    eps = 1e-9;
}

Explained_y_Cleavage::~Explained_y_Cleavage() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_y_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_y_Cleavages();
}

Explained_y_Cleavage::Explained_y_Cleavage() {
    m_feature_name = "Explained_y_Cleavage";
    eps = 1e-9;
}

Explained_b_cs1_Cleavage_Prop::~Explained_b_cs1_Cleavage_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_b_cs1_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_b_chg1_Cleavages()/(pPSM->get_Pep_Len()-1);
}

Explained_b_cs1_Cleavage_Prop::Explained_b_cs1_Cleavage_Prop() {
    m_feature_name = "Explained_b_cs1_Cleavage_Prop";
    eps = 1e-9;
}

Explained_b_cs1_Cleavage::~Explained_b_cs1_Cleavage() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_b_cs1_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_b_chg1_Cleavages();
}

Explained_b_cs1_Cleavage::Explained_b_cs1_Cleavage() {
    m_feature_name = "Explained_b_cs1_Cleavage";
    eps = 1e-9;
}

Explained_b_Cleavage::~Explained_b_Cleavage() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Explained_b_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_b_Cleavages();
}

Explained_b_Cleavage::Explained_b_Cleavage() {
    m_feature_name = "Explained_b_Cleavage";
    eps = 1e-9;
}

double Precursor_Charge::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getCharge();
}

Precursor_Charge::Precursor_Charge() {
    m_feature_name = "Precursor_Charge";
}

double Peptide_Length::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Pep_Len();
}

Peptide_Length::Peptide_Length() {
    m_feature_name = "Peptide_Length";
}

double Explained_cs1_Cleavage_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_chg1_Cleavages()/(pPSM->get_Pep_Len()-1);
}

Explained_cs1_Cleavage_Prop::Explained_cs1_Cleavage_Prop() {
    m_feature_name = "Explained_cs1_Cleavage_Prop";
    eps = 1e-9;
}

double Explained_cs1_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->getM_Num_Explained_chg1_Cleavages();
}

Explained_cs1_Cleavage::Explained_cs1_Cleavage() {
    m_feature_name = "Explained_cs1_Cleavage";
    eps = 1e-9;
}

double Explained_Cleavage::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Explained_Cleavages();
}

Explained_Cleavage::Explained_Cleavage() {
    m_feature_name = "Explained_Cleavage";
    eps = 1e-9;
}

Unassigned_Intensity_Prop::~Unassigned_Intensity_Prop() {
    // cout << "Releasing -- " << m_feature_name << endl;
}

double Unassigned_Intensity_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Unassigned_Peak_Intensity() / pPSM->get_total_intensity();
}

Unassigned_Intensity_Prop::Unassigned_Intensity_Prop() {
    m_feature_name = "Unassigned_Intensity_Prop";
}

UnassignedIntensity::~UnassignedIntensity() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double UnassignedIntensity::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Unassigned_Peak_Intensity();
}

UnassignedIntensity::UnassignedIntensity() {
    m_feature_name = "UnassignedIntensity";
}

Total_Intensity::~Total_Intensity() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Total_Intensity::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_total_intensity();
}

Total_Intensity::Total_Intensity() {
    m_feature_name = "Total_Intensity";
}

Matched_b_y_Intensity_Prop::~Matched_b_y_Intensity_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_b_y_Intensity_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Explained_Intensity()/pPSM->get_total_intensity();
}

Matched_b_y_Intensity_Prop::Matched_b_y_Intensity_Prop() {
    m_feature_name = "Matched_b_y_Intensity_Prop";
}

Matched_b_y_Intensity::~Matched_b_y_Intensity() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_b_y_Intensity::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Explained_Intensity();
}

Matched_b_y_Intensity::Matched_b_y_Intensity() {
    m_feature_name = "Matched_b_y_Intensity";
}

Matched_y_Intensity_Prop::~Matched_y_Intensity_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_y_Intensity_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Matched_Y_Intensity()/pPSM->get_total_intensity();
}

Matched_y_Intensity_Prop::Matched_y_Intensity_Prop() {
    m_feature_name = "Matched_y_Intensity_Prop";
}

Matched_y_Intensity::~Matched_y_Intensity() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_y_Intensity::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Matched_Y_Intensity();
}

Matched_y_Intensity::Matched_y_Intensity() {
    m_feature_name = "Matched_y_Intensity";
}

Matched_b_Intensity_Prop::~Matched_b_Intensity_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_b_Intensity_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Matched_B_Intensity()/pPSM->get_total_intensity();
}

Matched_b_Intensity_Prop::Matched_b_Intensity_Prop() {
    m_feature_name = "Matched_b_Intensity_Prop";
}

Matched_b_Intensity::~Matched_b_Intensity() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Matched_b_Intensity::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Matched_B_Intensity();
}

Matched_b_Intensity::Matched_b_Intensity() {
    m_feature_name = "Matched_b_Intensity";
}

Num_Of_Matched_Peaks_Prop_Calibrated::Num_Of_Matched_Peaks_Prop_Calibrated() {
    m_feature_name = "Num_of_Matched_Peaks_Prop_Calibrated";
}

#include <cmath>
double Num_Of_Matched_Peaks_Prop_Calibrated::calculate_feature(PSMFeature *pPSM) {
    double peplen= pPSM->get_Pep_Len()<=0? 1: pPSM->get_Pep_Len();
    double newnorm = 36*pow(peplen/20.0,log(3)/log(2)-1);
    return pPSM->get_Num_B_Y_Ions()/newnorm;
}

Num_Of_Matched_Peaks_Prop_Calibrated::~Num_Of_Matched_Peaks_Prop_Calibrated() {}

Num_Of_Matched_Peaks_Prop::~Num_Of_Matched_Peaks_Prop() {
}

double Num_Of_Matched_Peaks_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_B_Y_Ions()/pPSM->get_Peak_Num();
}

Num_Of_Matched_Peaks_Prop::Num_Of_Matched_Peaks_Prop() {
    m_feature_name = "Num_of_Matched_Peaks_Prop";
}

Num_Of_Matched_Peaks::~Num_Of_Matched_Peaks() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Num_Of_Matched_Peaks::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_B_Y_Ions();
}

Num_Of_Matched_Peaks::Num_Of_Matched_Peaks() {
    m_feature_name = "Num_of_Matched_Peaks";
}

Num_Of_y_ions_Prop::~Num_Of_y_ions_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Num_Of_y_ions_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Y_Ions()/pPSM->get_Peak_Num();
}

Num_Of_y_ions_Prop::Num_Of_y_ions_Prop() {
    m_feature_name = "Num_of_y_ions_Prop";
}

Num_Of_y_ions::~Num_Of_y_ions() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

double Num_Of_y_ions::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_Y_Ions();
}

Num_Of_y_ions::Num_Of_y_ions() {
    m_feature_name = "Num_of_y_ions";
}

Num_Of_b_ions::Num_Of_b_ions() {
    m_feature_name = "Num_of_b_ions";
}

double Num_Of_b_ions::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_B_Ions();
}

Num_Of_b_ions::~Num_Of_b_ions() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

Num_Of_b_ions_Prop::Num_Of_b_ions_Prop() {
    m_feature_name = "Num_of_b_ions_Prop";
}

double Num_Of_b_ions_Prop::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_B_Ions()/(pPSM->get_Peak_Num());
}

Num_Of_b_ions_Prop::~Num_Of_b_ions_Prop() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

Num_of_AA::Num_of_AA(char AA) {
    m_AA = AA;
    m_feature_name = "Num_of_AA_";
    m_feature_name.append(1, AA);
}

double Num_of_AA::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Num_AAs(m_AA);
}

Num_of_AA::~Num_of_AA() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

isChargeX::isChargeX(int charge) {
    m_charge = charge;
    m_feature_name = "is_Charge_"+to_string(m_charge);
}

double isChargeX::calculate_feature(PSMFeature *pPSM) {
    return (m_charge == pPSM->getCharge());
}

isChargeX::~isChargeX() {
    //cout << "Releasing -- " << m_feature_name << endl;
}

Num_of_Peaks::Num_of_Peaks() {
    m_feature_name="Num_of_Peaks";
}

double Num_of_Peaks::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_Peak_Num();
}

Num_of_Peaks::~Num_of_Peaks() {}

Mass_Error_Std::Mass_Error_Std() {
    m_feature_name="Mass_Error_Std";
}

double Mass_Error_Std::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_MassError_Std_Normalized();

}

Mass_Error_Std::~Mass_Error_Std() {}

Mass_Error_AbsAvg::Mass_Error_AbsAvg() {
    m_feature_name="Mass_Error_AbsAvg";
}

double Mass_Error_AbsAvg::calculate_feature(PSMFeature *pPSM) {
    return pPSM->get_MassError_AbsAvg_Normalized();

}

Mass_Error_AbsAvg::~Mass_Error_AbsAvg() {}
