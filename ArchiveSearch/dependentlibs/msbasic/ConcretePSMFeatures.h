#ifndef PSM_FEATURE_EXTRACTION_H
#define PSM_FEATURE_EXTRACTION_H


#include <map>
#include <vector>

using namespace std;
class PSMFeature;
class Feature {
public:
    string m_feature_name;
    double m_min; // the min max value for rescalling of the data values
    double m_max;
    Feature() {
        m_feature_name = "Basic_Feature";
    }
    virtual ~Feature() {
        //cout << "Releasing "  << m_feature_name << endl;
    };
    virtual double calculate_feature(PSMFeature *pPSM) { return 0; }
};

void Print_Ion_position_map(map<char, vector<int>> *ion_position_map);
class Num_Of_b_ions : public Feature {
public:
    Num_Of_b_ions();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_b_ions();
};
class Mass_Error_Std: public Feature{
public:
    Mass_Error_Std();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Mass_Error_Std();
};
class Mass_Error_AbsAvg: public Feature{
public:
    Mass_Error_AbsAvg();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Mass_Error_AbsAvg();
};

class Num_Of_b_ions_Prop : public Feature {
public:
    Num_Of_b_ions_Prop();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_b_ions_Prop();
};

class Num_Of_y_ions : public Feature {
public:
    Num_Of_y_ions();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_y_ions();
};

class Num_Of_y_ions_Prop : public Feature {
public:
    Num_Of_y_ions_Prop();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_y_ions_Prop();
};

class Num_Of_Matched_Peaks : public Feature {//feature_num_peaks_matched
public:
    Num_Of_Matched_Peaks();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_Matched_Peaks();
};

class Num_Of_Matched_Peaks_Prop : public Feature {//feature_num_peaks_matched
public:
    Num_Of_Matched_Peaks_Prop();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_Matched_Peaks_Prop();
};

class Num_of_Peaks : public Feature {//feature_num_peaks_matched
public:
    Num_of_Peaks();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_of_Peaks();
};

class Num_Of_Matched_Peaks_Prop_Calibrated : public Feature {//ratio calculated again
public:
    Num_Of_Matched_Peaks_Prop_Calibrated();
    double calculate_feature(PSMFeature *pPSM);
    virtual ~Num_Of_Matched_Peaks_Prop_Calibrated();
};

class Matched_b_Intensity : public Num_Of_b_ions {
public:
    Matched_b_Intensity();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_b_Intensity();
};

class Matched_b_Intensity_Prop : public Num_Of_b_ions {
public:
    Matched_b_Intensity_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_b_Intensity_Prop();
};

class Matched_y_Intensity : public Num_Of_y_ions {
public:
    Matched_y_Intensity();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_y_Intensity();
};

class Matched_y_Intensity_Prop : public Num_Of_y_ions {
public:
    Matched_y_Intensity_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_y_Intensity_Prop();
};

class Matched_b_y_Intensity : public Num_Of_Matched_Peaks {//feature_intensity_peaks_matched
public:
    Matched_b_y_Intensity();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_b_y_Intensity();
};

class Matched_b_y_Intensity_Prop : public Num_Of_Matched_Peaks {//feature_intensity_peaks_matched
public:
    Matched_b_y_Intensity_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Matched_b_y_Intensity_Prop();
};

class Total_Intensity : public Feature {
public:

    Total_Intensity();
    double calculate_feature(PSMFeature *pPSM);
    ~Total_Intensity();
};
class UnassignedIntensity : public Feature {
public:
    UnassignedIntensity();
    double calculate_feature(PSMFeature *pPSM);
    ~UnassignedIntensity();
};

class Unassigned_Intensity_Prop : public Feature {//feature_intensity_peaks_matched_ratio
public:
    Unassigned_Intensity_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Unassigned_Intensity_Prop();
};
class Explained_Cleavage : public Feature {
    double eps;
public:
    Explained_Cleavage();
    double calculate_feature(PSMFeature *pPSM);};

class Explained_cs1_Cleavage : public Feature {
    double eps;
public:
    Explained_cs1_Cleavage();
    double calculate_feature(PSMFeature *pPSM);
};

class Explained_cs1_Cleavage_Prop : public Feature {
    double eps;
public:
    Explained_cs1_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
};

class Peptide_Length : public Feature {
public:
    Peptide_Length();

    double calculate_feature(PSMFeature *pPSM);
};

class Precursor_Charge : public Feature {
public:
    Precursor_Charge();
    double calculate_feature(PSMFeature *pPSM);

};

class Explained_b_Cleavage : public Feature {
    double eps;
public:
    Explained_b_Cleavage();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_b_Cleavage();
};

class Explained_b_cs1_Cleavage : public Feature {
    double eps;
public:
    Explained_b_cs1_Cleavage();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_b_cs1_Cleavage();
};

class Explained_b_cs1_Cleavage_Prop : public Feature {
    double eps;
public:
    Explained_b_cs1_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_b_cs1_Cleavage_Prop();
};

class Explained_y_Cleavage : public Feature {
    double eps;
public:
    Explained_y_Cleavage();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_y_Cleavage();
};

class Explained_y_cs1_Cleavage : public Feature {
    double eps;
public:
    Explained_y_cs1_Cleavage();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_y_cs1_Cleavage();
};

class Explained_y_cs1_Cleavage_Prop : public Feature {
    double eps;
public:
    Explained_y_cs1_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_y_cs1_Cleavage_Prop();
};

class Explained_Cleavage_Prop : public Feature { // feature_bonds_explained_ratio
    double eps;
public:
    Explained_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_Cleavage_Prop();
};

class Explained_b_Cleavage_Prop : public Feature { // feature_bonds_explained_ratio
    double eps;
public:

    Explained_b_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_b_Cleavage_Prop();
};

class Explained_y_Cleavage_Prop : public Feature { // feature_bonds_explained_ratio
    double eps;
public:
    Explained_y_Cleavage_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Explained_y_Cleavage_Prop();
};
class Longest_ion_series : public Feature {
    double eps; // to be zero
public:

    Longest_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_ion_series();
};

class Longest_cs1_ion_series : public Feature {
    double eps; // to be zero
public:

    Longest_cs1_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_cs1_ion_series();
};

class Longest_cs1_ion_series_Prop : public Feature {
    double eps; // to be zero
public:

    Longest_cs1_ion_series_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_cs1_ion_series_Prop();
};

class Longest_b_ion_series : public Feature {
    double eps; // to be zero
public:
    Longest_b_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_b_ion_series();
};
class Longest_b_cs1_ion_series : public Feature {
    double eps; // to be zero
public:
    Longest_b_cs1_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_b_cs1_ion_series();
};

class Longest_b_cs1_ion_series_Prop : public Feature {
    double eps; // to be zero
public:
    Longest_b_cs1_ion_series_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_b_cs1_ion_series_Prop();
};
class Longest_y_ion_series : public Feature {
    double eps; // to be zero
public:
    Longest_y_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_y_ion_series();
};

class Longest_y_cs1_ion_series : public Feature {
    double eps; // to be zero
public:
    Longest_y_cs1_ion_series();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_y_cs1_ion_series();
};

class Longest_y_cs1_ion_series_Prop : public Feature {
    double eps; // to be zero
public:
    Longest_y_cs1_ion_series_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_y_cs1_ion_series_Prop();
};
class Longest_ion_series_Prop : public Feature {// feature_longest_bond_series_ratio (5)
    double eps; // to be zero
public:
    Longest_ion_series_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_ion_series_Prop();
};

class Longest_y_ion_sereis_Prop : public Feature {// feature_longest_yion_series_ratio (14)
    double eps; // to be zero
public:
    Longest_y_ion_sereis_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_y_ion_sereis_Prop();
};

class Longest_b_ion_series_Prop : public Feature {// feature_longest_bion_series_ratio (13)
    double eps; // to be zero
public:
    Longest_b_ion_series_Prop();
    double calculate_feature(PSMFeature *pPSM);
    ~Longest_b_ion_series_Prop();
};
#include <memory>
class classifier;
class FragmentationPatternScore : public Feature {
    bool m_ghost;
    double m_minInt;
    string m_model_file;
    string m_scoretype;
    string m_binaryPath;
//    string m_mgf_file;
    shared_ptr<classifier> m_lr;
public:
    FragmentationPatternScore(bool ghost, double minInt, string model_file, string scoretype, string binaryPath);
    double calculate_feature(PSMFeature *pPSM);
    ~FragmentationPatternScore();
};
class Num_of_AA : public Feature {
    char m_AA;
public:
    Num_of_AA(char AA = 'A');
    double calculate_feature(PSMFeature *pPSM);
    ~Num_of_AA();
};



//'ischarge1', 'ischarge2', 'ischarge3', 'ischarge4', 'ischarge5', 'ischarge6', 'ischarge7',
class isChargeX : public Feature {
    int m_charge;
public:
    isChargeX(int charge);
    double calculate_feature(PSMFeature *pPSM);
    ~isChargeX();
};
//void get_theoretical_spectrum(string peptide, int charge)
//{
//    predicted_spectra ps;
//    vector<Peaks> *pks = ps.predict(peptide,charge,30);
//
//}
#endif

/*
     * Features used in python
     * ['rank1feature_bonds_explained_ratio(1)', 'rank1feature_intensity_peaks_matched(3)',
                              'rank1feature_intensity_peaks_matched_ratio(4)', 'rank1feature_num_peaks_matched(2)',
                              'rank1feature_bonds_explained_by_b_y_cs1',
                              'rank1feature_bonds_explained_by_b_y_cs1_ratio',
                              'ischarge1'(6), 'ischarge2(7)', 'ischarge3'(8),
                              'ischarge4(9)', 'ischarge5'(10), 'ischarge6'(11),
                              'ischarge7(12)',
                              'rank1feature_longest_bion_series_ratio'(13),
                              'rank1feature_longest_yion_series_ratio'(14),
                              'rank1feature_longest_bond_series_ratio(5)',
                              'rank1feature_complementary_ion_num_ratio'

                              ]
     * */


