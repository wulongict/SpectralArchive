//
// Created by wulong on 4/11/17.
//

#ifndef MYTOOL_MGFREADER_HPP
#define MYTOOL_MGFREADER_HPP



#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <memory>



using namespace std;
class SpectraSTPeakList;
class Peptide;
class Feature;
class Progress;
class PeakInfo;
class CPeakPair;
class CMyMatrix;



// In essence, each spectrum contains a peaklist and corresponding precursor information
// We would work on each spectrum, or a set of spectra
// I guess we have done this for a lot of times
void collect_peak_info_for_single_spectrum(SpectraSTPeakList *pkl, vector<PeakInfo> &pifs);
void collect_charge1_b_y_ion_info(const vector<PeakInfo> &pifs, vector<PeakInfo> &base_type_no_NL_charge1_pifs);
void
get_each_b_y_ion_existence(const vector<PeakInfo> &pifs, vector<bool> &bionFound, vector<bool> &yionFound) ;
void beforePeakPairFeatureExtraction(string &pepstr, int &charge, vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                     vector<bool> &bionFound, vector<bool> &yionFound, Peptide *peptidePtr,
                                     SpectraSTPeakList *specpkl);

void
get_peakinfo_pair(vector<CPeakPair> &peakPairs, int psm_id, string peptideStr, int charge, PeakInfo &pif1,
                  PeakInfo &pif2);

void get_feature_for_by_ion_ghost(int psm_id, int charge, int k, int l, vector<CPeakPair> &peakPairs,
                                  const string &peptideStr,
                                  const vector<PeakInfo> &base_type_no_NL_charge1_pifs, const string& iontype);
void get_feature_of_ghost_peaks(vector<CPeakPair> &PPF, int psm_id, const string &pepstr, int charge,
                                const vector<PeakInfo> &base_type_no_NL_charge1_pifs,
                                const vector<bool> &bionFound,
                                const vector<bool> &yionFound) ;
void get_feature_for_one_spectrum(bool useGhostPeaks, vector<CPeakPair> &peakPairs, int psm_id, Peptide *ipeptide,
                                  SpectraSTPeakList *specpkl);
class MGFReaderX;
class MGFReader {
private:
    std::shared_ptr<MGFReaderX> m_lazy_mgfreader;
    vector<Peptide *> peptides;
    vector<int> m_annot_index;
    bool m_fixMz;
    map<string,vector<int>> m_pep2idx;
    bool m_isLowMassAcc;
    bool m_truthKnown;
    string m_filename;
    void initPep2Idx();
    bool isTruthPeptide(string &plain_peptide, bool verbose);
    string getTitle(int i );
public:
    MGFReader(const string& filename, bool isLowMassAcc, bool is_truth_known);

    void getWorkingIndex(const string& debug_peptide, int scan, vector<int> &workingIndex );

    ~MGFReader();
    int getSpectraNum();

public:
    SpectraSTPeakList * getSpectrumPtr(int i);
    void getSpectra(int spec_id, vector<tuple<double, double, string>> & onespec);
    string getMgfFileName();
    Peptide * getPeptidePtr(int i ){
        return peptides[i];
    }
//    void generate_theoretical_spectra_from_title();
    void annotate_with_title(bool verbosity, bool debug, int scan, const string& debug_peptide, bool isTrainingSet=false);
    void setPeptidePtr(int i, Peptide *p){
        peptides[i] = p;
    }
    void releasePeptidePtrs();
    //void annotation_with_search_result(MGFReader &mgfr, const string &decoyfile, int rank);
    void annotation_with_search_result(const string &decoyfile, int rank, bool debug, int scan, double fdr_threshold,
                                       const string& debug_peptide, bool isTrainingSet=false);

    void refineWorkingIndex(vector<tuple<int, string, int>> &workingindex);

    void print();
    void print_with_index(int i);
    void filter_min_mz(double min_mz);
    void output_new_mgf(bool verbosity);
    void export_fragmentation_pattern(const string& outfile, bool zero_intensity_fragment_ions);
    void output_pepnovo_training_set(bool verbosity = false);

    // todo list, what is the format of the posttrasloated modifications in PepNovo
    // todo What should I do on the pepNovo training step to get the modification in the correct foramt?

    void export_features_with_PSM(vector<Feature *> features, const string& outputfilename, bool debug, int scan, const string& debug_peptide);


    void annotate_one_spectrum(int i, const string &pepseq, int charge);

    void getPeakPairs(bool useGhostPeaks, vector<CPeakPair> &peakPairs);

private:
    void export_feature_for_one_spectrum(bool zero_intensity_fragment_ions, ofstream &fout, int i);
    void export_for_peakinfo_pair(ofstream &fout, int i, const string &pepstr, int charge, const PeakInfo &pif1,
                                  const PeakInfo &pif2) const;
    void export_feature_of_ghost_peaks(ofstream &fout, int i, const string &pepstr, int charge,
                                       const vector<PeakInfo> &base_type_no_NL_charge1_pifs, const vector<bool> &bion,
                                       const vector<bool> &yion) const;

    void export_feature_for_by_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                         const vector<PeakInfo> &base_type_no_NL_charge1_pifs, string iontype) const;
    void export_feature_for_b_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                        const vector<PeakInfo> &base_type_no_NL_charge1_pifs) const;
    void export_feature_for_y_ion_ghost(int i, int charge, int k, int l, ofstream &fout, const string &pepstr,
                                        const vector<PeakInfo> &base_type_no_NL_charge1_pifs) const;
};


#endif //MYTOOL_MGFREADER_HPP
