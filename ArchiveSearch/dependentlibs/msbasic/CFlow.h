//
// Created by wulong on 2/26/18.
//

#ifndef MYTOOL_CFLOW_H
#define MYTOOL_CFLOW_H

#include <memory>
#include <vector>
using namespace std;
class DataFile;
class CSpectrum;
class Feature;
class CPeakPairsImporter;
class PSMInfo;
class PeptideProphetParser;
class CTable;
class SimpleTimer;
class Progress;
class ICPepXMLParser;

class CException : public std::runtime_error{
public:
    CException(const char * msg);
};



class CROCPlot {
    struct SROC{
        double tpr;
        double fpr;
        double score;
        double fdr;
        int fp_num;
        int tp_num;
        SROC(double _fpr, double _tpr, double score_, int fp_num_, int tp_num_){
            tpr=_tpr;
            fpr=_fpr;
            score=score_;
            fp_num=fp_num_;
            tp_num=tp_num_;
            if(fp_num_ == 0) {
                fdr=0;
            }
            else if (tp_num_ == 0){
                fdr=1;
            }
            else {
                fdr= (double)fp_num_ / tp_num_;
            }
        }

        string getHeaderString(){
            return "FPR TPR Score FP TP FDR";
        }

        friend ostream & operator<<(ostream & os, const SROC& roc)
        {
            os << roc.fpr << " " << roc.tpr << " " << roc.score << " " << roc.fp_num << " " << roc.tp_num << " " << roc.fdr;
            return os;
        }


    };
    vector<tuple<double, double, double>> m_roc;
    vector<SROC> m_struct_roc;
public:
    CROCPlot(vector<double> &postiveScores, vector<double> &negativeScores);
    void saveROCtoTXT(const string& roc_data);
    double getAUC();
    void plotROCtoPNG(double auc, const string& roc_png_filename, const string& scorename);
};

class CTruth {
public:
    CTruth();
    virtual string getTruthFileName(){return "N/A";}
    virtual ~CTruth(); // virtual destruction
    virtual bool validate(int scan, string modified_pep);
    virtual string getTruth(int scan);
};

class CScanPepList: public CTruth {
    vector<tuple<int, string>> m_scanPepList;
    string m_truthfile;
public:
    CScanPepList(const string& scanpepListFile);
    bool validate(int scan, string modified_pep);
    string getTruth(int scan);
    string getTruthFileName(){return m_truthfile;}
};

class CTruthPepList: public CTruth{
    vector<string> m_pepList;

    string m_truthfile;
public:
    CTruthPepList(const string& pepListFile);
    bool validate(int scan, string modified_pep);
    string getTruthFileName(){return m_truthfile;}
};

shared_ptr<CTruth> CreateTruth(const string& filename, const string& method);

class HTMLReporter {
    string m_outfilename;
    ofstream m_fout;
public:
    HTMLReporter(string filename, const string& title = "PSM Validator");
    ~HTMLReporter();
    void addImage(const string& filename, const string& title, const string& description);
};

void load_tsv_file(const string &filename, vector<vector<string>> &data);

string getfragmodelfile(const string& inputfile, double m_minIntenFC, bool isTraining=false);

class CFlow {
    vector<CFlow *> following_steps;

public:
    string m_name;
    CFlow();
    virtual ~CFlow();
    void attach(CFlow *viewer);
    void subscribe(CFlow *previous_step);
    virtual void run();
    virtual void notify(const string &msg);
    virtual void update(string msg);
    virtual void update_with_key_value_pair(string key, string value){}
    string getname() const;
};

vector<Feature *>
createFeatureVector(const string& fragmodel, bool ghost, double minIntFC, const string& fragscoretype, const string& featurelistfile,
                    string binaryPath);

//void get_feature_from_spec(PSMInfo *psmInfo, CSpectrum *spec, vector<double> *OneSpecFeature, vector<Feature *> *features,
//                      int idx, string mzMLsourcefile, const bool fixMz, bool highMassAcc);

void
workOnOneBatch(PeptideProphetParser *ppp, const string& mzMLsourcefile, DataFile *df, vector<int> *index, vector<int> *batches,
               int i, vector<Feature *> *features, vector<vector<double>> *featuretable, shared_ptr<Progress> ps,
               bool fixMz, bool highMassAcc);

class ValidatePSM : public CFlow {
    string m_binaryPath;
    string m_pepxmlfile;
    string m_fragmodelfile;
    string m_fragmodelscoretype;
    double m_minIntFC;
    int m_threadNum;
    bool m_ghost;
    string m_basename;
    string m_validatorModelFile;
    vector<CFlow *> m_other_steps;
    string m_featurelistfile;
    bool fixMz;
    bool m_highMassAcc;
    bool m_useAlternativeProtein;
    string m_output_feature_file;
public:
    ValidatePSM(string pepxmlfile, const string& validatorModel, string fragmodelscoretype, double minInt, bool ghost,
                int threadNum, shared_ptr<CTruth> truth, int mtry, int ntree, string featurelistfile, int maxdepth,
                bool useAlternativeProt, string rangerbinary, string binaryPath);
    void getProcessed_i(DataFile *df, PeptideProphetParser &ppp,  vector<int> &index);
    void run();

    void exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                              const vector<vector<double>> &featuretable) const;

    void updatePsmTable(PeptideProphetParser &ppp, DataFile *df, const vector<int> &index, CTable &psmtable) const;
};

class ExtractFeaturesFromPepXML : public CFlow {
public:
    ExtractFeaturesFromPepXML(string pepxmlfile, const string& validatorModel, string fragmodelscoretype, double minInt,
                              bool ghost, int threadNum, string featurelistfile, bool useAlternativeProt,
                              string binaryPath, int hitrank);
    void run();
    ~ExtractFeaturesFromPepXML(){}
    void exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                              const vector<vector<double>> &featuretable) const;

    void updatePsmTable(ICPepXMLParser *pepxml, DataFile *df, const vector<int> &index, CTable &psmtable) const;

private:
    int m_hitrank;
    string m_binaryPath;
    string m_pepxmlfile;
    string m_fragmodelfile;
    string m_fragmodelscoretype;
    double m_minIntFC;
    int m_threadNum;
    bool m_ghost;
    string m_basename;
    string m_validatorModelFile;
//    vector<CFlow *> m_other_steps;
    string m_featurelistfile;
    bool fixMz;
    bool m_highMassAcc;
    bool m_useAlternativeProtein;
    string m_output_feature_file;
    char m_delimiter_psm_file;
};

class ExtractFeatures : public CFlow {
    string m_binaryPath;
    string m_pepxmlfile;
    string m_fragmodelfile;
    string m_fragmodelscoretype;
    double m_minIntFC;
    int m_threadNum;
    bool m_ghost;
    string m_basename;
    string m_validatorModelFile;
//    vector<CFlow *> m_other_steps;
    string m_featurelistfile;
    bool fixMz;
    bool m_highMassAcc;
    bool m_useAlternativeProtein;
    string m_output_feature_file;
    char m_delimiter_psm_file;
private:
    void getIndexPairs(DataFile *df, PeptideProphetParser &ppp, vector<pair<int, int>> &indexPair);
    void getProcessed_i(DataFile *df, PeptideProphetParser &ppp,  vector<int> &index);
    void updatePsmTable(PeptideProphetParser &ppp, DataFile *df, const vector<int> &index, CTable &psmtable) const;
    void exportTestingFeature(const vector<Feature *> &features, const string &feature_outfile,
                              const vector<vector<double>> &featuretable) const;
public:

    ExtractFeatures(string pepxmlfile, const string& validatorModel, string fragmodelscoretype, double minInt,
                    bool ghost, int threadNum, string featurelistfile, bool useAlternativeProt,
                    string binaryPath);

    void run();

    void printFeatureTable(const vector<Feature *> &features, const vector<vector<double>> &onefeaturetable) const;
};

void plot_FDR_curve(const string& outputfilename, const string& title, vector<tuple<double, double>> &fdr_counts);
void plot_ROC_with_score(CTable &psmtable, int column, const string& scoreName, string outprefix, bool onlytarget, HTMLReporter *reporter);

class resultAnalysis : public CFlow {
    string m_pepxmlfile;
    shared_ptr<CTruth> m_truth;
    string m_basename;
    bool m_useAlternativeProt;
public:
    resultAnalysis(string pepxmlfile, shared_ptr<CTruth> truth, string basename,bool useAlternativeProt=true);
    void run();
};

// TODO: there should be a table for all the random forest models
// each model is trained with different source files
// with different training samples
// with different parameters
// with different capabilities
// todo: finished the tasks above
class RangerFormat : public CFlow {
    string m_tsvfile;
    string m_positivefile;
    string m_negativefile;
    const int m_MIN_SAMPLE_NUM;// = 1000;
    const int m_MAX_SAMPLE_NUM;//=90000;
public:
    RangerFormat(const string &posFeature, const string &negFeature, const string &outPNsampleName, int trainingSampleSize);
    void run();
    int exportTrainingFeaturesToCSV(vector<vector<string>> &feature, const string &label, ofstream &fout, int start_col,
                                    int start_row, const int MAX_SAMPLE_NUM) const;
    void update_with_key_value_pair(string key, string value);
};


class AnnotationFormat{
    // an abstraction of every kind of annotation
    // iProphet result: pepXML
    // Search result: comet pepXML
    // Library spectra: sptxt
    // Note: feels that I have developed something like this...
};



class RFModelConfig{
    string m_tsvfile;
    bool m_isTraining;
    string m_rfModel;
    bool m_probPrediction;
    int m_mtry;
    int m_ntree;
    int m_maxDepth;
    string m_ranger_binary_path;
    int m_threadnum;
    string m_exportName;
public:
    const string &getMTsvfile() const;

    bool isMIsTraining() const;

    const string &getMRfModel() const;

    bool isMProbPrediction() const;

    int getMMtry() const;

    int getMNtree() const;

    int getMMaxDepth() const;

    const string &getMRangerBinaryPath() const;

    int getMThreadnum() const;

    const string &getMExportName() const;

public:
    RFModelConfig(string featuretsv, bool isTraining, string RFmodelfile, bool probPrediction, int mtry, int ntree,
    int maxdepth, string rangerbinary, string exportName);
    RFModelConfig(const string &configFileName);
    ~RFModelConfig();
    void write() ;
    void print();
};

class RangerWraper : public CFlow {
    string m_tsvfile;
    bool m_isTraining;
    string m_rfModel;
    bool m_probPrediction;
    string m_predictionprefix;
    int m_mtry;
    int m_ntree;
    int m_maxDepth;
    string m_ranger_binary_path;
    void train() const;
    void predict() const;
    int m_threadnum;
public:
    RangerWraper(RFModelConfig &rfconfig);
    void run();
    void update_with_key_value_pair(string key, string value);
    ~RangerWraper(){

    }
};

class FeatureWorkFlow : public CFlow {
    string m_inputfile;
    bool m_useGhostPeaks;
    double m_minInt; // bug fixed
    string m_fragmodelfile;
    string m_fragscoretype;
    string m_annotationMethod;  // Use title as sequence? or not
    string m_searchresult;
    bool m_debug;
    string m_feature_outputname;
    string m_featurelistfile;
    bool m_isLowMassAcc;
    string m_binaryPath;
public:
    // This is the workflow for generating model fragmentation
    FeatureWorkFlow(string inputfile, bool ghost, double minInt, string &fragmodelfile,
                    string fragscoretype, string annotationMethod, bool debug, const string& basename,
                    string featurelistfile, bool isLowMassAcc, string binaryPath);
    void run();
    void update_with_key_value_pair(string key, string value);
};

// todo: we are annotating with title
// in the future, we should use search result instead of title
class FragModel : public CFlow {
    string m_mgf_file;
    bool m_verbosity;
    double m_minIntenFC;
    bool m_overwritemodel;
    string m_modelFileName;
    bool m_exportFeatureToFile; // export feature and use binary train.
    bool m_testdecoy;
    bool m_use_ghost_peak;
    bool m_isLowMassAcc;
    string m_decoysearchfile;
    string m_binaryPath;
    void extractFeatures(bool verbosity, CPeakPairsImporter &positivePeakPairs, CPeakPairsImporter &negativePeakPairs);
public:
    FragModel(const string &mgf_input, bool verbosity, double minIntenFC, bool overwrite, string &fragPatternModel,
              bool outputfeature, bool useGhostPeak, bool isLowMassAcc, string binaryPath);
    ~FragModel();
    void run();
    void update_with_key_value_pair(string key, string value);
};

class FragPatternScoreFlow : public CFlow {
    string m_inputfile;
    string m_binaryPath;
    bool m_verbosity;
    double m_minIntenFC;
    bool m_overwritemodel;
    string m_modelfile;
    string m_searchresult; //  search result
    string m_scoretype;
    bool m_isLowMassAcc;
public:
    FragPatternScoreFlow(const string &inputfile, bool verbosity, double minIntenFC, bool overwrite, string outputmodelfile,
                         string scoretype, const string& searchresult, bool isLowMassAcc, string binaryPath);
    ~FragPatternScoreFlow();
    void run();
};

class CometSearch : public CFlow {
    string m_comet_binary;
    string m_database;
    string m_paramfile;
    string m_outputname;
    string m_specfile;
public:
    CometSearch(const string &cometbinary, const string &database, const string &paramfile, const string &outputname, const string &specfile);
    virtual ~CometSearch();
    virtual void run()  override;
};

class CometSearchTD : public CFlow {
    string m_comet_binary;
    string m_targetdb;
    string m_decoydb;
    string m_paramfile;
    string m_outputname;
    string m_specfile;
    CometSearch *m_pTarget;
    CometSearch *m_pDecoy;
public:
    CometSearchTD(string cometbinary, string targetdb, string decoydb, string paramfile, string specfile);
    ~CometSearchTD();
    void run();
    void update_with_key_value_pair(string key, string value);
};

// The entire workflow is here
// start from
// input:
//     MGF          -- input mass spectra file
//     database     -- input fasta file
//  for comet
//     cometinary   -- input comet binary path
//     paramfile    -- comet param
//
// workflow:
//      done:
//      Comet Search Target Decoy Separated
//              MGF --> target.pep.xml and decoy.pep.xml
//      Fragmentation Model train
//              MGF(Spectral Library) --> FragModel with title
//      Feature Extraction
//              MGF --> positive.feat(title); negative.feat(searchresult)
//      todo:
//      Final Model: validation of each PSM
//              MGF + Feature --> Score
class FlowAll : public CFlow {
    string m_inputspec;
    vector<CFlow *> m_pFlow;
public:
    FlowAll(const string& fragscoretype, bool useGhostPeak, bool outputfeature, bool overwrite, double minIntenFC, bool verbosity,
            string fragPatternModelFilename, const string& inputfile, string cometbinary, string targetdb, string decoydb,
            string paramfile, bool isTrainingRF, bool isRFProbOut, string writeRFModelTo, string validatorRFmodel,
            int mtry, int ntree, const string& featurelistfile, bool isLowMassAcc, int trainingSampleSize, int maxdepth,
            const string& rangerbinary, const string& binaryPath);
    void run();
    ~FlowAll();
};


#endif //MYTOOL_CFLOW_H
