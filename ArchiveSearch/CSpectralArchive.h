//
// Created by wulong on 8/6/18.
//

#ifndef MYTOOL_CINDEXMANAGER_H
#define MYTOOL_CINDEXMANAGER_H

#include <string>
#include <fstream>
#include <map>
using namespace std;
class CAnnotationDB;
class ICQuery;
class CPQParam;
class CDataBaseManager;
class CPValueCalculator;
class CMzFileReader;
class CMultiIndices;
class ICMzFactory;
class MzSpecInfo;
class ICMzFile;
class CAnnSpectra;
class PSMInfo;
class PeptideProphetParser;
class CMzSpec;
class CArxivSearchResult;
class SPsmAnnotation;
class SLinearRegressionModel;
class Chromatogram;


// thread safe.
struct CSocketServerSummary{
    int num_of_search_done;
    double time_used_for_current_search;
    double total_time_used_for_search;

    CSocketServerSummary();
    void print();
    void update(double time_used_in_sec);
};

struct SArchiveMatch {
    CAnnSpectra * m_annNodes;
    vector<SPsmAnnotation> m_gtInDB;
    int validation(const PSMInfo &searchpsminfo, float &dist, float &evalue, SPsmAnnotation * &matched_gtinfo_ptr);
    void retrieveAnnotation(CAnnotationDB & annoDB);
    int size() const;
    const SPsmAnnotation & getGt(int i) const;
    bool isGtSignificant(int i) const;
    bool isSsmSignificant(int i, double pvaluethreshold) const;
    double getPvalue(int i) const;
    double getDist(int i) const;
    double getDP(int i) const;

    long getQueryIdx() const;
    bool validationAllNeighbors(float &dist, float &evalue, SPsmAnnotation * &matched_gtinfo_ptr,  string gtPeptide, vector<string>& promisingAlternativeSolutionPeptides, int &ranking);

};

class CMzFileSearchResults {
    double m_pvalue_th;
    vector<SArchiveMatch> res;
    CAnnotationDB &m_AnnotationDB;
    int m_topHitsNum;
public:
    CMzFileSearchResults(CAnnotationDB &annodb, int topHitsNum);
    void retrieveAnnotation(vector<CAnnSpectra *> &annOfQueries);
    void exportTsv(string outfilename);
    void validation(string searchfile, CMzFileReader &mzfile);
private:
    void toTsv(const string &outfile, bool withheader, bool wrapLines);
    void createPvalueList(vector<double> &p_values, bool topHitOnly);
    void calcFDRThresholdPvalue(const string &outfile, bool topHitOnly);
    bool isSignificantID(const SArchiveMatch &arcMatch, const MzSpecInfo &scaninfo) const;
};

double intTol2double(int tol);

struct SAnnGTSummary
{
    long totalnum;
    long correctnum;
    map<long, int> m_queryindex2scan;
    map<long, string> m_queryindex2filename;
    shared_ptr<PeptideProphetParser> ppp;
    bool m_recallOfTrueNeighbor;
    int m_recallTNNtopK;
    double m_recallTNNminDP;
    ofstream m_fout;
    string outputfile;

    void setRecallTNN(bool recallTNN, int topK, double recallTNNminDP)
    {
        m_recallOfTrueNeighbor = recallTNN;
        m_recallTNNtopK = topK;
        m_recallTNNminDP = recallTNNminDP;
        if(m_recallOfTrueNeighbor){
            m_fout.open(outputfile.c_str(), ios::out|ios::app);
        }

    }

    void outputTofile(string &str){
        if(m_fout) m_fout << str ;
    }

    SAnnGTSummary();
    ~SAnnGTSummary(){}
    void print();
    void increase(bool iscorrect, long queryindex=-1);

    void setvalidationfile(string peptideprophetpepxmlfile);

    void setqueryindexMap(CMzFileReader &querySpectra);
};
string convertChromatogramToSVG(shared_ptr<Chromatogram> chr, string filename, string output_title="");
class CSpectralArchive {
private:
    SAnnGTSummary agtsummary;
    int m_dim;
    string m_mzXMLListFileName;
    string m_pepxmlFileName;
    string m_indexFileName;
    bool m_removeprecursor;
    bool m_useflankingbins;
    shared_ptr<CAnnotationDB> m_AnnotationDB; // SQL db
    const int PeakNumPerSpec;
    shared_ptr<CPValueCalculator> m_pc; // pvalue
    shared_ptr<CMzFileReader> m_csa;      // Precise dp calc. MZ file;
    shared_ptr<ICMzFile> m_pScorer;
	shared_ptr<CMultiIndices> m_indices; // six FAISS index
    shared_ptr<ICMzFactory> m_scorerFactory;
    bool m_verbose;

	int m_tol;
	int m_minPeakNum;
	bool m_usegpu;
public:
    CSpectralArchive(string mzXMLList, string pepxml, string indexfile, bool removeprecursor, bool useflankingbins,
                     int tol, int minPeakNum, bool myOwnIndex, CPQParam option, string indexstrings, bool usegpu,
                     bool rebuildsqldb, int seedpvalue, const int topPeakNum, bool createfilenameBlackList,
                     bool saveBackgroundScore, bool verbose);

    ~CSpectralArchive();

    void searchNeighborsWithin(double min_dp, long start = 0, long end = -1);

    string getPlatform();
    void update(string new_experimental_data, string new_search_result, string new_search_result_list,
                string new_experimental_datalist);

    void updateListOfSpectra(string spectraList);

    void searchQuery(long query_index, string &jsonstring, int topN, int calcEdge, int nprobe, vector<uint16_t> &query,
                     bool visualize);
    void addRemark(long query_index, string &remarks);
    void getRemark(long query_index, string &remarks);
    void searchMzFileInBatch(CMzFileReader &querySpectra, long first, long last, string validationfile,
                             int topHitsNum, int numOfpValues, bool recalltrueneighbor, int batchSize,
                             int bgspec_seed, int recallTNNtopK, double recallTNNminDP,
                             bool skipBackgroundScoreCalc, bool useFlankingBins);
    string searchPeptide(string peptide);
    string searchFileName(string filename, string startscan, string endscan);
    void setnProbe(int nprobe);
    long size();
    CMzSpec getMzSpec(long queryindex);


    void getRawSpec(long queryindex, vector<double> &mz, vector<double> &intensity);

    void getChromatograms(long queryindex, bool getChromatogram, bool updateChromSVG, double mzTolerance,
                          int charge, map<int, shared_ptr<Chromatogram>> &ChromatogramOfIsotopes,string &mzrtmap);


private:
	void appendFileName(const string &mzXMLfile);
	void addListOfRawData(const string &new_experimental_datalist);
	void addpepxmlfile(string pepxmllist, string new_gt_file);
    void addRawData(string mzXMLfile);
    void addSearchResult(string pepxmlfile);
    void addListOfSearchResults(string pepxmlfilelist);
	void calcPvalue(int query_index, int tol, string outputbasename, bool normalize) ;
    void createIndices(bool myOwnIndex, shared_ptr<CPQParam> option, string &indexstrings);
    void createMzFileObject();
    void do_pvalue_onfraction(vector<CAnnSpectra *> &vqrSimple, ICQuery &q, long query_index,
                              bool plot_figure, vector<int> *ptrDPs, SLinearRegressionModel *ptrLR);
    void do_neighbor_pvlaue_on_all(long query_index, vector<CAnnSpectra *> &annSpecList, CAnnSpectra *R1);
    void do_pairwise_distance(int calcEdge, long query_index, vector<CAnnSpectra *> &vecAnns, CArxivSearchResult &annResults);
    void do_pvalue_onall(long query_index, vector<CAnnSpectra*> &vqrSimple);
    void doAddSearchResult(string gtfile);
    void exportResponse(int query_index, vector<CAnnSpectra*> &annSpecList, string &jsonstring, vector<int> *ptrDPs, SLinearRegressionModel *ptrLRModel);
    void filterWithMinPeakNum(bool verbose, vector<long> &retIdx) ;
    void getAccurateTopNeighbor(ICQuery &query, vector<long> &topIdx, bool isLowMassAcc);
    void getkTrueNearestNeighbor(ICQuery &query, vector<vector<long>> &topKTrueNN, bool isLowMassAcc, int topK,
                                 int dp_UInt);
    int getDim();
    void getScorerFactoryPtr();
    void massdiff(ICQuery &q, CArxivSearchResult &annRes);
    void populateGroundTruthTable() ;
	void searchICQuery(int topN, ICQuery &query, int tolerance, bool verbose, CArxivSearchResult &archiveRes);
    void searchICQuery(ICQuery &query, int tolerance, bool verbose, CArxivSearchResult &archiveRes);
    void visualization_of_score_distribution(long query_index, vector<CAnnSpectra*> &annSpecList, bool tworound, int topN);
    void visualizeRawScores(const string &outputbasename, vector<int> &all_scores, bool noZeros) const;
    void visualizeCDF(string outputbasename, vector<int> &all_scores, bool noZero, bool logCDF) ;
    void updateIndex(bool verbosity);
    void syncIndicesWithSpecFileTable();
};


#endif //MYTOOL_CINDEXMANAGER_H
