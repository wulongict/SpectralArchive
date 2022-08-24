//
// Created by wulong on 4/12/19.
//

#ifndef MYTOOL_CANNOTATIONDB_H
#define MYTOOL_CANNOTATIONDB_H

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;
class CDataBaseManager;
class SPsmAnnotation;
class DataFile;
//class specfileinfo;
class ICGtInfoUpdate;
class CGtUpdater;
class PSMInfo;
class CArxivSearchResult;
class CDBEntry;


struct specfileinfo
{
    string filename;
    long start;
    long end;
    int fileid;
    specfileinfo();
    explicit specfileinfo(vector<string>& result);
    explicit specfileinfo(CDBEntry &dbentry);

    specfileinfo(const specfileinfo &other);

    specfileinfo & operator=(const specfileinfo &other);
    void init();
    void init(vector<string> & result);
    void init(CDBEntry &dbentry);
    specfileinfo(string f, long s, long e, int id);
    void display() const;
    int inValidFileID();
    bool isGood();
};

class CBatchSQL
{
    int m_batchsize;
    vector<string> manySQLs;
    shared_ptr<CDataBaseManager> m_dbmanager;
    bool m_verbose;
public:
    CBatchSQL(int batchsize, shared_ptr<CDataBaseManager> dbmanager, bool verbose);

    void append(string sql);
    void done();
    ~CBatchSQL();
};

class CBlackList{
    vector<long> m_blackList;
    map<long, bool> m_blackListmap;
    unordered_map<long,bool> m_blackListUnorderedMap;
    unordered_set<long> m_blackListUnorderedSet;

public:
    enum dataType {MAP, UNORDERED_MAP, UNORDERED_SET};
    CBlackList(){}
    int size(){return m_blackList.size();}
    void insert(vector<long> & blackList){
        m_blackList.insert(m_blackList.end(), blackList.begin(), blackList.end());
    }
    bool isBlack(dataType dt, long id);
    void push_back(long id){
        // todo: attention: only m_blackList is update???
        m_blackList.push_back(id);
    }

private:
    bool isBlackUnorderedMap(long id);
    bool isBlackUnorderedSet(long id);
    bool isBlackMap(long id);

};


class CAnnotationDB
{
    shared_ptr<CDataBaseManager> m_dbmanager;
    map<string, string> m_createTableSql;
    map<string, vector<string>> m_tablename2header;
    CBlackList m_cbl;
    bool m_createFileNameBlackList;

public:
    void deleteLastNFile(int n);
    string toJsonNode(CDBEntry &dbentry,int i);

    string toJsonNodes(CDBEntry &dbentry);

    CAnnotationDB(bool createFileNameBlackList);
    ~CAnnotationDB();

// return the list of raw files, used for refresh the mzxmlfile list.
    vector<string> getListOfSpecFiles();
 
    // the code to be updated!
    //    void set(shared_ptr<CDataBaseManager> dbmanager);
    void setDB(const string& filename);
    void addNeighborInfo(CArxivSearchResult &archiveRes);
    void createBlackListWithCE(bool verbose);
    void createBlackListWithFileName(const string &datafilename, bool verbose);
    // init database, create if not exist
    // creating empty tables in database.
    void createDatabase(bool rebuild, string dbfilename, bool verbose);

    // create table if not exist.
    bool createTable(const string& tableName, bool overwrite, bool verbosity);

    // create tables if not exist. call createTable
    void createTables(bool rebuild, bool verbose);
    shared_ptr<CDataBaseManager>  getDB();

    bool retrieveGtinfo(long residx, SPsmAnnotation &gtinfo);
    void fixChargeState(SPsmAnnotation &gtinfo);
    void fixAllChargeState();


    void clearTable(const string& table_name);

    shared_ptr<CDBEntry> searchPeptide(const string& peptide);
    void searchGTWithFileName(const string& filename, string startscan, string endscan,
                              shared_ptr<CDBEntry> &dbentry);

    void
    searchSigPeptide(double precursor_neutral_mass, int charge, const string& peptide, vector<vector<string>> &results,
                     const string& excludefile);

    void searchGtFileName(const string& basename, vector<string> &filenames);

    string getGroundTruthFileWithSameName(const string &specfile);

    void searchSpecFileName(const string& specfile, vector<string> &filenames);
    bool isSpecFileExist(const string& specfile);
    void getSpecFileRows(const string &raw_file_name, vector<vector<string>> &results);
    bool inSearchSpace(PSMInfo &x, const string& excludefilename);

private:
    void addFileNameToResultsOfGTRecords(vector<vector<string>> &results);
    void addFileNameToResultsOfGTRecords(CDBEntry &dbentry);
    void getGtinfoRow(long idx, CDBEntry &entry);
public:
    bool isTableExist(string table_name);
    void filterwithblacklist(bool verbose, vector<long> &retIdx);

    string getSpecFileNameByID(int fileid);
    int getTotalFileNum();

    long getTotalSpecNum();

    int getFileIDWithEndIdx(long index_total);

    specfileinfo getSpecFileInfo(int fileid);

    void populateGtfilesTable(const string& gtfile);

    shared_ptr<CBatchSQL> buildBatchSQL(int batchSize, bool verbose);

    shared_ptr<CGtUpdater> buildGtUpdater(const string& gtfile);

    void addNewGtFile(string gtfile);

    void insertTnnInfo(long queryindex, int topk, double mindp, int indexnum, long archivesize, int nprobe, string date, vector<long> &tnnidx, vector<double> &tnndp);

    bool getTNN(long queryindex, int topk, double mindp, int indexnum, long archivesize, int nprobe, vector<long> &tnnidx, vector<double> &tnndp);

    void add(specfileinfo &sfinfo);

    void updateSpecRemarks(long spec_id, string &remarks);

    void searchRemarks(long spec_id, vector<vector<string>> &remarks);

    int getLastFileIdInSpecFileTable();

    specfileinfo getLastSpecFile();

    // add a new raw data.
    void appendNewDataFile(DataFile &df);

    void appendNewDataFile(DataFile &df,  ICGtInfoUpdate *updater);

    void addColumnCEIfNotExist(bool verbosity);

    bool isColumnExists(string colName,string table, bool verbosity);

    void addColumnNeighborIfNotExist(bool verbosity);

    void addColumnAlterPepIfNotExist(bool verbosity);

    void addColumnRFScoreIfNotExist(bool verbosity);

    void appendGtFileList(const vector<string> &annotationFiles);

    specfileinfo addToGtTable(DataFile &df, ICGtInfoUpdate *updater);

    void appendRawDataList(const string &m_mzXMLListFileName);
};



class CGtUpdater{
    CAnnotationDB* m_annotationdb;
    string m_gtfile;
public:
    CGtUpdater(CAnnotationDB* annodb, string gtfile);
    ~CGtUpdater();
private:
    int getFileId(const string& eachfile);
    string createUpdateSqlOfRow(long i, ICGtInfoUpdate &gt);
    void update(const string& rawfilename, ICGtInfoUpdate &gt);

};

#endif //MYTOOL_CANNOTATIONDB_H
