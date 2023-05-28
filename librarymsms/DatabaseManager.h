//
// Created by wulong on 11/4/18.
//

#ifndef MYTOOL_DATABASEMANAGER_H
#define MYTOOL_DATABASEMANAGER_H

#include <sqlite3.h>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <map>
using namespace std;
class CTable;

class CDBEntry{
    map<string, vector<string>> m_data;
    bool m_colname_fixed;
public:
    int size(){
//        cout << "m_data size " << m_data.size() << endl;
//        for(auto & pair : m_data){
//            cout << pair.first << "\t" << pair.second.size() << endl;
//        }
        int num = 0;
        if(not m_data.empty()){
            num =m_data.begin()->second.size();
        }
        return num;
        //return m_data.begin()->second.size();
    }

    void shrinkto(int n)
    {
        if (n < 0)
        {
            cout << "invalid n " << n << endl;
        }
        else if (n > this->size())
        {
            cout << "n is larger than size " << n << " " << this->size() << endl;
        }
        else
        {
            for (auto &pair : m_data)
            {
                pair.second.resize(n);
            }
        }
    }

    bool empty(){
        return m_data.empty() or m_data.begin()->second.empty();
        }
    CDBEntry(){m_colname_fixed = false;}
    CDBEntry(const vector<string>& headers);
    void setValue(string header, int record_id, string value){
        if(m_data.count(header)==0 or record_id >= m_data[header].size()){
            throw runtime_error("Error: header and record out of range");
        }
        m_data[header][record_id] = value;
    }
    void extendTable(string header){
        if(m_data.find(header)!=m_data.end()){
            // found new header already exist. do nothing.
        }else{
            int n = size();
            m_data[header]=vector<string>(n,"");
            // fill in blank strings
        }
    }
    // add data
    void add(const string& header, const string& value);
    // add another dbentry
    void add(const CDBEntry& entry){
        for(auto & pair : entry.m_data){
            for(auto & value : pair.second){
                add(pair.first, value);
            }
        }
    }
    string get(string header, int record_id){
        string result;
        if(m_data.find(header)==m_data.end()){
//            cout << "Invalid header: " << header << endl;
            cout << "Invalid column --" << header << endl;
            cout << "valid header are " << endl;
            for(auto &x: m_data){
                cout << x.first << endl;
            }
        }
        else if(m_data[header].size()<=record_id) {
            cout << "not enough records: record id out of range " << record_id << " range: " << m_data[header].size() << endl;
        }else{
            result = m_data[header].at(record_id);
        }
        return result;
    }
    long getInt(string header, int record_id);
    double getFloat(string header, int record_id);

    void print();
};


// Here is how a callback function looks like
//typedef int (*sqlite3_callback)(
//        void*,    /* Data provided in the 4th argument of sqlite3_exec() */
//        int,      /* The number of columns in row */
//        char**,   /* An array of strings representing fields in the row */
//        char**    /* An array of strings representing column names */
//);

static int get_row_callback(void *data, int numCol, char **argv, char **azColName){
    // the columns of a gt table
    // ID = 0
    //FILEID = 0
    //MS2COUNTS = 0
    //PEPTIDE = UNKNOWN
    //SCORE = -1.0
    //SCAN = 11
    //CTERM = 0.0
    //NTERM = 0.0
    //MODIFICATION = UNMODIFIED
    //PRECURSOR = 462.151
    //CHARGE = 1
    //RT = 5.42103
    //PEPTIDEPROPHETPROB = 0.0
    //IPROPHETPROB = 0.0
    //RFSCORE = 0.0
    //ISDECOY = 0
    //SIGNIFICANCE = 0
    //PROTEIN =
    //CE =
    //ALTERPEPTIDE = NULL
    //NEIGHBOR =
    int i;
    vector<string> *pdata = static_cast<vector<string>*>(data);
    pdata->assign(numCol,"");
    for(i=0; i<numCol; i++){
//        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");

        (*pdata)[i] = argv[i]?argv[i]:"NULL";
    }
    //printf("\n");
    return 0;
}


static int get_row_callback_dbEntry(void *data, int numCol, char **argv, char **azColName){
    CDBEntry *pdata = static_cast<CDBEntry*>(data);

    for(int i=0; i<numCol; i++){
        pdata->add(azColName[i],argv[i]?argv[i]:"NULL");
    }
    return 0;
}

// vector<vector<string> > * data 
static int get_multiple_rows_callback(void *data, int numCol, char **argv, char **azColName){
    int i;
    vector<vector<string>> *pdata = (vector<vector<string>>*)data;

    pdata->emplace_back(vector<string>(numCol,""));
    int row_num = pdata->size()-1;
    for(i=0; i<numCol; i++){
        //printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
        (*pdata)[row_num][i] = argv[i]? argv[i]:"NULL";
    }
    //printf("\n");
    return 0;
}



class CDataBaseManager
{
    sqlite3 *m_db;

    void openDB(const string& databaseFileName);
public:
    CDataBaseManager(const string& databaseFileName);

    ~CDataBaseManager();

    void getRow(vector<string> &result, const string& idname, int idvalue, const string& tablename, bool verbose);
    void getRow_deprecated(vector<string> &result, const string& sql, bool verbose);
    void getRow(CDBEntry &result, const string& sql, bool verbose);
    void getMultipleRows(vector<vector<string>> &result, const string& sql, bool verbose);


    void execAsTransaction(const string& sql, bool verbose);

    bool tableExists(const string& tablename, bool verbosity);

    long getTotalRows(string tablename);

    void batchSQL(bool run, vector<string> &manySQLs, const string& sql, int batchsize, bool verbose);


};


void ConvertTable();


// This two functions are use to create the sql table from begining. However, we may want to update it gradually.
// SpecFiles table
// id, fileid, start, end
void batchInsertSpecfiles(CTable &specfiles, CDataBaseManager &cbm, int start, int end);

// groundtruth table
// id, fildid, ms2counts, peptide, score, scan, cterm, nterm, modification, precursor, charge, rt, peptide prophet prob, iprophet prob.
void batchInsertGroundTruth(CTable &groundtruth, CDataBaseManager &cbm, int start, int end);

#endif //MYTOOL_DATABASEMANAGER_H
