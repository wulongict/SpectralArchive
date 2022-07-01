//
// Created by wulong on 11/4/18.
//

#include "DatabaseManager.h"
#include "Util.h"


namespace GT_TABLE {
    enum HEADER {
        MS2IDX,
        FILLENAME,
        FILEID,
        MS2COUNTS,
        PEPTIDE,
        SCORE,
        SCAN,
        CTERM,
        NTERM,
        MODIFICATION,
        PRECURSOR,
        CHARGE,
        RT,
        PEPTIDEPROPHETPROB,
        IPROPHETPROB
    };

}

namespace SPECFILES_TABLE {
    enum HEADER {
        SPEC_FILENAME, SPEC_START, SPEC_END
    };
}


void createmzXMLTable() {
    string tablefile = "/data/wulong/data/honeybee/all_mzXML.txt";
    CTable mzxmlfiles(tablefile, '\t', false, 0);
    string databasename = "MassSpectra.db";
    CDataBaseManager cbm(databasename);

    /* Create SQL statement */
    string sql = "CREATE TABLE RAWDATAFILE("  \
         "ID                 INTEGER PRIMARY KEY AUTOINCREMENT ," \
         "FILENAME           TEXT    NOT NULL);";
//    cbm.createTableInDatabase(sql);
    cbm.execAsTransaction(sql, true);

    /* Insert items into table */
    string manySQLs;//"BEGIN; \n";
    for (int i = 0; i < mzxmlfiles.m_row; i++) {
        string filename = mzxmlfiles.getEntry(i, 0);
        sql = "INSERT INTO RAWDATAFILE (FILENAME) VALUES ('" + filename + "'); ";


        sql = "INSERT INTO RAWDATAFILE (FILENAME) \n" \
            "SELECT '" + filename + "' \n" \
            "WHERE NOT EXISTS(SELECT 1 FROM RAWDATAFILE WHERE FILENAME = '" + filename + "');";
        manySQLs += sql;
        manySQLs += "\n";

    }
    //manySQLs += "COMMIT;";
    cout << "running " << manySQLs << endl;
    cbm.execAsTransaction(manySQLs, true);


//    sql = "INSERT INTO COMPANY (ID,FILENAME) "  \
//         "VALUES (1, 'Paul', 32, 'California', 20000.00 ); " \
//         "INSERT INTO COMPANY (ID,NAME,AGE,ADDRESS,SALARY) "  \
//         "VALUES (2, 'Allen', 25, 'Texas', 15000.00 ); "     \

}


void batchInsertSpecfiles(CTable &specfiles, CDataBaseManager &cbm, int start, int end) {
    string sql;

    string manySQLs;
    for (int i = start; i < specfiles.m_row and i < end; i++) {
        // to use it
        string values = to_string(i) + ",";
        values += "'" + specfiles.getEntry(i, SPECFILES_TABLE::HEADER::SPEC_FILENAME) + "',";
        values += specfiles.getEntry(i, SPECFILES_TABLE::HEADER::SPEC_START) + ",";
        values += specfiles.getEntry(i, SPECFILES_TABLE::HEADER::SPEC_END);


        sql = "INSERT INTO SPECFILES (FILE_ID, FILENAME,START, END) "
              "VALUES (" + values + ");";
//        cout << sql << endl;
        manySQLs += sql;
        manySQLs += "\n";

    }

    cout << "running " << endl << "---" << start << "---" << end << "---" << endl;
    cbm.execAsTransaction(manySQLs, true);
}


void batchInsertGroundTruth(CTable &groundtruth, CDataBaseManager &cbm, int start, int end) {
    string sql;

    string manySQLs;
    for (int i = start; i < groundtruth.m_row and i < end; i++) {
        // to use it
        string values;
        values += groundtruth.getEntry(i, GT_TABLE::MS2IDX) + ",";
        // values += "'"+groundtruth.getEntry(i,GT_TABLE::FILLENAME)+"',";
        values += groundtruth.getEntry(i, GT_TABLE::FILEID) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::MS2COUNTS) + ",";
        values += "'" + groundtruth.getEntry(i, GT_TABLE::PEPTIDE) + "',";
        values += groundtruth.getEntry(i, GT_TABLE::SCORE) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::SCAN) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::CTERM) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::NTERM) + ",";
        values += "'" + groundtruth.getEntry(i, GT_TABLE::MODIFICATION) + "',";
        values += groundtruth.getEntry(i, GT_TABLE::PRECURSOR) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::CHARGE) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::RT) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::PEPTIDEPROPHETPROB) + ",";
        values += groundtruth.getEntry(i, GT_TABLE::IPROPHETPROB);

        sql = "INSERT INTO GROUNDTRUTH (ID,FILEID,MS2COUNTS,PEPTIDE,SCORE,SCAN, "
              "CTERM,NTERM, MODIFICATION,PRECURSOR,CHARGE,RT,PEPTIDEPROPHETPROB, IPROPHETPROB) "
              "VALUES (" + values + "); ";

        manySQLs += sql;
        manySQLs += "\n";

    }

    cout << "running " << endl << "---" << start << "---" << end << "---" << endl;
    cbm.execAsTransaction(manySQLs, true);
}


void creategroundtruthTable() {
    string gttablefile = "/data/wulong/data/honeybee/all_mzXML.txt_gt.tsv";

    string databasename = "MassSpectra.db";
    CDataBaseManager cbm(databasename);

    /* Create SQL statement */
    /*
     * IDX, FILENAME, FILEID, MS2COUNTS, PEPTIDE, SCORE, SCAN,
     * CTERM, NTERM,MODIFICATION, PRECURSOR, CHARGE, RT, PEPTIDEPROPHETPROB, IPROPHETPROB
     * */
    string sql = "CREATE TABLE GROUNDTRUTH("  \
         "ID                 INTEGER PRIMARY KEY AUTOINCREMENT ," \
         "FILENAME           TEXT    NOT NULL, " \
         "FILEID                INT,"
                 "MS2COUNTS             INT,"
                 "PEPTIDE               CHAR(100),"
                 "SCORE                 FLOAT,"
                 "SCAN                  INT,"
                 "CTERM                 FLOAT,"
                 "NTERM                 FLOAT,"
                 "MODIFICATION          TEXT,"
                 "PRECURSOR             FLOAT,"
                 "CHARGE                INT,"
                 "RT                    FLOAT,"
                 "PEPTIDEPROPHETPROB    FLOAT,"
                 "IPROPHETPROB          FLOAT);";
//    cbm.createTableInDatabase(sql);
    cbm.execAsTransaction(sql, true);


    CTable groundtruth(gttablefile, '\t', false, 0);
    /* Insert items into table */

    int batchsize = 10000;
    for (int i = 0; i < groundtruth.m_row; i += batchsize) {
        batchInsertGroundTruth(groundtruth, cbm, i, i + batchsize);
    }

    vector<string> result;
    cbm.getRow(result, "ID", 0, "GROUNDTRUTH", false);
    for (const auto& x: result) cout << x << endl;
}


void ConvertTable() {
//    createmzXMLTable();
    creategroundtruthTable();
}

// Open database for use
// create if not exist!
void CDataBaseManager::openDB(const string& databaseFileName) {
    if (m_db != nullptr) {
        cout << "Can not open another database before closing the previous one" << endl;
        exit(0);
    }
    char *zErrMsg = 0;
    int rc;
    cout << "[Info] Openting database " << databaseFileName << endl;
    rc = sqlite3_open(databaseFileName.c_str(), &m_db); // create database!

    if (rc) {
        fprintf(stderr, "Error: %s: %s\n", databaseFileName.c_str(), sqlite3_errmsg(m_db));
        exit(0);
    } else {
        fprintf(stderr, "Opened database successfully\n");
    }
}

void CDataBaseManager::execAsTransaction(const string& sql, bool verbose) {

//    cout << sql << endl;
    char *zErrMsg = nullptr;
    int rc;

    /* Execute SQL statement */
    rc = sqlite3_exec(m_db, "BEGIN;", nullptr, 0, &zErrMsg);
    rc = sqlite3_exec(m_db, sql.c_str(), nullptr, 0, &zErrMsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        if(verbose) cout << sql << endl;
        sqlite3_free(zErrMsg);
        sqlite3_exec(m_db, "ROLLBACK;", nullptr, 0, &zErrMsg);
    } else {
        sqlite3_exec(m_db, "COMMIT;", nullptr, 0, &zErrMsg);
//        if (verbose)
//            fprintf(stdout, "sql run successfully\n");
    }

}

CDataBaseManager::~CDataBaseManager() {
    if (m_db != nullptr) {
        cout << "Closing database" << endl;
        sqlite3_close(m_db); // close database
        m_db = nullptr;
    }
}


// new DBentry
void CDataBaseManager::getRow(CDBEntry &result, const string& sql, bool verbose) {
    char *zErrMsg = nullptr;
    int rc;

    rc = sqlite3_exec(m_db, sql.c_str(), get_row_callback_dbEntry, &result, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg); 
    } else {
        if (verbose)
            fprintf(stdout, "sql run successfully\n");
    }
}

void CDataBaseManager::getRow_deprecated(vector<string> &result, const string& sql, bool verbose) {
    char *zErrMsg = nullptr;
    int rc;

    rc = sqlite3_exec(m_db, sql.c_str(), get_row_callback, &result, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg); 
    } else {
        if (verbose)
            fprintf(stdout, "sql run successfully\n");
    }
}

void CDataBaseManager::getRow(vector<string> &result, const string& idname, int idvalue, const string& tablename, bool verbose) {
    char *zErrMsg = nullptr;
    int rc;

    string sql = "select * from " + tablename + " where " + idname + "=" + to_string(idvalue) + ";";

    /* Execute SQL statement */
    rc = sqlite3_exec(m_db, sql.c_str(), get_row_callback, &result, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    } else {
        if(verbose) fprintf(stdout, "sql run successfully\n");
    }
}

bool CDataBaseManager::tableExists(const string& tablename, bool verbosity) {

    char *zErrMsg = nullptr;
    int rc;

    CDBEntry dbentry({"type","name"});
    string sql = "select type, name from sqlite_master where type='table' and name='" + tablename + "';";

    /* Execute SQL statement */
    rc = sqlite3_exec(m_db, sql.c_str(), get_row_callback_dbEntry, &dbentry, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    } else {
        //fprintf(stdout, "sql run successfully\n");
    }
    if (dbentry.empty()) {
        if(verbosity)cout << "Table " << tablename << " does not exist" << endl;
        return false;
    } else {
        if(verbosity)cout << "Table " << tablename << " exists" << endl;
        return true;

    }

}


void CDataBaseManager::getMultipleRows(vector<vector<string>> &result, const string& sql, bool verbose) {
    char *zErrMsg = nullptr;
    int rc;

    /* Execute SQL statement */
    rc = sqlite3_exec(m_db, sql.c_str(), get_multiple_rows_callback, &result, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    } else {
        if (verbose)
            fprintf(stdout, "sql run successfully\n");
    }

}

// create database for use
// create if not exist
CDataBaseManager::CDataBaseManager(const string& databaseFileName) {
    m_db = nullptr;
    openDB(databaseFileName);
}

long CDataBaseManager::getTotalRows(string tablename) {
    // never used
    string sql="select count(*) as num_of_rows from " + tablename;

    CDBEntry dbentry;
    getRow(dbentry, sql, false);
    long row_num = dbentry.getInt("num_of_rows",0);

    return row_num;
}

void CDataBaseManager::batchSQL(bool run, vector<string> &manySQLs, const string& sql, int batchsize, bool verbose) {
    manySQLs.push_back(sql);
    if(batchsize<1) batchsize = 1;
    if(run or batchsize == manySQLs.size())
    {
        string sqls;
        for(const auto& eachsql: manySQLs) sqls += eachsql+"\n";
        if(!sqls.empty()) execAsTransaction(sqls, verbose);
        vector<string>().swap(manySQLs); // todo: better way to clear
    }
}

CDBEntry::CDBEntry(const vector<string>& headers):CDBEntry() {
    if(headers.empty()) {
        cout << "Error: invalid header. Empty header is not allowed in CDBEntry object" << endl;
    }    else {
        m_colname_fixed = true;
    }
    for(const auto& header: headers){
        m_data[header] = vector<string>();
        m_data[header].reserve(1024);
    }
}

void CDBEntry::add(const string& header, const string& value) {
    if(m_data.find(header)!=m_data.end())    {
        m_data[header].push_back(value);
    } else if ( not m_colname_fixed){
        m_data[header] = vector<string>();
        m_data[header].reserve(1024);
        m_data[header].push_back(value);
    }
    
    
}

void CDBEntry::print() {
    int size = m_data.begin()->second.size();
    for(int i = 0; i < size; i ++){
        cout << "record i = " << i << endl;
        for(auto x: m_data){
            cout << x.first << " = " << x.second[i] << endl;
        }
        cout << "---------- i = " << i << endl;
    }

}

long CDBEntry::getInt(string header, int record_id) {
    return strtol(get(std::move(header),record_id).c_str(),nullptr, 10);
}

double CDBEntry::getFloat(string header, int record_id) {
    return strtod(get(std::move(header),record_id).c_str(),nullptr);
}
