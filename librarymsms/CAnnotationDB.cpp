//
// Created by wulong on 4/12/19.
//


#include "DatabaseManager.h"
#include "ProteomicsDataTypes.h"
#include "ICGtInfoUpdate.h"
#include "CAnnotationDB.h"
#include "Util.h"
#include "CPSMAnnotation.h"
#include "XMLFileParser.h"
#include "CArchiveSearchReply.h"
#include <string>
#include <map>
#include <spdlog/spdlog.h>

// someone might attack the database by using sql injection
// so we need to check the input string
// if it contains invalid symbol, we return a empty string.
// those includes: ; % ' " and so on
string toSqlSafe(const string &inputString){
    if(inputString.find_first_of(";%\'\"")!=string::npos){
        // found invalid symbol
        return "";
    }
    return inputString;
}

CBatchSQL::CBatchSQL(int batchsize, shared_ptr<CDataBaseManager> dbmanager, bool verbose) {
    m_dbmanager = dbmanager;
    m_batchsize = batchsize;
    m_verbose = verbose;
}

void CBatchSQL::append(string sql) {
    m_dbmanager->batchSQL(false, manySQLs, sql, m_batchsize, m_verbose);
}

void CBatchSQL::done() {
    m_dbmanager->batchSQL(true, manySQLs, "", m_batchsize, m_verbose);
}

CBatchSQL::~CBatchSQL() {
    if (!manySQLs.empty()) done();
}



specfileinfo::specfileinfo(const specfileinfo &other) {
    cout << "copy constructor" << endl;
    fileid = other.fileid;
    filename = other.filename;
    start = other.start;
    end = other.end;
    display();
}

specfileinfo &specfileinfo::operator=(const specfileinfo &other) {
    cout << "assign constructor " << endl;
    this->fileid = other.fileid;
    this->filename = other.filename;
    this->start = other.start;
    this->end = other.end;
    display();
    return *this;
}

void specfileinfo::init() {
    filename = "";
    start = 0;
    end = 0;
    fileid = inValidFileID();
}

void specfileinfo::init(vector<string> &result) {
    CSqlSpecfileTableRow specrow(result);
    filename = specrow.toString(CSqlSpecfileTableRow::FILENAME);
    fileid = specrow.toInt(CSqlSpecfileTableRow::FILE_ID);
    start = specrow.toLong(CSqlSpecfileTableRow::START);
    end = specrow.toLong(CSqlSpecfileTableRow::END);
}

void specfileinfo::init(CDBEntry &dbentry) {
    filename = dbentry.get("FILENAME",0);
    fileid = dbentry.getInt("FILE_ID",0);
    start = dbentry.getInt("START",0);
    end = dbentry.getInt("END",0);

}

specfileinfo::specfileinfo(string f, long s, long e, int id) {
    filename = f;
    start = s;
    end = e;
    fileid = id;
}


void specfileinfo::display() const {
    cout <<"fileid\tstart\tend\tfilename" << endl;
    cout << fileid << "\t" << start << "\t" << end << "\t" << filename << endl; 
    // printf("fileid\tstart\tend\tfilename\n%d\t%d\t%d\t")
    // cout << "== spectra file information == start ==" << endl;
    // cout << "fileid: " << fileid << endl
    //      << "filename: " << filename << endl
    //      << "start: " << start << endl
    //      << "end: " << end << endl;
    // cout << "== spectra file information == end ==" << endl;
}

specfileinfo::specfileinfo(vector<string> &result) {
    if(result.empty()){
        init();
    }else{
        init(result);
    }
}

specfileinfo::specfileinfo() {init();}

specfileinfo::specfileinfo(CDBEntry &dbentry) {
//    cout << "using dbentry to init specfileinfo" << endl;
    if(dbentry.empty()){
        init();
    }else{
        init(dbentry);
    }

}

bool specfileinfo::isGood() {
    return inValidFileID() != fileid;
}

int specfileinfo::inValidFileID() {return -1;}


void CAnnotationDB::appendNewDataFile(DataFile &df) {
    RawDataFile rdf(df);
    appendNewDataFile(df, &rdf);
}


void CAnnotationDB::appendNewDataFile(DataFile &df, ICGtInfoUpdate *updater) {
    specfileinfo new_file = addToGtTable(df, updater);
    add(new_file);
}

specfileinfo CAnnotationDB::addToGtTable(DataFile &df, ICGtInfoUpdate *updater) {
    int ms2counts= 0;
    specfileinfo sfinfo = getLastSpecFile();
    shared_ptr<CBatchSQL> batchSQL = buildBatchSQL(1000, false);
    //Progress ps(df.getSpectrumNum(), "add raw data spectrum to groundtruth table ");
    for (int j = 0; j < df.getSpectrumNum(); j++) {
        //cout << "source file " << df.getSourceFileName() << endl;
        //ps.increase();
        CSpectrum *spec = df.getSpectrum(j);
        if (spec->getMSLevel() != 2) {
            continue;
        }
        SPsmAnnotation gtinfo;

        gtinfo.set(spec->getScanNum(), spec->m_parentMz, spec->m_precursor_charge, spec->m_rt_in_sec, sfinfo.fileid+1,
                   sfinfo.end + ms2counts, df.getSourceFileName(), ms2counts);
        updater->updateGtInfo(gtinfo); // do nothing
        batchSQL->append(gtinfo.createInsertSql());
        ms2counts++;
    }
    specfileinfo ret_info={df.getSourceFileName(), sfinfo.end, sfinfo.end + ms2counts, sfinfo.fileid + 1};
    return ret_info;

}

// todo: it is safer to get the value based on column names. As we have added several new columns now
CAnnotationDB::CAnnotationDB(bool createFileNameBlackList) {
    m_verbose = false;
    m_createFileNameBlackList = createFileNameBlackList;
//    m_verboseMsger = make_shared<CVerboseMessage>();
    m_dbmanager = nullptr;
    m_tablename2header["GROUNDTRUTH"] = {"ID","FILEID","MS2COUNTS","PEPTIDE","SCORE","SCAN","CTERM","NTERM","MODIFICATION","PRECURSOR","CHARGE",
                                         "RT","PEPTIDEPROPHETPROB","IPROPHETPROB","RFSCORE","ISDECOY","SIGNIFICANCE","PROTEIN","CE","ALTERPEPTIDE","NEIGHBOR","RESCUEDPEPTIDE"};

    m_createTableSql["GROUNDTRUTH"] = "CREATE TABLE GROUNDTRUTH("
                                      "ID                 INTEGER PRIMARY KEY AUTOINCREMENT ,"
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
                                      "IPROPHETPROB          FLOAT, "
                                      "RFSCORE          FLOAT, "
                                      "ISDECOY               INT, "
                                      "SIGNIFICANCE           INT, "
                                      "PROTEIN               varchar, "
                                      "CE                     varchar, "
                                      "ALTERPEPTIDE          varchar, "
                                      "RESCUEDPEPTIDE          varchar);"
                                      "CREATE INDEX GROUNDTRUTH_PEPTIDE ON GROUNDTRUTH (PEPTIDE);";


    m_tablename2header["TNN"] = {"ID", "QUERY_ID", "TOPK", "MINDP", "INDEXNUM", "ARCHIVESIZE", "NPROBE", "DATE", "TNN_IDX", "TNN_DP"};
    m_createTableSql["TNN"] = "CREATE TABLE TNN("
                                      "ID                 INTEGER PRIMARY KEY AUTOINCREMENT ,"
                                      "QUERY_ID                INT,"
                                      "TOPK             INT,"
                                      "MINDP               FLOAT,"
                                      "INDEXNUM                 INT,"
                                      "ARCHIVESIZE                  INT,"
                                      "NPROBE                 INT,"
                                      "DATE                 CHAR(100),"
                                      "TNN_IDX               varchar, "
                                      "TNN_DP                varchar);"
                                      "CREATE INDEX TNN_QUERY_ID ON TNN (QUERY_ID);";



    m_tablename2header["SPECFILES"] = {"FILE_ID","FILENAME","START","END"};
    m_createTableSql["SPECFILES"] = "CREATE TABLE SPECFILES("
                                    "FILE_ID INTEGER PRIMARY KEY AUTOINCREMENT, "
                                    "FILENAME TEXT NOT NULL, "
                                    "START INT NOT NULL, "
                                    "END INT NOT NULL);";

    m_tablename2header["GTFILES"] = {"FILE_ID","FILENAME"};
    m_createTableSql["GTFILES"] = "CREATE TABLE GTFILES("
                                  "FILE_ID INTEGER PRIMARY KEY AUTOINCREMENT, "
                                  "FILENAME TEXT NOT NULL); ";

    m_tablename2header["SPECREMARKS"] = {"ID","SPEC_ID","REMARKS"};
    m_createTableSql["SPECREMARKS"] = "CREATE TABLE SPECREMARKS("
                                      "ID INTEGER PRIMARY KEY AUTOINCREMENT, "
                                      "SPEC_ID INT, "
                                      "REMARKS TEXT NOT NULL);";
}

CAnnotationDB::~CAnnotationDB()
{
    
}

vector<string> CAnnotationDB::getListOfSpecFiles()
{
    // before release the database.
    vector<string> datafilelist;
    if (m_dbmanager)
    {

        int fileNum = this->getTotalFileNum();
        for (int i = 0; i < fileNum; i++)
        {
            string filename = this->getSpecFileNameByID(i);
            // cout << "File name in DB: " << i << " " << filename << endl;
            datafilelist.push_back(filename);
        }
        // try to refresh the mzxmlfile list
    }
    return datafilelist;
}

bool CAnnotationDB::createTable(const string& tableName, bool overwrite,bool verbosity) {
    bool renewed = false;
    // Invalid table name
    if (m_createTableSql.count(tableName) == 0) {
        cout << "Invalid SQL table name: " << tableName << endl;
        cout << "Supported table names are : " << endl;
        for (auto &item: m_createTableSql) {
            cout << item.first << "\t";
        }
        cout << endl;
        return renewed;
    }

    // Remove table first
    if (overwrite) {
        clearTable(tableName);
    }

    if (not m_dbmanager->tableExists(tableName, verbosity)) {
        cout << "[Info] Creating table " << tableName << endl;
        m_dbmanager->execAsTransaction(m_createTableSql[tableName], verbosity);
        renewed = true;
    }
    return renewed;
}

void CAnnotationDB::clearTable(const string& table_name) {
    string sql = "delete from " + toSqlSafe(table_name) + ";";
    m_dbmanager->execAsTransaction(sql, false);
    spdlog::get("A")->info("{} table deleted from db, table status = {}", table_name,
                           m_dbmanager->tableExists(table_name, false));

}

void CAnnotationDB::searchGtFileName(const string& basename, vector<string> &filenames) {

    string sql = "select * from GTFILES where FILENAME like  \"%" + toSqlSafe(basename) + "%\" limit 10;";
    CDBEntry dbentry;
    m_dbmanager->getRow(dbentry,sql, false);
    filenames.resize(dbentry.size());
    for(int i = 0; i < dbentry.size(); i ++){
        filenames[i] = dbentry.get("FILENAME",i);
    }
}


string CAnnotationDB::getGroundTruthFileWithSameName(const string &specfile) {
    string basename = CPath(specfile).m_filename_no_ext;
    vector<string> gtfiles;
    searchGtFileName(basename, gtfiles);

    string groundtruthfile;
    if (!gtfiles.empty()) {
        groundtruthfile = gtfiles[0];
        cout << "Got GT file: " << groundtruthfile << endl;
    }
    return groundtruthfile;
}

void CAnnotationDB::searchSpecFileName(const string& specfile, vector<string> &filenames) {
    // to be used in all the sql cases.
    if(string::npos!=specfile.find_first_of("\'\":,%$#*^!=\\;?><")){
        cout << "In valid string found in pattern " << specfile << endl;
        return;
    }
    string sql = "select * from SPECFILES where FILENAME like  \"%" + toSqlSafe(specfile) + "\" limit 10;";
    CDBEntry dbentry;
    m_dbmanager->getRow(dbentry,sql, false);
    filenames.resize(dbentry.size());
    for(int i = 0 ;i < dbentry.size(); i ++){
        filenames[i] = dbentry.get("FILENAME",i);
    }
}

// name is returned by the reference.
void CAnnotationDB::addFileNameToResultsOfGTRecords(CDBEntry &dbentry) {
    dbentry.extendTable("FILENAME");
    for(int i = 0; i < dbentry.size(); i ++){
        int fileid=dbentry.getInt("FILEID",i);
        string filename = getSpecFileNameByID(fileid);
        dbentry.setValue("FILENAME",i,File::CFile(filename).filename);
    }
}
void CAnnotationDB::addFileNameToResultsOfGTRecords(vector<vector<string>> &results) {
    map<string, string> fileid2Name;
    for (auto & res : results) {
        CSqlGtTableRow gt_row(res);
        string fileid = gt_row[CSqlGtTableRow::FILEID];
        if (fileid2Name.count(fileid) == 0) {
            string filenameX = getSpecFileNameByID(stringTo<int>(fileid));
            fileid2Name[fileid] = filenameX;
        }

        File::CFile fileObj(fileid2Name[fileid]);
        res.push_back(fileObj.filename);
    }
    cout << "File name added at the end of the results" << endl;
}


// single row
void CAnnotationDB::getGtinfoRow(long idx, CDBEntry &entry) {
    string sql = "select * from GROUNDTRUTH  where ID=" + to_string(idx);
    m_dbmanager->getRow(entry,sql,false);
}

string CAnnotationDB::getSpecFileNameByID(int fileid) {
    return getSpecFileInfo(fileid).filename;
}

specfileinfo CAnnotationDB::getSpecFileInfo(int fileid) {
    CDBEntry dbentry;
    string sql = "select * from SPECFILES where FILE_ID=" + to_string(fileid);
    m_dbmanager->getRow(dbentry,sql,false);

    if (dbentry.empty()) {
        cout << "[Warning]: Not found any file with ID: " << fileid << endl;
    }
    specfileinfo sfinfo(dbentry);

    return sfinfo;
}

// the GTFILES table is created but never used!
void CAnnotationDB::populateGtfilesTable(const string& gtfile) {
    if (File::isExist(gtfile,true)) {
        vector<string> annotationFiles = readlines(gtfile);
        appendGtFileList(annotationFiles);
    }
}

void CAnnotationDB::appendGtFileList(const vector<string> &annotationFiles) {
    int batchSize = 1;
    shared_ptr<CBatchSQL> batchSqls = buildBatchSQL(batchSize, true);

    for (const auto & annotationFile : annotationFiles) {
        vector<string> results;
        searchGtFileName(File::CFile(annotationFile).basename, results);
        if (results.empty()) {
            addNewGtFile(annotationFile);
//            string sql = "INSERT INTO GTFILES (FILENAME) VALUES('" + annotationFile + "');";
//            batchSqls->append(sql);
        }
    }
}

void CAnnotationDB::addNewGtFile(string gtfile) {
    string sql = "INSERT INTO GTFILES (FILENAME) VALUES('" + toSqlSafe(gtfile) + "');";
    cout << "running: " << sql << endl;
    m_dbmanager->execAsTransaction(sql, m_verbose);
}

void CAnnotationDB::add(specfileinfo &sfinfo) {
    if(m_verbose){
        sfinfo.display();
    }
    vector<string> filenames;
    searchSpecFileName(sfinfo.filename, filenames);
    if (filenames.empty()) {
        ostringstream oss;
        oss << sfinfo.fileid << ",'" << sfinfo.filename << "'," << sfinfo.start << "," << sfinfo.end << endl;
        string sql = "INSERT INTO SPECFILES (FILE_ID, FILENAME,START, END) "
                     "VALUES (" + oss.str() + ");";
        m_dbmanager->execAsTransaction(sql, m_verbose);
    } else {
        // If we skip this file, we should not add the corresponding peaks/annotations also.
        // TODO: conflicts caused by duplicate file in database: Sep 2020
        cout << "File already exist in sqlite3 database!  skip step" << endl;
        spdlog::get("A")->error("duplicate file detected! {}", sfinfo.filename);
    }

}

// single row
int CAnnotationDB::getLastFileIdInSpecFileTable() {
    int last_fileid = -1;
    if (getTotalFileNum() > 0) {
        string sql = "select max(FILE_ID) as max_file_id from SPECFILES;";
        CDBEntry dbentry;
        m_dbmanager->getRow(dbentry,sql,m_verbose);
        if(not dbentry.empty()){
            last_fileid = dbentry.getInt("max_file_id",0);
        }

    }
    return last_fileid;
}

specfileinfo CAnnotationDB::getLastSpecFile() {
    int last_fileid = getLastFileIdInSpecFileTable();
    if (last_fileid == -1) return specfileinfo();
    else return getSpecFileInfo(last_fileid);
}

// search for rows with the given filename and scan range. 
void
CAnnotationDB::searchGTWithFileName_new(const string& filename, string startscan, string endscan,
                                    shared_ptr<CDBEntry> &dbentry) {
    SimpleTimer st("searching for annotation");
    dbentry = make_shared<CDBEntry>();


    string filename_sql = "pragma case_sensitive_like=1;select * from SPECFILES where FILENAME like \"%" + toSqlSafe(filename) + "%\" limit 10";
    CDBEntry dbentry_filename;
    m_dbmanager->getRow(dbentry_filename,filename_sql,false);
    if(dbentry_filename.empty()){
        cout << "no file found with name: " << filename << endl;
        return;
    }
    // dbentry_filename.print();
    int totalnum = 0;
    for(int i = 0; i < dbentry_filename.size(); i ++){
        string fileid = dbentry_filename.get("FILE_ID",i);
        string start_idx = dbentry_filename.get("START",i);
        string end_idx = dbentry_filename.get("END",i);
        string current_filename = dbentry_filename.get("FILENAME",i);

        // cout << "dbentry_filename: " << current_filename << " " << start_idx << " " << end_idx << endl;

        string gq_sql = "pragma case_sensitive_like=1;select * from GROUNDTRUTH inner join SPECFILES on SPECFILES.FILE_ID=GROUNDTRUTH.FILEID where ID >= "+ start_idx +" and ID <=  "+ end_idx+ " and FILEID = " + fileid + " and SCAN >= " + startscan + " and SCAN <= "
                    + endscan + " limit 100 ;";
        CDBEntry dbentry_gq;
        m_dbmanager->getRow(dbentry_gq,gq_sql,false);
        if(dbentry_gq.empty()){
            cout << "no groundtruth found with name: " << fileid << " " << current_filename << endl;
        
        }
        totalnum += dbentry_gq.size();
        cout << dbentry_gq.size() << "/" << totalnum << " groundtruth found with name: " << filename  << " in new approach"<< endl;

        dbentry->add(dbentry_gq);
        if(dbentry->size() > 1000){
            break;
        }
        
    }
    // cout << "total number of ground truth found " << totalnum  << " new_dbentry size " << new_dbentry->size()<< endl;
    dbentry->shrinkto(100);

    // clean the file name remove path.
    for(int i = 0; i < dbentry->size(); i ++){
        File::CFile fullname(dbentry->get("FILENAME",i));
        dbentry->setValue("FILENAME",i,fullname.filename);
    }
    // cout << "total number of ground truth found " << totalnum  << " new_dbentry size " << new_dbentry->size()<< endl;

    // checking each entry.
    bool compare_two_version = false;
    if (compare_two_version)
    {
        shared_ptr<CDBEntry> old_dbentry = make_shared<CDBEntry>();
        searchGTWithFileName(filename, startscan, endscan, old_dbentry);
        for (int i = 0; i < dbentry->size(); i++)
        {
            // compare with the old one.
            string id = dbentry->get("ID", i);
            string scan = dbentry->get("SCAN", i);
            string fileid = dbentry->get("FILEID", i);
            string filename = dbentry->get("FILENAME", i);
            string charge = dbentry->get("CHARGE", i);
            string peptide = dbentry->get("PEPTIDE", i);
            string protein = dbentry->get("PROTEIN", i);
            string score = dbentry->get("SCORE", i);
            string neighbor = dbentry->get("NEIGHBOR", i);
            string isdecoy = dbentry->get("ISDECOY", i);

            // the old one
            string old_id = old_dbentry->get("ID", i);
            string old_scan = old_dbentry->get("SCAN", i);
            string old_fileid = old_dbentry->get("FILEID", i);
            string old_filename = old_dbentry->get("FILENAME", i);
            string old_charge = old_dbentry->get("CHARGE", i);
            string old_peptide = old_dbentry->get("PEPTIDE", i);
            string old_protein = old_dbentry->get("PROTEIN", i);
            string old_score = old_dbentry->get("SCORE", i);
            string old_neighbor = old_dbentry->get("NEIGHBOR", i);
            string old_isdecoy = old_dbentry->get("ISDECOY", i);

            // compare
            if (id != old_id)
            {
                cout << "id not equal " << id << " " << old_id << endl;
            }
            if (scan != old_scan)
            {
                cout << "scan not equal " << scan << " " << old_scan << endl;
            }
            if (fileid != old_fileid)
            {
                cout << "fileid not equal " << fileid << " " << old_fileid << endl;
            }
            if (filename != old_filename)
            {
                cout << "filename not equal " << filename << " " << old_filename << endl;
            }
            if (charge != old_charge)
            {
                cout << "charge not equal " << charge << " " << old_charge << endl;
            }
            if (peptide != old_peptide)
            {
                cout << "peptide not equal " << peptide << " " << old_peptide << endl;
            }
            if (protein != old_protein)
            {
                cout << "protein not equal " << protein << " " << old_protein << endl;
            }
            if (score != old_score)
            {
                cout << "score not equal " << score << " " << old_score << endl;
            }
            if (neighbor != old_neighbor)
            {
                cout << "neighbor not equal " << neighbor << " " << old_neighbor << endl;
            }
            if (isdecoy != old_isdecoy)
            {
                cout << "isdecoy not equal " << isdecoy << " " << old_isdecoy << endl;
            }
            cout << "-----------------------" << endl;
            cout << "entry i " << i << " are the same in the two dbentry object." << endl;
        }
    }
}




// search for rows with the given filename and scan range. 
void
CAnnotationDB::searchGTWithFileName(const string& filename, string startscan, string endscan,
                                    shared_ptr<CDBEntry> &dbentry) {
    SimpleTimer st("searching for annotation");
    // string sql = "pragma case_sensitive_like=1;select * from GROUNDTRUTH where SCAN >= " + startscan + " and SCAN <= "
    //              + endscan + " and FILEID in (select FILE_ID from SPECFILES where FILENAME like \"%"
    //              + filename + "%\") limit 1000 ;";


    string sql = "pragma case_sensitive_like=1;select * from GROUNDTRUTH inner join SPECFILES on SPECFILES.FILE_ID=GROUNDTRUTH.FILEID where SCAN >= " 
    + startscan + " and SCAN <= " + endscan + " and FILENAME LIKE \"%"
    + toSqlSafe(filename) + "%\" limit 100";

    cout << "searching database with sql : " << sql << endl;

    dbentry = make_shared<CDBEntry>();

    m_dbmanager->getRow(*dbentry, sql, false);
    // clean the file name remove path.
    for(int i = 0; i < dbentry->size(); i ++){
        File::CFile fullname(dbentry->get("FILENAME",i));
        dbentry->setValue("FILENAME",i,fullname.filename);
    }
}
void CAnnotationDB::fixChargeState(SPsmAnnotation &gtinfo){
    if(gtinfo.fixChargeStates())
    {
        string sql = gtinfo.createChargeUpdateSql();
        this->m_dbmanager->execAsTransaction(sql,false);
    }
}


bool CAnnotationDB::retrieveGtinfo(long residx, SPsmAnnotation &gtinfo) {
    bool ret = false;
    vector<string> headers = {"ID","NEIGHBOR"};
    headers = m_tablename2header["GROUNDTRUTH"];
    CDBEntry dbEntry(headers);
    getGtinfoRow(residx, dbEntry);
    if (not dbEntry.empty()) {
        ret = true;
        try {
            gtinfo.initWithRow(dbEntry);
        } catch (...) {
            cout << "error inside init" << endl;
            throw runtime_error("init");
        }
        string full_path_name = getSpecFileNameByID(gtinfo.fileid);
        gtinfo.mzxml_filename = CPath(full_path_name).m_filename;
    }
    return ret;
}


void CAnnotationDB::searchSigPeptide(double precursor_neutral_mass, int charge, const string& peptide,
                                     vector<vector<string>> &results,
                                     const string& excludefile) {
    File::CFile fileObj(std::move(excludefile));
    string sql = "select * from GROUNDTRUTH where PEPTIDE=\"" + toSqlSafe(peptide) + "\" and SIGNIFICANCE=1 and CHARGE=" +
                 to_string(charge) +
                 " and fileid not in (select file_id from SPECFILES where filename like '%" + toSqlSafe(fileObj.basename) + "%')" +
                 " limit 1;";

    m_dbmanager->getMultipleRows(results, sql, false);
}

void CAnnotationDB::updateSpecRemarks(long spec_id, string &remarks) {
    ostringstream oss;
    oss << spec_id << "," << "'" << toSqlSafe(remarks) << "'";
    string sql = "INSERT INTO SPECREMARKS( SPEC_ID, REMARKS) VALUES(" + oss.str() + ");";
    m_dbmanager->execAsTransaction(sql, true);
}

// now get all the remarks for each spectra
// ID
// SPEC_ID
// REMARKS
void CAnnotationDB::searchRemarks(long spec_id, vector<vector<string>> &remarks) {
    string sql = "select * from SPECREMARKS WHERE SPEC_ID=" + to_string(spec_id) + ";";
    m_dbmanager->getMultipleRows(remarks, sql, true);
}


// search for proper HCD energy
// the data file ~/speclib/human_hcd_selected.sptxt contains different energy.
void CAnnotationDB::createBlackListWithCE(bool verbose) {
    string sql = R"(select ID,PEPTIDE,CE,CHARGE from GROUNDTRUTH where CE is not NULL and CE!="";)";

    CDBEntry dbentry;
    m_dbmanager->getRow(dbentry,sql,false);

    float currentBest = 0;
    string prevPeptide;
    int prevcharge = 0;
    long black_ID_candidate = -1;

    auto isBetterCE =[](float newCE, float oldCE)->bool{
        const float best_CE = 35; // eV
        return fabs(newCE-best_CE) < fabs(oldCE-best_CE);
    };

    for(int k = 0; k < dbentry.size(); k ++){
        long id = dbentry.getInt("ID", k);
        string pep = dbentry.get("PEPTIDE",k);
        string cestr = dbentry.get("CE",k);
        int charge = dbentry.getInt("CHARGE",k);

        int ce_len = cestr.length();
        float ce = 0.0f;
        if (ce_len > 2) {
            cestr = cestr.substr(0, ce_len - 2);
            ce = stringTo<float>(cestr);
        }

        if (prevcharge != charge or prevPeptide != pep) {
            currentBest = ce;
            prevcharge = charge;
            prevPeptide = pep;
            black_ID_candidate = id;
            continue;
        } else if (isBetterCE(ce, currentBest)) {
            // try blocklist. not black list.
            // same peptide, charge, and collision energy(CE) closer to 35eV, put the previous tto blacklist.
            m_cbl.push_back(black_ID_candidate);
            // update , peptide, charge are the same.
            black_ID_candidate = id;
            currentBest = ce;
        } else {
            // same peptide, charge, and previous CE is closer to 35eV, put the current one to blacklist.
            m_cbl.push_back(id);
        }
    }

    if(verbose)cout << "black ID list created of size: " << m_cbl.size() << endl;

}

// Create a black list of unique IDs to avoid reporting candiates/neighbors from the same file
void CAnnotationDB::createBlackListWithFileName(const string &datafilename, bool verbose) {
    if (not m_createFileNameBlackList) {
        if(verbose)cout << "Attention: file name black list will not be created" << endl;
        return;
    }
    if(verbose) cout << typeid(this).name() << " Creating black list with file name: " << datafilename << endl;
    File::CFile fileObj(datafilename);

    string sql = "select * from SPECFILES where filename like \"%" + toSqlSafe(fileObj.filename) + "%\"";
    vector<vector<string>> results;
    m_dbmanager->getMultipleRows(results, sql, false);
    if (not results.empty()) {
        cout << typeid(this).name() << " Found " << datafilename << endl;
        for (auto res: results) {
            CSqlSpecfileTableRow spec_row(res);
            int start_idx = stringTo<int>(spec_row.getResults()[CSqlSpecfileTableRow::START]);
            int end_idx = stringTo<int>(spec_row.getResults()[CSqlSpecfileTableRow::END]);
            string filename = spec_row.getResults()[CSqlSpecfileTableRow::FILENAME];
            cout << filename << "\t" << start_idx << "\t" << end_idx << endl;
            vector<long> file_name_black_list(end_idx - start_idx, 0);
            iota(file_name_black_list.begin(), file_name_black_list.end(), start_idx);
            cout << "m_blacklist will be extended " << *file_name_black_list.begin() << " to "
                 << *file_name_black_list.rbegin() << endl;
            cout << m_cbl.size() << endl;
            m_cbl.insert(file_name_black_list);

            cout << m_cbl.size() << endl;


        }
    }
}

// Not used: to be deleted!
shared_ptr<CDataBaseManager> CAnnotationDB::getDB() {
    return m_dbmanager;
}
bool CAnnotationDB::isColumnExists(string colName, string table, bool verbosity) {
    string sql = R"(select sql from sqlite_master where type='table' and name like "%)" + toSqlSafe(table) + R"(%" and sql like '%, )" + toSqlSafe(colName) + R"( %';)";
    CDBEntry dbEntry;
    m_dbmanager->getRow(dbEntry, sql, verbosity);
    return not dbEntry.empty();
}

void CAnnotationDB::addColumnCEIfNotExist(bool verbosity) {
    if(isColumnExists("CE","GROUNDTRUTH",verbosity)){
        if (verbosity){
            cout << "column CE already exist in table GROUNDTRUTH; will NOT add again" << endl;
        }
    } else {
        string sql = R"(ALTER TABLE GROUNDTRUTH ADD CE VARCHAR)";
        cout << "to add one more column to our table " << sql << endl;
        //ALTER TABLE table_name
        //  ADD new_column_name column_definition;
        m_dbmanager->execAsTransaction(sql, verbosity);
    }

}

void CAnnotationDB::addColumnNeighborIfNotExist(bool verbosity) {
    string sql = R"(select sql from sqlite_master where type='table' and name like "%GROUNDTRUTH%" and sql like '%, NEIGHBOR %';)";


    CDBEntry dbEntry;
    m_dbmanager->getRow(dbEntry, sql, verbosity);
    if (dbEntry.empty()) {
        // cout << "to add one more column NEIGHBOR to our table " << sql << endl;
        //ALTER TABLE table_name
        //  ADD new_column_name column_definition;
        sql = R"(ALTER TABLE GROUNDTRUTH ADD NEIGHBOR VARCHAR)";
        m_dbmanager->execAsTransaction(sql, verbosity);
    } else if (verbosity){
        cout << "column NEIGHBOR already exist in table GROUNDTRUTH; will NOT add again" << endl;
    }

}


void CAnnotationDB::addColumnRFScoreIfNotExist(bool verbosity) {
    string sql = R"(select sql from sqlite_master where type='table' and name like "%GROUNDTRUTH%" and sql like '%, RFSCORE %';)";

    CDBEntry dbEntry;
    m_dbmanager->getRow(dbEntry, sql, verbosity);
    if (dbEntry.empty()) {
        cout << "[Info] to add column RFSCORE to the table GROUNDTRUTH " << sql << endl;
        sql = R"(ALTER TABLE GROUNDTRUTH ADD RFSCORE FLOAT)";
        m_dbmanager->execAsTransaction(sql, verbosity);
    } else if (verbosity){
        cout << "column RFSCORE already exist in table GROUNDTRUTH; will NOT add again" << endl;
    }

}

void CAnnotationDB::addColumnRescuedPepIfNotExist(bool verbosity) {
    string sql = R"(select sql from sqlite_master where type='table' and name like "%GROUNDTRUTH%" and sql like '%RESCUEDPEPTIDE %';)";

    CDBEntry dbEntry;
    m_dbmanager->getRow(dbEntry, sql, verbosity);
    if (dbEntry.empty()) {
        sql = R"(ALTER TABLE GROUNDTRUTH ADD RESCUEDPEPTIDE VARCHAR)";
        m_dbmanager->execAsTransaction(sql, verbosity);
    } else if (verbosity){
        cout << "column RESCUEDPEPTIDE already exist in table GROUNDTRUTH; will NOT add again" << endl;
    }

}

void CAnnotationDB::addColumnAlterPepIfNotExist(bool verbosity) {
    string sql = R"(select sql from sqlite_master where type='table' and name like "%GROUNDTRUTH%" and sql like '%, ALTERPEPTIDE %';)";

    CDBEntry dbEntry;
    m_dbmanager->getRow(dbEntry, sql, verbosity);
    if (dbEntry.empty()) {
        sql = R"(ALTER TABLE GROUNDTRUTH ADD ALTERPEPTIDE VARCHAR)";
        m_dbmanager->execAsTransaction(sql, verbosity);
    } else if (verbosity){
        cout << "column ALTERPEPTIDE already exist in table GROUNDTRUTH; will NOT add again" << endl;
    }

}

long CAnnotationDB::getTotalSpecNum() {
    return getLastSpecFile().end;

}

bool CAnnotationDB::isTableExist(string table_name) {
    return m_dbmanager->tableExists(std::move(table_name), false);
}

int CAnnotationDB::getTotalFileNum() {
    vector<string> results;
    string sql = "select count(*) from SPECFILES";
    m_dbmanager->getRow_deprecated(results, sql, false);
    int total_row_num = stringTo<long>(results[0]);
    return total_row_num;
}

int CAnnotationDB::getFileIDWithEndIdx(long index_total) {
    vector<string> results;
    m_dbmanager->getRow(results, "END", index_total, "SPECFILES", false);
    int start_id = 0;

    if (not results.empty() ) {
        start_id = strtol(results[0].c_str(), nullptr, 10) + 1;
    } else if (index_total >0){
        cout << "Corruption detected inside index/spec_table pair! Rebuild the index to fix the problem. " << endl;
        exit(0);
    }
    return start_id;
}

// open database file for use 
// create database if not exist
void CAnnotationDB::setDB(const string& filename) {
    m_dbmanager = make_shared<CDataBaseManager>(filename);

}

// search for raw file with name like the `raw_file_name`
void CAnnotationDB::getSpecFileRows(const string &raw_file_name, vector<vector<string>> &results) {
    string name_only = File::CFile(raw_file_name).basename;
    // todo: add protection here.
    string sql = "select * from SPECFILES where FILENAME like \"%" + toSqlSafe(name_only) + "%\";";
    
    m_dbmanager->getMultipleRows(results, sql, false);
}

shared_ptr<CBatchSQL> CAnnotationDB::buildBatchSQL(int batchSize, bool verbose) {
    return make_shared<CBatchSQL>(batchSize, m_dbmanager, verbose);
}

bool CAnnotationDB::isSpecFileExist(const string& specfile) {
    vector<string> filenames;
    searchSpecFileName(specfile, filenames);
    return not filenames.empty();
}

shared_ptr<ICGtInfoUpdate> factoryMethodGtInfoUpdate(string ground_truth_file, DataFile &df) {
    if (ground_truth_file.find(".pep.xml") != string::npos and
        ground_truth_file.find("interact-") != string::npos) {
        return make_shared<PeptideProphetParser>(ground_truth_file);
    } else if (ground_truth_file.find(".pep.xml") != string::npos or ground_truth_file.find(".pepXML") != string::npos) {
        return make_shared<CometPepXMLParser>(ground_truth_file);
    } else if (ground_truth_file.find(".sptxt") != string::npos or ground_truth_file.find(".mgf") != string::npos) {
        return make_shared<AnnotationDataFile>(df);
    } else {
        if (not ground_truth_file.empty()) {
            cout << "[Warning] the following file is not supported! " << ground_truth_file << endl;
        }
        return make_shared<RawDataFile>(df);
    }
}

void CAnnotationDB::filterwithblacklist(bool verbose, vector<long> &retIdx) {
    auto it = remove_if(retIdx.begin(), retIdx.end(),[&](const long &x){return m_cbl.isBlack(CBlackList::UNORDERED_SET, x);});
    if(verbose) cout << "Num of spectra removed as black list " << std::distance(it, retIdx.end()) << endl;
    retIdx.erase(it,retIdx.end());
}

// Create tables if not exist
// Table: SPECFILES SPECREMARKS GROUNDTRUTH
void CAnnotationDB::createTables(bool rebuild, bool verbose) {
//    cout << "Create Tables " << endl;
    // create all the tables
    for(auto &eachTable : m_createTableSql){
        bool tableRenewed = createTable(eachTable.first,rebuild, verbose);
        if(tableRenewed){
           if(verbose) cout << "[Info] table " << eachTable.first << " created!" << endl;
        } else {
            if(verbose) cout << "[Info] table " << eachTable.first << " already exists!" << endl;
        }
    }
//    if(verbose)
//    bool verbosity = true;
    addColumnCEIfNotExist(verbose);
    addColumnAlterPepIfNotExist(verbose);
    addColumnRFScoreIfNotExist(verbose);
    addColumnNeighborIfNotExist(verbose);
    addColumnRescuedPepIfNotExist(verbose);
    cout << "Tables created" << endl;


}

void CAnnotationDB::appendRawDataList(const string &m_mzXMLListFileName) {// IT IS THE GROUNDTRUTH TABLE JUST CREATED!
// TODO: detect duplicates.
    if (m_mzXMLListFileName.empty()){
        cerr << "the filename for list of mzxml is empty" << endl;
        throw runtime_error("empty mzxml list file name");
    }
    //cout << "add raw data to sql db : " << m_mzXMLListFileName << endl;
    CTable mzXMLTable(m_mzXMLListFileName, '\t', false, 0);
    //Progress ps(mzXMLTable.m_row, "Add raw data files");
    int skipped_file_num = 0;
    for (int i = 0; i < mzXMLTable.m_row; i++) {
        string specfile = mzXMLTable.getEntry(i, 0);
        //ps.increase();
        if(isSpecFileExist(specfile)){
            skipped_file_num ++;
//            spdlog::get("A")->info("Processing file {} / {} : {} already exist, skip it",
//                                   i + 1, mzXMLTable.m_row, specfile);
            continue;
        }
        spdlog::get("A")->info("\nProcessing file {} / {} : {} ", i + 1, mzXMLTable.m_row, specfile);
        DataFile expDF(specfile, 0, -1);
        string ground_truth_file = getGroundTruthFileWithSameName(specfile);
        shared_ptr<ICGtInfoUpdate> p = factoryMethodGtInfoUpdate(ground_truth_file, expDF);
        appendNewDataFile(expDF, p.get());
    }
    cout << skipped_file_num << " raw data files already exist and are skipped. " << endl;
}

shared_ptr<CGtUpdater> CAnnotationDB::buildGtUpdater(const string& gtfile) {
    return make_shared<CGtUpdater>(this, gtfile);
}

bool CAnnotationDB::inSearchSpace(PSMInfo &x, const string& excludefilename) {
    // this function looks for the exact peptide in database.
    // x.searchhits[0].m_peptide does not contain modifications
    vector<vector<string>> results;
    searchSigPeptide(x.precursorNeutMass, x.charge, x.searchhits[0]->m_peptide, results, excludefilename);
    return not results.empty();
}

void CAnnotationDB::addNeighborInfo(CArxivSearchResult &archiveRes) {
    // requirement: the NEIGHBOR column should exist.
    shared_ptr<CBatchSQL> batch_sql = make_shared<CBatchSQL>(1000, m_dbmanager, false);
    for(int i = 0; i < archiveRes.size(); i ++)  {
        CAnnSpectra *p = archiveRes.get(i);
        batch_sql->append(p->createUpdateNeighborSQL());

    }
}

string CAnnotationDB::toJsonNode(CDBEntry &dbentry, int i) {
//    m_tablename2header["GROUNDTRUTH"] = {"ID","FILEID","MS2COUNTS","PEPTIDE","SCORE","SCAN","CTERM","NTERM",
//    "MODIFICATION","PRECURSOR","CHARGE", "RT","PEPTIDEPROPHETPROB","IPROPHETPROB","RFSCORE",
//    "ISDECOY","SIGNIFICANCE","PROTEIN","CE","ALTERPEPTIDE","NEIGHBOR"};

    ostringstream oss;
    oss << "{"
        << R"("id": ")" << dbentry.get("ID",i)
        << R"(", "peptide": ")" << dbentry.get("PEPTIDE",i) //m_results[PEPTIDE]
        << R"(", "filename": ")" << dbentry.get("FILENAME",i) //m_results[FILENAME_APPENDED]
        << R"(", "scan": ")" << dbentry.get("SCAN",i) //m_results[SCAN]
        << R"(", "precursor": )" << dbentry.get("PRECURSOR",i) // m_results[PRECURSOR]
        << R"(, "charge": )" << dbentry.get("CHARGE",i) //m_results[CHARGE]
        << R"(, "cterm": ")" << dbentry.get("CTERM",i) // m_results[CTERM]
        << R"(", "nterm": ")" << dbentry.get("NTERM",i) //m_results[NTERM]
        << R"(", "othermod": ")" << dbentry.get("MODIFICATION",i) // m_results[MODIFICATION]
        << R"(", "rt": )" << dbentry.get("RT",i) // m_results[RT]
        << R"(, "score": ")" << dbentry.get("SCORE",i) // m_results[SCORE]
        << R"(", "pProb": )" << dbentry.get("PEPTIDEPROPHETPROB",i) // m_results[PEPTIDEPROPHETPROB]
        << R"(, "iProb": )" << dbentry.get("IPROPHETPROB",i) // m_results[IPROPHETPROB]
        << R"(, "isDecoy": )" << dbentry.get("ISDECOY",i) // m_results[ISDECOY]
        << R"(, "protein": ")" << dbentry.get("PROTEIN",i) // m_results[PROTEIN]
        << R"(", "significance": )" << dbentry.get("SIGNIFICANCE",i) // m_results[SIGNIFICANCE]
        << "}" << endl;
    string jsonnode = oss.str();
    return jsonnode;
}

string CAnnotationDB::toJsonNodes(CDBEntry &dbentry) {
    string message="";
    for(int i = 0; i < dbentry.size(); i ++){
        if(not message.empty())
        message += ",\n";

        message += toJsonNode(dbentry,i);
    }
    return message;
}

void CAnnotationDB::fixAllChargeState() {
    SPsmAnnotation gtinfo;
    long totalSpecNum = getTotalSpecNum();
    auto batchSQL = buildBatchSQL(1000,true);
    string sql = "select ID from GROUNDTRUTH where charge=0 and peptide!='UNKNOWN'";
    CDBEntry dbentry;
    m_dbmanager->getRow(dbentry,sql, false);

    if(dbentry.empty()) {
        return;
    }

    cout << "charge states to be fixed: " << dbentry.size() << endl;

    Progress ps(dbentry.size(),"fixing charge states");
    CountFrequency cf;
    for (long i = 0; i < dbentry.size(); i ++ ){
        ps.increase();
        retrieveGtinfo(dbentry.getInt("ID",i),gtinfo);
        if(gtinfo.fixChargeStates()){
            cf.add_data(gtinfo.charge);
            batchSQL->append(gtinfo.createChargeUpdateSql());
        }
    }
    cf.print();

}

// exact search, instead of searching for pattern.
shared_ptr<CDBEntry> CAnnotationDB::searchPeptide(const string &peptide) {
    string sql = "select * from GROUNDTRUTH where PEPTIDE = \"" + toSqlSafe(peptide) + "\" limit 10;";
    cout << "[Info] searching database with sql: " << sql << endl;

    shared_ptr<CDBEntry> dbentry = make_shared<CDBEntry>();
    m_dbmanager->getRow(*dbentry,sql,false);

    addFileNameToResultsOfGTRecords(*dbentry);
    return dbentry;
}

void CAnnotationDB::deleteLastNFile(int n) {
    if (n <= 0) {
        return;
    }
    int FileNum = getTotalFileNum();
    if (n > FileNum) {
        n = FileNum;
    }
    int lastFileId = FileNum - n;


    string sql = "delete from GROUNDTRUTH where FILEID >= " + to_string(lastFileId) + ";";
    m_dbmanager->execAsTransaction(sql, true);
    sql = "delete from SPECFILES where FILE_ID >= " + to_string(lastFileId) + ";";
    m_dbmanager->execAsTransaction(sql, true);

    cout << "DB updated" << endl;
}

void CAnnotationDB::connectDatabase(bool rebuild, string dbfilename, bool verbose) {
    setDB(dbfilename);
    createTables(rebuild, verbose);

}
/*
 *
    m_tablename2header["TNN"] = {"ID", "QUERY_ID", "TOPK", "MINDP", "INDEXNUM", "ARCHIVESIZE", "NPROBE", "DATE", "TNN_IDX", "TNN_DP"};
    m_createTableSql["TNN"] = "CREATE TABLE TNN("
                                      "ID                 INTEGER PRIMARY KEY AUTOINCREMENT ,"
                                      "QUERY_ID                INT,"
                                      "TOPK             INT,"
                                      "MINDP               FLOAT,"
                                      "INDEXNUM                 INT,"
                                      "ARCHIVESIZE                  INT,"
                                      "NPROBE                 INT,"
                                      "DATE                 CHAR(100),"
                                      "TNN_IDX               varchar, "
                                      "TNN_DP                varchar);"
                                      "CREATE INDEX TNN_QUERY_ID ON TNN (QUERY_ID);";*/
void CAnnotationDB::insertTnnInfo(long queryindex, int topk, double mindp, int indexnum, long archivesize, int nprobe, string date, vector<long> &tnnidx, vector<double> &tnndp) {
    ostringstream oss;
    oss << queryindex << "," << topk << "," << mindp << ",'" << indexnum
        << "'," << archivesize << "," << nprobe << ",'" << date << "',";
    oss << "'";
    for (int i=0; i < tnnidx.size(); i ++){
        oss << tnnidx[i];
        if(i!=tnnidx.size()-1) oss << " ";
    }
    oss << "','";
    for (int i=0; i < tnndp.size(); i ++){
        oss << tnndp[i];
        if(i!=tnnidx.size()-1) oss << " ";
    }
    oss << "'";


//    << nterm_mass
//        << ",'" << modificationstr << "'," << precursorMz << "," << charge << "," << retentiontimeinsec
//        << "," << pProb << "," << iProb << "," << isDecoy << "," << significance << ",'" << protein
//        << "','" << m_collision_energy << "'," << rfscore;

    string sql = "INSERT INTO TNN (QUERY_ID, TOPK, MINDP, INDEXNUM, ARCHIVESIZE, NPROBE, DATE, TNN_IDX, TNN_DP) "
                 "VALUES (" + oss.str() + "); ";
    m_dbmanager->execAsTransaction(sql, true);
//    return sql;
}

bool CAnnotationDB::getTNN(long queryindex, int topk, double mindp, int indexnum, long archivesize, int nprobe,
                           vector<long> &tnnidx, vector<double> &tnndp) {
    CDBEntry dbentry;
    ostringstream  oss;
    oss << "select * from TNN where QUERY_ID = " << queryindex
    // the following parameters does not change the TNN
//    << " and INDEXNUM =" << indexnum << " and nprobe = " << nprobe

    << " and ARCHIVESIZE = " << archivesize
    << " and TOPK >=  " << topk << " and MINDP <=" << mindp  <<";" ;
    string sql = oss.str();
//    cout << sql << endl;
    m_dbmanager->getRow(dbentry, sql, false);
    if (dbentry.empty()) {
//        cout << "No results found in db " << endl;
        return false;
    }
    else{
        string tnnidxstr = dbentry.get("TNN_IDX",0);
        string tnndpstr = dbentry.get("TNN_DP", 0);

        istringstream  iss(tnnidxstr);
        while(iss){
            int tnnid;
            iss >> tnnid;
            tnnidx.push_back(tnnid);
        }

        istringstream  iss2(tnndpstr);
        while(iss2){
            double dp;
            iss2 >> dp;
            tnndp.push_back(dp);
        }

        // now cut the dp threshold.
        double MAX_SCORE=42925;
        for(int i = 0; i < tnndp.size(); i ++){
            if(tnndp[i]/MAX_SCORE>=mindp and i<topk){

                continue;
            }else{

                tnndp.erase(tnndp.begin() + i, tnndp.end());
                tnnidx.erase(tnnidx.begin() + i, tnnidx.end());
                break;
            }
        }

    }
    return true;
}


CGtUpdater::CGtUpdater(CAnnotationDB *annodb, string gtfile) {
    m_annotationdb = annodb;
    m_gtfile = gtfile;
}

int CGtUpdater::getFileId(const string& eachfile) {
    cout << "to be implemented" << endl;
    return -1;
}

// This function looks weird...
string CGtUpdater::createUpdateSqlOfRow(long i, ICGtInfoUpdate &gt) {
    SPsmAnnotation gtinfo;
    string updateSql;
    if (m_annotationdb->retrieveGtinfo(i, gtinfo) and gt.updateGtInfo(gtinfo)) {
        updateSql = gtinfo.createUpdateSql();
    }
    return updateSql;
}


void CGtUpdater::update(const string& rawfilename, ICGtInfoUpdate &gt) {
    vector<vector<string>> results;
    m_annotationdb->getSpecFileRows(rawfilename, results);

    if (results.empty() or results.size() > 1) {
        cout << "No raw file found with name: " << rawfilename << endl;
    } else if (results.size() > 1) {
        cout << results.size() <<" raw files found."<< endl;
    } else {
        specfileinfo sfinfo(results[0]);
        int batchSize = 10000;
        shared_ptr<CBatchSQL> batchSql = m_annotationdb->buildBatchSQL(batchSize, false);

        int nan_counts = 0;
      //  Progress ps(sfinfo.end - sfinfo.start, string("Updating ground truth table with ") + m_gtfile);
        for (long i = sfinfo.start; i < sfinfo.end; i++) {
          //  ps.increase();
            string update_sql = createUpdateSqlOfRow(i, gt);
            if (update_sql.empty()) {
                nan_counts++;
                continue;
            }
            batchSql->append(update_sql);
        }

        spdlog::get("A")->info("{}, PSM annotation updated:  {} / {} ~ {:.4f} %", rawfilename, sfinfo.end - sfinfo.start - nan_counts, sfinfo.end - sfinfo.start, 100 - nan_counts * 100.0 / (sfinfo.end - sfinfo.start));

    }

}

// repeated code compared to the factoryMethod
CGtUpdater::~CGtUpdater() {
    string ext = File::CFile(m_gtfile).ext;
    if (ext == "xml" and m_gtfile.find("interact-") != string::npos and
        m_gtfile.substr(m_gtfile.length() - 7).find("pep.xml") != string::npos) {
        cout << "[Info] Found xinteract result, PeptideProphet or iProphet result are both supported" << endl;
        PeptideProphetParser ppp(m_gtfile);

        for (const auto& each_file : ppp.m_allSourceFiles) {
            update(each_file, ppp);
        }

    } else if (m_gtfile.substr(m_gtfile.length() - 7).find("pep.xml") != string::npos
               or ext == "pepXML") {

        CometPepXMLParser crp(m_gtfile);

        string mzxml_name = File::CFile(m_gtfile).basename;

        if (ext == "pepXML") {
            // cout << "Got PepXML file" << endl;
            mzxml_name = m_gtfile;
        }
        update(mzxml_name, crp);
    } else if (ext == "mgf" or ext == "sptxt") {
        DataFile lib_file(m_gtfile, 0, -1);
        AnnotationDataFile adf(lib_file);
        update(m_gtfile, adf);
    } else if (ext == "all" and m_gtfile.find(".psm.all") != string::npos) {
        // support tsv table.
        cout << "Find new gt file type" << endl;
        // find csv file
        CTable psm(m_gtfile, ',', true, 0);
        // cout << "table loaded" << endl;
        CsvAnnotation anno(psm);
        unordered_set<string> s;
        for (int i = 0; i < psm.m_row; i++) {
            s.insert(psm.getEntry(i, psm.getColByHeader("filename")));
        }
        for (const auto & it : s) {
            //TODO: the following position to be changed.
            string filename_prefix = "/nasbackup/wulong_backup/pxd000561/";
            string rawfilename = filename_prefix + File::CFile(it).filename;
            cout << "[Info] processing file " << rawfilename << endl;
            update(rawfilename, anno);

        }
    } else if (ext == "tsv" and m_gtfile.find(".spectroscape.tsv") != string::npos) {
        // support new spectroscape.tsv format.
        // support tsv table.
        // cout << "Find new gt file type" << endl;
        // find csv file
        CTable psm(m_gtfile, '\t', true, 0);
        // cout << "table loaded" << endl;
        TsvAnnotation anno(psm);
        unordered_set<string> s;
        for (int i = 0; i < psm.m_row; i++) {
            s.insert(psm.getEntry(i, psm.getColByHeader("filename")));
        }
        for (const auto & it : s) {
            //TODO: the following position to be changed.
            // string filename_prefix = "/nasbackup/wulong_backup/pxd000561/";
            string rawfilename =  File::CFile(it).filename;
            cout << "[Info] processing file " << rawfilename << endl;
            update(rawfilename, anno);

        }
    }else {
        cout << "Unknown ground truth format! doing nothing " << endl;
    }
}

bool CBlackList::isBlack(CBlackList::dataType dt, long id) {
    bool ret = false;
    switch (dt) {
        case MAP:
            ret = isBlackMap(id);
            break;
            // the two function below saves time!
        case UNORDERED_MAP:
            ret = isBlackUnorderedMap(id);
            break;
        case UNORDERED_SET:
            ret = isBlackUnorderedSet(id);
            break;
        default:
            cout << "[isBlack] Error data type " << dt << endl;

    }
    return ret;
}

bool CBlackList::isBlackUnorderedSet(long id) {
    if (m_blackListUnorderedSet.empty()) {
        for_each(m_blackList.begin(), m_blackList.end(), [&](const long &x) { m_blackListUnorderedSet.insert(x); });
    }
    return m_blackListUnorderedSet.find(id) != m_blackListUnorderedSet.end();
}

bool CBlackList::isBlackUnorderedMap(long id) {
    if (m_blackListUnorderedMap.empty()) {
        for (auto const &x: m_blackList) {
            m_blackListUnorderedMap[x] = true;
        }
    }
    return m_blackListUnorderedMap.find(id) != m_blackListUnorderedMap.end();
}

bool CBlackList::isBlackMap(long id) {
    if (!m_blackList.empty() and m_blackListmap.empty()) {
        for (long & i : m_blackList) {
            m_blackListmap[i] = true;
        }
    }

    return 0 != m_blackListmap.count(id);
}
