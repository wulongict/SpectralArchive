#include <iostream>
#include "../librarymsms/ProteomicsDataTypes.h"
#include "../librarymsms/Util.h"
#include "../librarymsms/DatabaseManager.h"
#include "../ArchiveSearch/dependentlibs/randomSpecParser.h"
#include "gtest/gtest.h"
#include "CAnnotationDB.h"

using namespace std;

TEST(CASE1, EQUAL_INT){
    ASSERT_EQ(1,1);
    cout << "Hello world" << endl;
}

TEST(MSMS_TEST, PEAK_REMOVAL){
    // getAllPeaks have several different functionalities.
    // 1. remove low intense peaks. noise peaks.
    // 2. remove isotopic peaks.
    // 3. remove TMT reporter ions
    // 4. remove immonium ions. 

    CSpectrum cs;
    cs.addOnePeak(100,200);
    cs.addOnePeak(161,200);
    vector<double> mz, intensity;
    cs.getAllPeaks(mz, intensity, true, true, true, 10);
    
    ASSERT_EQ(mz.size(), 1);
    
}

TEST(SPTXT_PARSER, RANDOM_ACCESS){
    string filename = "../tests/data/abc.sptxt";
    vector<double> mz, intensity;
    getPeakList(filename, 0,mz, intensity);
    cout << "Done" << endl;
    ASSERT_FLOAT_EQ(102.0552, mz[1]);
    ASSERT_FLOAT_EQ(	232.7415, intensity[1]);

}

TEST(MGFPARSER, RANDOM_ACCESS){
    string filename = "../tests/data/iPRG2012_first_3scans.mgf";
    vector<double> mz, intensity;
    getPeakList(filename, 2,mz, intensity);
    cout << "peak number " << mz.size() << " " << intensity.size() << endl;
    cout << mz[0] << "\t" << intensity[0] << endl;
    cout << "Done" << endl;
//    989.6261597 72.00098419
    ASSERT_FLOAT_EQ(989.6261597,mz[1]);
    ASSERT_FLOAT_EQ(72.00098419, intensity[1]);

}

TEST(DB, NEW_COL){
    cout << "testing new db" << endl;
    shared_ptr<CAnnotationDB> m_AnnotationDB = make_shared<CAnnotationDB>(false);
    m_AnnotationDB->connectDatabase(false, "abc.sqlite3db",true);
}


// creat a test of CTable to read a text file in data folder. 
// the file name is new_search_result_file_list.txt
TEST(CTABLE, READ_TXT){
    string pepxmlfilelist = "../tests/data/new_search_result_file_list.txt";
    CTable pepxmllist(pepxmlfilelist, '\t', false, 0);
    
    ASSERT_EQ(1,pepxmllist.m_row);
    
}


// memory footprint of data types. 
TEST(MEMORY_USAGE, DATATYPES){
    cout << "size of int " << sizeof(int) << endl;
    cout << "size of long " << sizeof(long) << endl;
    cout << "size of float " << sizeof(float) << endl;
    cout << "size of double " << sizeof(double) << endl;
    cout << "size of long long " << sizeof(long long) << endl;
    cout << "size of long double " << sizeof(long double) << endl;
    cout << "size of size_t " << sizeof(size_t) << endl;
    cout << "size of shared_ptr<int> " << sizeof(shared_ptr<int>) << endl;
    cout << "size of unsigned int " << sizeof(unsigned int) << endl;
    cout << "size of unsigned long " << sizeof(unsigned long) << endl;
    cout << "size of unsigned long long " << sizeof(unsigned long long) << endl;
    
    

}
// testing speed of file search
TEST(FILE_SEARCH, SPEED_TEST){
    // this test is disabled now. It will run only when the file is available.
    // string dbfile = "../tests/data/arxiv_pxd000561.sqlite3db";
    
    string dbfile = "../tests/data/_arxiv_pxd000561.sqlite3db";
    if(!File::isExist(dbfile)){
        cout << "file not exist" << endl;
        return;
    }
    shared_ptr<CAnnotationDB> m_AnnotationDB = make_shared<CAnnotationDB>(false);
    m_AnnotationDB->connectDatabase(false, dbfile,true);
    shared_ptr<CDBEntry> dbentry=nullptr;
    string filename;
    cout << "Please input a filename: ('Q' to quit): " <<endl << flush;
    cin >> filename;
    while(filename !="Q") {
         m_AnnotationDB->searchGTWithFileName_new(filename,"1","1000",dbentry);
        // dbentry->print();
        cout << "---> entries found in current method" << dbentry->size() << endl;;
        cout << "---> Please input a filename: ('Q' to quit)" << endl << flush;
        cin >> filename;
    }
   



}

