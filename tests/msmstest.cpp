#include <iostream>
#include "../librarymsms/ProteomicsDataTypes.h"
#include "../librarymsms/Util.h"
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


