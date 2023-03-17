#include <iostream>
#include "../librarymsms/ProteomicsDataTypes.h"
#include "../ArchiveSearch/dependentlibs/randomSpecParser.h"
#include "gtest/gtest.h"
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
    string filename = "/data/wulong/data/iPRG/NIST_human_hcd_all_2020_TargetDecoy.sptxt";
    vector<double> mz, intensity;
    getPeakList(filename, 0,mz, intensity);
    cout << "Done" << endl;

}

TEST(MGFPARSER, RANDOM_ACCESS){
    string filename = "/data/wulong/data/iPRG/iPRG2012.mgf";
    vector<double> mz, intensity;
    getPeakList(filename, 2,mz, intensity);
    cout << "peak number " << mz.size() << " " << intensity.size() << endl;
    cout << mz[0] << "\t" << intensity[0] << endl;
    cout << "Done" << endl;

}
