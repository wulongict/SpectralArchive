#include <iostream>
#include "../librarymsms/ProteomicsDataTypes.h"
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
