//
// Created by wulong on 3/13/21.
//

#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <vector>
#include <memory>
#include "CFlow.h"
#include "CDebugMode.h"

using namespace std;

TEST(ROC_,AUC_){
    initlog("x.log", "A");
    std::vector<double> positiveScores, negativeScores;
    positiveScores.push_back(0.1);
    negativeScores.push_back(0.2);
    positiveScores.push_back(0.2);
    negativeScores.push_back(0.3);
    //
    //  pos     0.1     0.2
    //  neg             0.2     0.3
    //  --------------------------------
    //  tpr    0.5      1.0     1.0
    //  fpr    0        0.5     1.0
    //  -------------------------------
    //
    /* AUC 0.875
     * ROC curve
     *                0.5______1.0
     *      1.0 |   /   |      |
     *          |_/     |      |
     *      0.5 |       |      |
     *          |       |      |
     *        0  -----------------
     *          0     0.5     1.0
     * */
    //
    //


    CROCPlot roc( negativeScores,positiveScores),rocSwitch( positiveScores,negativeScores);

    ASSERT_EQ(roc.getAUC(),7/8.0);
    ASSERT_EQ(rocSwitch.getAUC(),1/8.0);


}