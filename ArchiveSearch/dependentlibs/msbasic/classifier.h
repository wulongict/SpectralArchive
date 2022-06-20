//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_LIBLINEARTOOL_H
#define MYTOOL_LIBLINEARTOOL_H

#include <string>
#include <mutex>
#include "CPeakPair.h"
using namespace std;

class CPeakPairsImporter;
class CFragScore;
class model;
class feature_node;
class COneHotEncodingAA;



class classifier {
private:
    COneHotEncodingAA m_oheAA;
    bool m_doTest;
    string m_modelFileName;
    string m_trainingBinaryPath;
    string m_predictBinaryPath;
    string m_binaryPath;
    bool m_featureToFile;
    model *m_model_ptr;
    mutex lock_model_ptr;
    void set_linear_regression_tool_path();
public:
    explicit classifier(string modelFileName, string binaryPath);
    ~classifier();
    void setOutputMethod(bool outputFile);
    void setTestingMode(bool testingMode);
    void trainModel(CPeakPairsImporter &trainingdataset, bool overwrite);
    CFragScore *predict(CPeakPairsImporter &peakPairs);
    void predictReady();
private:
    // predict single instance.
    double predict(model *model_, feature_node *x, vector<double> &probability) const;
};



#endif //MYTOOL_LIBLINEARTOOL_H
