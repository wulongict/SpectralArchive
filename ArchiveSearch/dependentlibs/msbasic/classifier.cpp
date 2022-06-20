//
// Created by wulong on 11/25/17.
//

#include <iostream>
#include <utility>
#include <spdlog/spdlog.h>
#include "CPeakPairsImporter.h"
#include "CFragScore.h"
#include "PeakPairScore.h"
#include "liblinearWrapper.h"
#include "classifier.h"
#include "../../../librarymsms/Util.h"
#ifdef LIBLINEAR_MT
#include "../../../External/liblinear-multicore-2.42/linear.h"
#else
#include "../../../External/liblinear-2.11/linear.h"
#endif

#define INF HUGE_VAL



classifier::classifier(string modelFileName, string binaryPath)  {
    m_model_ptr = nullptr;
    m_binaryPath = binaryPath;
    cout << "binary path is " << binaryPath << " set to : " << m_binaryPath <<endl;
    set_linear_regression_tool_path();
    m_featureToFile = false;
    m_doTest = false;
    m_modelFileName = std::move(modelFileName);
    cout <<"[Info] Model name: " <<  m_modelFileName << endl;
}

// to be updated
void classifier::set_linear_regression_tool_path() {
#ifdef LIBLINEAR_MT
    m_trainingBinaryPath = m_binaryPath + "/" + "train_mt";
    m_predictBinaryPath = m_binaryPath + "/" + "predict_mt";
#else
    m_trainingBinaryPath = m_binaryPath+"/"  + "train";
    m_predictBinaryPath =m_binaryPath+"/"  + "predict";
#endif
    cout << "[Info] liblinear version : " << liblinear_version << endl;
    cout << "[Info] path of training and testing tools: \n" << m_trainingBinaryPath << endl << m_predictBinaryPath << endl;
}

classifier::~classifier() {
    if(m_model_ptr == nullptr){
        cout << "free model" << endl;
        free_and_destroy_model(&m_model_ptr);
    }
};

void classifier::trainModel(CPeakPairsImporter &trainingdataset, bool overwrite) {
    if (File::isExist(m_modelFileName, true) and not overwrite) {
        spdlog::get("A")->info("Fragmentation model {} exist! Skip training step!", m_modelFileName);
    } else if (m_featureToFile)  {
        // Not used now!
        cout << "[Info] training model " << m_modelFileName << endl;
        string trainingfile = trainingdataset.getOutFileStr() + "train.txt";
        cout << "[Info] Training data will be exported into " << trainingfile << endl;
        trainingdataset.exportTofile(trainingfile);
        string cmdline = m_trainingBinaryPath;
        if (liblinear_version == 242){
            cmdline += " -m 10 ";
        }
        cmdline = cmdline + " -s 0 ";
        cmdline = cmdline + trainingfile;
        cmdline = cmdline + " " + m_modelFileName;
        cout << "[Info]  $ " << cmdline << endl;
        int status = system(cmdline.c_str());
        cout << "Status: " << status << endl;
    } else {
        // todo: split training set into training and validation part.
        //  do cross validation
        cout << "[Info] prepare to train liblinear model" << endl;
        parameter param;
        make_default_param_liblinear(param);
        problem prob;
        initialize_problem(&prob, trainingdataset.size());
        trainingdataset.toLiblinearProb(&prob, true, m_oheAA);
        const char *error_msg;
        error_msg = check_parameter(&prob, &param);
        if (error_msg) {
            fprintf(stderr, "ERROR: %s\n", error_msg);
            exit(1);
        }
        cout << "Start to train the model of fragmentation " << endl;
        model * model_ = train(&prob, &param);
        if (model_ == nullptr) {
            cout << "invalid pointer model_: " << __FILE__ << ":" << __LINE__ << endl;
            exit(0);
        }
        save_model(m_modelFileName.c_str(), model_);
        cout << "model saved to : " << m_modelFileName << endl;
        free_and_destroy_model(&model_);
        destroy_param(&param);
        free_problem(&prob);
    }
}


double classifier::predict(model *model_, feature_node *x, vector<double> &probability) const{
    int nr_class = model_->nr_class;
    double *prob_estimates = nullptr;
    prob_estimates = (double *) malloc(nr_class * sizeof(double));
    if (prob_estimates == nullptr) {
        cout << "malloc fail for prob_estimates: " << __FILE__ << ":" << __LINE__ << endl;
        exit(0);
    }
    double predict_label = predict_probability(model_, x, prob_estimates);
    for (int i = 0; i < model_->nr_class; i++) {
        probability.push_back(prob_estimates[i]);
    }
    free(prob_estimates); // remember: do not use delete here!
    return predict_label;
}


CFragScore *classifier::predict(CPeakPairsImporter &peakPairs) {
    if(nullptr == m_model_ptr){
        predictReady();
    }
    if (m_model_ptr == nullptr) {
        cout << "model_ pointer is invalid " << __FILE__ << ":" << __LINE__ << endl;
        exit(0);
    }
    problem prob;
    initialize_problem(&prob, peakPairs.size());
    peakPairs.toLiblinearProb(&prob, false, m_oheAA);
    CFragScore *result = new CFragScore();

    if (nullptr == result) {
        cout << "invlid pointer result: " << __FILE__ << ":" << __LINE__ << endl;
        exit(0);
    }

    int len = prob.l;
    double accuracy = 0;
    for (int i = 0; i < len; i++) {
        vector<double> probability;
        double label = predict(m_model_ptr, prob.x[i], probability);

        PeakPairScore pps(label, probability.at(0), probability.at(1));
        result->add(pps);
        if ((label * prob.y[i]) > 0) {
            accuracy += 1;
        }
    }


    free_problem(&prob);

    if (len > 0 && m_doTest) {
        accuracy /= len;
        cout << "Testing on: " << peakPairs.getOutFileStr() << endl;
        cout << setprecision(4) << "[Info] accuracy is : " << accuracy * 100 << "%" << endl;
    }
    // todo: in fact, we do not need the pepxml file. we only need a list of peptide,
    //  which can comes from any search engine or post processing tool. e.g. percolator
    // todo: we need a pepxml and a threshold to tell us the truth.

    return result;
}


void classifier::setOutputMethod(bool outputFile) {
    m_featureToFile = outputFile;
}

void classifier::setTestingMode(bool testingMode) {
    m_doTest = testingMode;
}

void classifier::predictReady() {
    // load model
    lock_model_ptr.lock();
    if(nullptr == m_model_ptr){
        m_model_ptr = load_model(m_modelFileName.c_str());

    }
    lock_model_ptr.unlock();

}
