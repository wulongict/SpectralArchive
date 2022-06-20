//
// Created by wulong on 11/25/17.
//

#include "FragmentationScore.h"
#include "gnuplot_functions.h"
#include "classifier.h"
#include "CPeakPairsImporter.h"
#include "CFragScore.h"
#include "PeakPairScore.h"
#include "PSMScoreCalculator.h"
#include "CPeakPair.h"

void calculatePSMscore(const string &scoretype, CPeakPairsImporter &sample_Data, classifier &lr) {
    cout << "[Info] Score type: " << scoretype << endl;
    sample_Data.print();
    CFragScore * result = lr.predict(sample_Data);
    FragmentationScore PSMscore(result);
    string scorefilename = PSMscore.compute(scoretype, sample_Data);
}

FragmentationScore::FragmentationScore(CFragScore *predictionFile) {
    string inputsamplefile = predictionFile->getM_filename();
    m_samplefilename = inputsamplefile.substr(0,inputsamplefile.find_last_of('.'));
    m_predicted_result = inputsamplefile;
    m_psmid_filename = m_samplefilename + ".psm_id";
    cout << m_predicted_result << endl << m_psmid_filename << endl;
    load_ground_truth();  // groundtruth
    load_psmid();  // psm id and sample id
    m_predictionFile = predictionFile;
}

FragmentationScore::FragmentationScore(CPeakPairsImporter &data, string &scoreType, classifier &lr) {
    m_samplefilename = data.getOutFileStr();
    m_predictionFile = lr.predict(data);
    for(int i = 0; i < data.size(); i ++)    {
        CPeakPair x = data.get_by_index(i);
        m_spectrum_id.push_back(x.m_PSM_ID);
        m_SampleID.push_back(x.m_sample_id);
        m_groundtruthlabel.push_back(x.m_y >0 ? 1.0: -1.0);
    }
    groupSampleByPSM(data);
    calcScore(scoreType);
}

FragmentationScore::FragmentationScore(const FragmentationScore &other) {
    cout << "Error we are using undefined function" << endl;
    exit(0);
}

FragmentationScore &FragmentationScore::operator=(const FragmentationScore &other) {
    cout << "Error we are using undefined function" << endl;
    exit(0);
}

void FragmentationScore::load_psmid() {
    ifstream fin;
    fin.open(m_psmid_filename, ios::in);
    if (!fin) {
        cout << "Fail to read psm id file: " << m_psmid_filename << endl;
        exit(0);
    }
    int psmid, sample_id;
    while (fin >> psmid >> sample_id) {
        m_spectrum_id.push_back(psmid);
        m_SampleID.push_back(sample_id);
    }
    cout << "Number of PSM loaded : " << m_spectrum_id.size() << endl;
    fin.close();
}

void FragmentationScore::load_ground_truth() {
    ifstream fin;
    fin.open(m_samplefilename, ios::in);
    if (!fin) {
        cout << "Fail to read ground truth file: " << m_samplefilename << endl;
        exit(0);
    }
    string line;
    while (getline(fin, line)) {
        double label = atof(line.substr(0, line.find_first_of(' ')).c_str());
        m_groundtruthlabel.push_back(label);
    }
    cout << "Number of ground truth label loaded " << m_groundtruthlabel.size() << endl;
    fin.close();
}


void FragmentationScore::print() {
    cout << m_samplefilename << endl;
    cout << "PSM num " << m_PSMscore.size() << endl;
    if (m_PSMscore.empty()) {
        cout << "Empty score" << endl;
        exit(0);
    }
    m_PSMscore[0].print();
}


void FragmentationScore::calcScore(const string& scoretype) {
    for (auto & i : m_PSMscore) {
        i.compute(scoretype);
    }

    vector<double> scores;
    for (auto eachpsm : m_PSMscore) scores.push_back(eachpsm.getcombinedscore());
    bool plotScoreIntenFC=false;
    if(plotScoreIntenFC)
    {
        cout << "plot score and intensity fc score" << endl;
        for (int i = 0; i < m_PSMscore.size() && i < 1; i++) {
            vector<vector<double>> fc_score;
            m_PSMscore[i].collect(fc_score);
            scatterplot(scoretype + ": " + m_samplefilename, fc_score);
        }
    }
}


void FragmentationScore::groupSampleByPSM(CPeakPairsImporter &input) {
    int i = 0;
    int spectrum_id = -1;
    while (i < m_groundtruthlabel.size()) {
        if (m_spectrum_id.at(i) != spectrum_id) {
            spectrum_id = m_spectrum_id.at(i);
            PSMScoreCalculator tmp;
            m_PSMscore.push_back(tmp);
        }

        if (m_PSMscore.empty()) {
            cout << "Empty score vector " << endl;
            exit(0);
        }

        m_PSMscore[m_PSMscore.size() - 1].setPSMId(spectrum_id);
        PeakPairScore score = m_predictionFile->getScore(i);

        CPeakPair ts = input.linearSearch(m_SampleID.at(i));
        m_PSMscore[m_PSMscore.size() - 1].append(score, m_groundtruthlabel.at(i), ts);

        if (ts.getSampleID() == m_SampleID.at(i)) {
        } else {
            cout << "Error, incorrect index" << endl;
            exit(0);
        }
        i++;
    }
}


FragmentationScore::~FragmentationScore() {
    if (m_predictionFile != nullptr)    {
        delete m_predictionFile;
        m_predictionFile = nullptr;
    }
}

double FragmentationScore::getscore(int psmid) {
    double score = 0;
    for (auto & i : m_PSMscore)    {
        if(i.getPSMId() == psmid)        {
            score = i.getcombinedscore();
            break;
        }
    }
    return score;
}

string FragmentationScore::compute(const string &scoretype, CPeakPairsImporter &SampleData) {
    groupSampleByPSM(SampleData);
    calcScore(scoretype);
    string scorefilename = exportScore(scoretype);
    return scorefilename;
}

string FragmentationScore::exportScore(const string& scoretype) {
    string outfilename = m_samplefilename + "_" + scoretype + ".psmscore";
    cout << "[Info] Exporting PSM scores to file: \n" << outfilename << endl;
    ofstream fout;
    fout.open(outfilename.c_str(), ios::out);
    if (!fout) {
        cout << "Error: fail to open file: " << outfilename << endl;
        exit(0);
    }
    for (auto x: m_PSMscore) fout << x.getPSMId() << " " << x.getcombinedscore() << endl;

    fout.close();
    return outfilename;
}





