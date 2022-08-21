//
// Created by wulong on 6/19/19.
//

#include "CMyIndex.h"
#include "CKMeans.h"
#include <iostream>

using namespace std;

//CMyIndex::CMyIndex(CPQParam option) {
//    cpq = nullptr;
//    m_option = option;
//    cout << "inside index wrapper " << m_option._option << endl;
//}

CMyIndex::CMyIndex(shared_ptr<CPQParam> optionptr) {
    cpq = nullptr;
    m_option = optionptr;
    cout << "inside index wrapper " << m_option->_option << endl;
}

//CMyIndex::CMyIndex() :CMyIndex(CPQParam()){
//
//}

CMyIndex::~CMyIndex() {}

void CMyIndex::read(string filename) {
    cpq = make_shared<CProductQuantization>();
    cpq->load(filename);
}

void CMyIndex::add(long newspecnum, float *vec) {
    cpq->add(newspecnum, vec, dim());

}

void CMyIndex::display() {
    cpq->display(false);
}

void CMyIndex::write(string filename) {
    cpq->save(filename);
}

long CMyIndex::total() {
    return cpq->size();
}

void CMyIndex::search(int num, float *queries, int ret_num, vector<float> &dist, vector<long> &ind) {
    int d = dim();
    CPQResult res(ret_num);
    cpq->search(queries, num, res);

    res.flattern(dist, ind, d);
}

void CMyIndex::setnProb(int nprobe) { cpq->setnProbe(nprobe); }

bool CMyIndex::empty() { return nullptr == cpq; }

bool CMyIndex::istrained() { return cpq->isTrained(); }

void CMyIndex::createEmptyIndex(int dim, string indexstr) {
    //cout << "Creating Empty Index " << indexstr << endl;
    int found = indexstr.find(",");
    if (found != string::npos) {
        string ivf = indexstr.substr(0, found), pqn = indexstr.substr(found + 1);
        int clusterNum = atoi(ivf.substr(3).c_str());
        int subspaceNum = atoi(pqn.substr(2).c_str());
        cout << "PQ-" << subspaceNum << "; K = " << clusterNum << endl;
        cpq = make_shared<CProductQuantization>(subspaceNum, clusterNum, dim);
    } else {
        cout << "Error in index str: " << indexstr << endl;
    }

}

void CMyIndex::train(int batchSize, float *vecBatch) {
    int d = dim();
    int K = cpq->getClusterNum();

    CDataFloatPtr dfp(vecBatch, batchSize, d);
    cpq->train(&dfp, *m_option);
}

int CMyIndex::dim() { return cpq->getDim(); }

void CMyIndex::toCPU() {}

void CMyIndex::toGPU() {}

void CMyIndex::removeIds(vector<long> &idx) {
    cout << "under construction " << endl;
}
