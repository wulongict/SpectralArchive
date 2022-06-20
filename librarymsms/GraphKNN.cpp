//
// Created by wulong on 2/18/19.
//

#include "GraphKNN.h"
#include "Util.h"

#include <iostream>

using namespace std;

GraphKNN::GraphKNN(int numCandidates) : ncandidates(numCandidates) {
    m_knn = nullptr;
}

GraphKNN::~GraphKNN() {
    if (m_knn) delete[] m_knn;
}

bool GraphKNN::empty() {
    return m_knn == nullptr;
}

long *GraphKNN::getRow(long row) {
    return m_knn + row * ncandidates;
}

void GraphKNN::load_knn_graph_from_file(string knngraphfile) {
    long nSize = 0;
    File::getfilesize(knngraphfile, nSize);
    FILE *pfile = fopen(knngraphfile.c_str(), "rb");
    m_knn = new long[nSize / sizeof(long)];
    long numBytes = fread((char *) m_knn, sizeof(char), nSize, pfile);
    cout << "bytes loaded: " << numBytes << endl;
    fclose(pfile);
}
