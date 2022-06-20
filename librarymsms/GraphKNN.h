//
// Created by wulong on 2/18/19.
//

#ifndef MYTOOL_GRAPHKNN_H
#define MYTOOL_GRAPHKNN_H

#include <string>


class GraphKNN {
    long *m_knn;
public:
    const int ncandidates;
    explicit GraphKNN(int numCandidates);
    ~GraphKNN();
    bool empty();
    long * getRow(long row);
    void load_knn_graph_from_file(std::string knngraphfile);
};


#endif //MYTOOL_GRAPHKNN_H
