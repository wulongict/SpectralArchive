//
// Created by wulong on 10/20/19.
//

#ifndef MYTOOL_RANDOMSPECPARSER_H
#define MYTOOL_RANDOMSPECPARSER_H

//
// Created by wulong on 10/20/19.
//


#include <vector>
#include <string>
using namespace std;


//bool getPeakListFromSplib(string filename, int scannum, vector<double> &mz, vector<double> & intensity);
bool getPeakList(string filename, int scannum, vector<double> &mz, vector<double> &intensity);


#endif //MYTOOL_RANDOMSPECPARSER_H
