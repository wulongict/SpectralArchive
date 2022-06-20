//
// Created by wulong on 8/29/18.
//

#ifndef MYTOOL_COMMONFUN_H
#define MYTOOL_COMMONFUN_H

#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>
#include <vector>

void displayParam(boost::program_options::variables_map &vm);
void uniqueVector(std::vector<long> &x);
void stableUniqueVector_deprecated(std::vector<long> &x, bool verbose);
void stableUniqueVectorSet(std::vector<long> &x, bool verbose,const int sizeOfAnnByEachIndex);
std::string argToStr(int argc, char * argv[]);
#endif //MYTOOL_COMMONFUN_H
