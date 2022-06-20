//
// Created by wulong on 8/29/18.
//

#include "commonfun.h"
#include <iostream>
using namespace std;

void displayParam(boost::program_options::variables_map &vm) {
    namespace po = boost::program_options;
    spdlog::get("A")->info("-------------Parameter-----Start---------");
    for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++) {
        if (it->second.value().type() == typeid(int)) {
            spdlog::get("A")->info("{} = {}", it->first, boost::any_cast<int>(it->second.value()));
        } else if (it->second.value().type() == typeid(long)) {
            spdlog::get("A")->info("{} = {}", it->first, boost::any_cast<long>(it->second.value()));
        } else if (it->second.value().type() == typeid(double)) {
            spdlog::get("A")->info("{} = {}", it->first, boost::any_cast<double>(it->second.value()));
        } else if (it->second.value().type() == typeid(string)) {
            spdlog::get("A")->info("{} = {}", it->first, boost::any_cast<string>(it->second.value()));
        } else if (it->second.value().type() == typeid(bool)) {
            spdlog::get("A")->info("{} = {}", it->first, boost::any_cast<bool>(it->second.value()));
        } else {
            spdlog::get("A")->info("{} : {} ", it->first, "an option of an unknown type!");
        }
    }
    spdlog::get("A")->info("-------------Parameter-----End---------\n");
}

void uniqueVector(vector<long> &x) {
    sort(x.begin(), x.end());
    vector<long>::iterator it = std::unique(x.begin(), x.end());
    x.resize(std::distance(x.begin(), it));
}

void stableUniqueVector_deprecated(vector<long> &x, bool verbose) {
    vector<int> sizeOfIdx; //
    sizeOfIdx.push_back(x.size()); // total number of ids from N indexes

    map<long, bool> exist={{-1, true}};
    vector<long> y;
    int count = 0;
    for(int i = 0; i < x.size(); i ++)   {
        if(0==exist.count(x[i]))   {
            exist[x[i]] = true;
            y.push_back(x[i]);
        }
        if(i>= count+1024){
            sizeOfIdx.push_back(y.size());
            count += 1024;
        }
    }
    x.swap(y);

    sizeOfIdx.push_back(x.size());

    if(verbose){
        cout << "MultipleIndex";
        for(auto eachSize: sizeOfIdx){
            cout << "\t" << eachSize;
        }
        cout << endl;
    }
}

void stableUniqueVectorSet(vector<long> &x, bool verbose, const int sizeOfAnnByEachIndex) {
    set<long> s = {-1};
    vector<int> sizeOfIdx; //
    vector<long> y;
    sizeOfIdx.push_back(x.size());

    for (int k = 0; k < x.size(); k++) {
        if (s.find(x[k]) == s.end()) {
            y.push_back(x[k]);
            s.insert(x[k]);
        }
        if ( (k+1) % sizeOfAnnByEachIndex == 0 ) {
            sizeOfIdx.push_back(y.size());
        }
    }
//    sizeOfIdx.push_back(y.size());
    x.swap(y);

    if (verbose) {
        cout << "MultipleIndex";
        for (auto eachSize: sizeOfIdx) {
            cout << "\t" << eachSize;
        }
        cout << endl;
    }
}

string argToStr(int argc, char **argv) {
    string abc(argv[0]);
    for(int i = 1; i < argc; i ++)  {
        abc += " ";
        abc += argv[i];
    }
    return abc;
}
