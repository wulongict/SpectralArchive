//
// Created by wulong on 11/25/17.
//

#ifndef MYTOOL_PEAK_ANNOTATION_H
#define MYTOOL_PEAK_ANNOTATION_H


#include <string>
#include <vector>
#include <iostream>
using namespace std;

class peak_annotation {
public:
    peak_annotation();
    peak_annotation(const peak_annotation &other);

    void print() const;

    static void print(vector<peak_annotation> &pka);

    string ion_base_type;
    string ion_NL;
    bool isotopic;
    int pos;
    int charge;
    double masserror;
};


void parse_annotation_as_struct(string annotation, vector<peak_annotation> &pka);
void parse_annotation_as_structX(string annotation, vector<peak_annotation> &pkas);


#endif //MYTOOL_PEAK_ANNOTATION_H
