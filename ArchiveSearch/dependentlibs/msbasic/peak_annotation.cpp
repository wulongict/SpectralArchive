//
// Created by wulong on 11/25/17.
//

#include "peak_annotation.h"
#include <cmath>

peak_annotation::peak_annotation() {
    ion_base_type = "";
    ion_NL = "0";
    isotopic = false;
    pos = -1;
    charge = 1;
    masserror = 0;
}

void peak_annotation::print() const {
    cout << "type: " << ion_base_type << " nl:" << ion_NL
         << " p: " << pos << " c: " << charge << "+"
         << " iso: " << isotopic << " Δm: " << masserror << endl;
}

void peak_annotation::print(vector <peak_annotation> &pka) {
    for (auto each: pka) each.print();
}

peak_annotation::peak_annotation(const peak_annotation &other) {
    ion_base_type = other.ion_base_type;
    ion_NL = other.ion_NL;
    isotopic = other.isotopic;
    pos = other.pos;
    charge = other.charge;
    masserror = other.masserror;
}

void parse_annotation_as_struct(string annotation, vector<peak_annotation> &pka) {
    int i = 0;
    while (i < annotation.length()) {
        peak_annotation peak_anno;
        if (annotation[i] == 'b') {
            peak_anno.ion_base_type = "b";
        } else if (annotation[i] == 'y') {
            peak_anno.ion_base_type = "y";
        } else if (annotation[i] == '?') {
            break;
        } else {
            peak_anno.ion_base_type = "others";
        }
        i++;
        int pos = atoi(annotation.substr(i, 1).c_str());
        i++;
        while (annotation[i] >= '0' and annotation[i] <= '9') {
            pos *= 10;
            pos += atoi(annotation.substr(i, 1).c_str());
            i++;
        }

        peak_anno.pos = pos;
        if (annotation[i] == '-' or annotation[i] == '+')// start NL
        {
            int nl = 0;
            int j = 1;
            while (annotation[i + j] >= '0' and annotation[i + j] <= '9') {
                j++;
            }
            peak_anno.ion_NL = annotation.substr(i, j);
            i += j;
        }

        if (annotation[i] == 'i') {
            peak_anno.isotopic = true;
            i++;
        }
        if (annotation[i] == '^') {
            i++;
            peak_anno.charge = atoi(annotation.substr(i, 1).c_str());
            i++;
        }

        if (annotation[i] == '/') {
            i++;
        }
        int errorlen = 1;
        while (i + errorlen < annotation.length() and annotation[i + errorlen] != ',') {
            errorlen++;
        }
        string masserror = annotation.substr(i, errorlen);
        i += errorlen;
        peak_anno.masserror = atof(masserror.c_str());
        pka.push_back(peak_anno);

        i++;
    }

}

void parse_annotation_as_structX(string annotation, vector<peak_annotation> &pkas) {
    char *pEnd = (char *)annotation.c_str();
    if(not annotation.empty()) pkas.emplace_back();

    while(*pEnd!='\0'){
        switch(*pEnd){
            case 'b': pEnd += 1; pkas.back().ion_base_type = "b"; pkas.back().pos = strtod(pEnd, &pEnd); break;
            case 'y': pEnd += 1; pkas.back().ion_base_type = "y"; pkas.back().pos = strtod(pEnd, &pEnd);break;
            case 'c':
            case 'z':
            case 'x':
            case 'm':
            case 'p':
            case 'a': pEnd +=1 ; pkas.back().ion_base_type = "others";pkas.back().pos = strtod(pEnd, &pEnd); break;
            case ':': pEnd +=1; pkas.back().pos = strtod(pEnd, &pEnd); break;
            case 'i': pEnd +=1; pkas.back().isotopic = true; break;
            case '^': pEnd += 1; pkas.back().charge=strtod(pEnd,&pEnd); break;
            case '/': pEnd += 1;  pkas.back().masserror = strtod(pEnd,&pEnd); break;
            case ',': pEnd += 1; pkas.emplace_back(); break; // start a new one
            case '-': pkas.back().ion_NL = to_string(strtol(pEnd,&pEnd,10)); break;
            case '?': pkas.resize(0); pEnd += 1; break;

            default:
                do{
                    pEnd += 1;
                } while(*pEnd >='A' and *pEnd <='Z');
                // other ion type
//                cout << pEnd << " Error " << endl << annotation  << endl;
                pkas.back().ion_base_type = "others";pkas.back().pos = strtod(pEnd, &pEnd);

                break;
        }

    }
}
