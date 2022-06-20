//
// Created by wulong on 11/25/17.
//

#include "CPeakPair.h"
#include "../../../librarymsms/Util.h"
#include <fstream>

using namespace std;


CPeakPair::CPeakPair() {
    m_sample_id = -1;
}

CPeakPair::~CPeakPair() {

}


CPeakPair::CPeakPair(const CPeakPair &other) {
    m_PSM_ID = other.m_PSM_ID;
    m_peptide = other.m_peptide;
    m_precursor_charge = other.m_precursor_charge;
    m_intensity_of_peakA = other.m_intensity_of_peakA;
    m_left_aa_of_peakA = other.m_left_aa_of_peakA;
    m_right_aa_of_peakA = other.m_right_aa_of_peakA;
    m_Aby = other.m_Aby;
    m_Apos = other.m_Apos;
    m_intensity_of_peakB = other.m_intensity_of_peakB;
    m_left_aa_of_peakB = other.m_left_aa_of_peakB;
    m_right_aa_of_peakB = other.m_right_aa_of_peakB;
    m_Bby = other.m_Bby;
    m_Bpos = other.m_Bpos;
    m_y = other.m_y;
    m_sample_id = other.m_sample_id;
    m_AAcounts = other.m_AAcounts;
}

CPeakPair &CPeakPair::operator=(const CPeakPair &other) {
    m_PSM_ID = other.m_PSM_ID;
    m_peptide = other.m_peptide;
    m_precursor_charge = other.m_precursor_charge;
    m_intensity_of_peakA = other.m_intensity_of_peakA;
    m_left_aa_of_peakA = other.m_left_aa_of_peakA;
    m_right_aa_of_peakA = other.m_right_aa_of_peakA;
    m_Aby = other.m_Aby;
    m_Apos = other.m_Apos;
    m_intensity_of_peakB = other.m_intensity_of_peakB;
    m_left_aa_of_peakB = other.m_left_aa_of_peakB;
    m_right_aa_of_peakB = other.m_right_aa_of_peakB;
    m_Bby = other.m_Bby;
    m_Bpos = other.m_Bpos;
    m_y = other.m_y;
    m_sample_id = other.m_sample_id;
    m_AAcounts = other.m_AAcounts;
    return *this;
}


int CPeakPair::getSampleID() {
    return m_sample_id;
}

void CPeakPair::setSampleID(int sample_id) {
    if (m_sample_id == -1) {
        m_sample_id = sample_id;
    }
}

double CPeakPair::getFoldChange() {
    if (fabs(m_intensity_of_peakA) < epsilon and fabs(m_intensity_of_peakB) < epsilon )    {
        cout << "Two zero" << endl;
        return 1;
    } else if(fabs(m_intensity_of_peakB) < epsilon)    {
        return DIVIDED_BY_ZERO_INF;
    }
    return m_intensity_of_peakA / m_intensity_of_peakB;
}

void CPeakPair::print() {
    cout << "PSM ID: " << m_PSM_ID << "    Index=" << m_sample_id << endl;
    cout << "Peptide: " << m_peptide << endl;
    cout << m_left_aa_of_peakA << " "
         << m_right_aa_of_peakA << " "
         << m_left_aa_of_peakB << " "
         << m_right_aa_of_peakB << " "
         << m_Apos << " "
         << m_Bpos << " "
         << m_Aby << " "
         << m_Bby << " "
         << m_y << " " << endl;
}

// could be optimized
int CPeakPair::getAAcounts(char aa) {
    if (m_AAcounts.size() == 0) {
        m_AAcounts.assign(26, 0);
        for (int i = 0; i < m_peptide.size(); i++) {
            int pos = m_peptide[i] - 'A';
            m_AAcounts[pos]++;
        }
    }
    int pos = aa - 'A';
    return m_AAcounts[pos];
}

ostream &operator<<(ostream &fout, CPeakPair &peakPair) {
    fout << "\nm_PSM_ID: " << peakPair.m_PSM_ID
         << "\nm_peptide: " << peakPair.m_peptide
         << "\nm_precursor_charge: " << peakPair.m_precursor_charge
         << "\nm_Aintensity: " << peakPair.m_intensity_of_peakA
         << "\nm_Aleft: " << peakPair.m_left_aa_of_peakA
         << "\nm_Aright: " << peakPair.m_right_aa_of_peakA
         << "\nm_Aby: " << peakPair.m_Aby
         << "\nm_Apos: " << peakPair.m_Apos
         << "\nm_Bintensity: " << peakPair.m_intensity_of_peakB
         << "\nm_Bleft: " << peakPair.m_left_aa_of_peakB
         << "\nm_Bright: " << peakPair.m_right_aa_of_peakB
         << "\nm_Bby: " << peakPair.m_Bby
         << "\nm_Bpos: " << peakPair.m_Bpos
         << "\nm_y : " << peakPair.m_y << endl;

    return fout;
}

istream &operator>>(istream &fin, CPeakPair &ts) {
    fin >> ts.m_PSM_ID >> ts.m_peptide >> ts.m_precursor_charge
        >> ts.m_intensity_of_peakA
        >> ts.m_left_aa_of_peakA >> ts.m_right_aa_of_peakA >> ts.m_Aby
        >> ts.m_Apos >> ts.m_intensity_of_peakB >> ts.m_left_aa_of_peakB
        >> ts.m_right_aa_of_peakB >> ts.m_Bby >> ts.m_Bpos
        >> ts.m_y;

    return fin;
}

void getOneHotEncodingAA(map<char, vector<int>> &aaToOneHotEncoding) {
    string AA = "ACDEFGHIKLMNPQRSTVWY";
    for (int i = 0; i < AA.length(); i++) {
        vector<int> x(AA.length(), 0);
        x[i] = 1;
        aaToOneHotEncoding[AA[i]] = x;
    }

}

void export_trainingsample(const vector<CPeakPair> &sample, const string &outfile) {
    ofstream fout;
    string sep = "\t";
    fout.open(outfile.c_str(), ios_base::out);
//    fout.close();

    SimpleTimer st("Export training samples");
    map<char, vector<int>> aa_to_binary;
    getOneHotEncodingAA(aa_to_binary);

    Progress ps(sample.size());

    for (auto x: sample) {
        ps.increase();
        fout << x.m_y << sep;
        fout << (x.m_Aby == 'y') << sep;
        fout << x.m_Apos << sep;
        for (auto y : aa_to_binary[x.m_left_aa_of_peakA]) { fout << y << sep; }
        for (auto y : aa_to_binary[x.m_right_aa_of_peakA]) { fout << y << sep; }

        fout << (x.m_Bby == 'y') << sep;
        fout << x.m_Bpos << sep;
        for (auto y : aa_to_binary[x.m_left_aa_of_peakB]) { fout << y << sep; }
        for (auto y : aa_to_binary[x.m_right_aa_of_peakB]) { fout << y << sep; }
        fout << endl;

    }
    fout.close();
}

COneHotEncodingAA::COneHotEncodingAA() {
    // U is selenocysteine
    //  according to : https://en.wikipedia.org/wiki/Selenocysteine
    m_legalAAs = "ACDEFGHIKLMNPQRSTUVWY";

    m_invalidAA.assign(m_legalAAs.length(), 0);
    for (int i = 0; i < m_legalAAs.length(); i++) {
        vector<int> x(m_legalAAs.length(), 0);
        x[i] = 1;
        aaToOneHotEncoding[m_legalAAs[i]] = x;
    }

}

vector<int>& COneHotEncodingAA::getEncoding(char c, bool verbose) {
    if(aaToOneHotEncoding.count(c)==0){
        if(verbose) cout << "[Error] Invalid amino acid symbol: " << c << endl;
        if(m_invalid_attempts.count(c)==0){
            m_invalid_attempts[c]=0;
        }else{
            m_invalid_attempts[c] ++;
        }
        return m_invalidAA;
    }
    return aaToOneHotEncoding[c];
}

string COneHotEncodingAA::getlegalAAs() {return m_legalAAs;}
