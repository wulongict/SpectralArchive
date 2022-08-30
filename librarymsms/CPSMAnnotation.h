//
// Created by wulong on 3/10/21.
//

#ifndef MYTOOL_CPSMANNOTATION_H
#define MYTOOL_CPSMANNOTATION_H

#include "Util.h"
#include <string>
#include <vector>
using namespace std;

class CDBEntry;
struct SPsmAnnotation {
    // new members
    long idx;
    int ms2idx;
    int fileid;

    int isDecoy;
    int significance;
    string protein;

    string peptideseq;// = "UNKNOWN";
    double score;// = -1;
    double cterm_mass;// = 0;
    double nterm_mass;// = 0;
    string modificationstr;// = "UNMODIFIED";
    double precursorMz;// = -1;
    int charge;// = -1;
    int ms2_scan;
    string mzxml_filename;
    double retentiontimeinsec;// = 0;
    double pProb;// = 0;
    double iProb;// = 0;
    double rfscore;//=0
    double precursorNeutralMass() const    {
        double pmass = 1.007276;
        return precursorMz * charge - charge * pmass;
    }
    string m_collision_energy; // collision energy
    int peaknum;
    string m_neighbors;

    SPsmAnnotation();

    friend  ostream &operator<<(ostream & out, const SPsmAnnotation & gt);


    SPsmAnnotation(int scan, double precursor, int chg, double rt);

    void set(int scan, double precursor, int chg, double rt, int fileId, long total_idx, string specfile, int ms2counts);
    void setCollisionEnergy(string collisionEnergy);

    // estimate charge state from mz and mass.
    bool fixChargeStates();

    void initWithRow(CDBEntry &results);

    string createInsertSql() const;
    string createUpdateSql() const;
    string createChargeUpdateSql() const;

    void outputToFile(ofstream &fout) const    {
        fout << idx << "\t" << mzxml_filename << "\t" << fileid << "\t" << ms2idx << "\t" << peptideseq
             << "\t" << score << "\t" << ms2_scan << "\t" << cterm_mass << "\t" << nterm_mass
             << "\t" << modificationstr << "\t" << precursorMz << "\t" << charge << "\t" << retentiontimeinsec
             << "\t" << pProb << "\t" << iProb << "\t" << isDecoy << "\t" << significance << "\t" << protein << endl;
    }

    void toOstringStreamNoId(ostringstream &oss) const;

    void initMembers();

    string getModifiedPeptide(bool verbose) const;
    bool isSig()const {return significance == 1;}

    // vector of (idx, distance) pairs
    struct IdxDistPair{
        struct idxDist{
            double dist;
            long idx;
        };
        vector<idxDist> m_data;

    };


    IdxDistPair  neighborStrToIdxDistPair( char delimitor_1 = ';', char delimitor_2='@') const{
        IdxDistPair idxDistPairs;
        if(m_neighbors!="NULL"){
            vector<string> distId;
            split_string(m_neighbors, distId, delimitor_1);
            for(const auto& eachDistId: distId) {
                vector<string> tmp;
                split_string(eachDistId, tmp, delimitor_2);
                long newId = strtol(tmp[1].c_str(),nullptr, 10);
                double newDist = strtod(tmp[0].c_str(),nullptr);
                idxDistPairs.m_data.push_back(IdxDistPair::idxDist{newDist, newId});
            }
        }
        return idxDistPairs;
    }
};
string getmodificationfrompeptidestring(string peptidestr, modification &mod);
void getmodificationfrompeptidestring(string peptidestr, SPsmAnnotation &gtinfo);

#endif //MYTOOL_CPSMANNOTATION_H
