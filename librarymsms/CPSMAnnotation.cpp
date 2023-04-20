//
// Created by wulong on 3/10/21.
//

#include "CPSMAnnotation.h"
#include "DatabaseManager.h"

SPsmAnnotation::SPsmAnnotation() {
    initMembers();
}

void SPsmAnnotation::initMembers() {
    idx = -1;
    ms2idx = -1;
    fileid = -1;

    peptideseq = "UNKNOWN";
    m_rescued_peptide = "";
    score = -1;
    cterm_mass = 0;
    nterm_mass = 0;
    modificationstr = "UNMODIFIED";
    precursorMz = 0;
    charge = 0;  // 10

    ms2_scan = 0;
    mzxml_filename = "UNKNOWN";

    retentiontimeinsec = 0;
    pProb = 0;
    iProb = 0;
    rfscore = 0;

    isDecoy = 0;
    significance = 0;
    protein = "";
    m_collision_energy = ""; // 20

    peaknum = 0;

    m_neighbors = "";

}

void SPsmAnnotation::initWithRow(CDBEntry &results) {
    // todo: the results might contain the same information in different order.
//    cout << "Warning: the results might contain the same information in different order" << endl;
    // ID  FILEID MS2COUNTS  PEPTIDE SCORE  SCAN CTERM  NTERM MODIFICATION  PRECURSOR CHARGE  RT
// PEPTIDEPROPHETPROB  IPROPHETPROB RFSCORE  ISDECOY SIGNIFICANCE  PROTEIN CE  ALTERPEPTIDE NEIGHBOR
    idx = results.getInt("ID",0);
    ms2idx = (int)results.getInt("MS2COUNTS",0);
    fileid = (int)results.getInt("FILEID",0);

    peptideseq = results.get("PEPTIDE",0);
    score = results.getFloat("SCORE",0);
    cterm_mass = results.getFloat("CTERM",0);
    nterm_mass = results.getFloat("NTERM",0);

    modificationstr = results.get("MODIFICATION",0);
    // gtrow.toString(CSqlGtTableRow::MODIFICATION);
    precursorMz = results.getFloat("PRECURSOR", 0);
    // gtrow.toDouble(CSqlGtTableRow::PRECURSOR);
    charge = (int)results.getInt("CHARGE",0);
    // gtrow.toInt(CSqlGtTableRow::CHARGE);
    ms2_scan = (int)results.getInt("SCAN",0);
    // gtrow.toInt(CSqlGtTableRow::SCAN);

    retentiontimeinsec = results.getFloat("RT",0);// gtrow.toDouble(CSqlGtTableRow::RT);
    pProb =results.getFloat("PEPTIDEPROPHETPROB",0);// gtrow.toDouble(CSqlGtTableRow::PEPTIDEPROPHETPROB);
    iProb =results.getFloat("IPROPHETPROB",0);// gtrow.toDouble(CSqlGtTableRow::IPROPHETPROB);
    rfscore = results.getFloat("RFSCORE",0);// gtrow.toDouble(CSqlGtTableRow::RFSCORE);

    isDecoy = (int)results.getInt("ISDECOY",0);// gtrow.toInt(CSqlGtTableRow::ISDECOY);
    significance =  (int)results.getInt("SIGNIFICANCE",0);// gtrow.toInt(CSqlGtTableRow::SIGNIFICANCE);
    protein = results.get("PROTEIN",0);// gtrow.toString(CSqlGtTableRow::PROTEIN);

    m_neighbors = results.get("NEIGHBOR",0);
    m_rescued_peptide = results.get("RESCUEDPEPTIDE",0);


}

string SPsmAnnotation::getModifiedPeptide(bool verbose) const{
    if(nullptr == this ) {
        return "UNKNOWN";
    }
    ostringstream oss;
    // no modification
    if(fabs(nterm_mass) < EPSILON ) // nterm_mass == 0
    {
        // no nterm modification
    } else {
        // there is nterm modification
        oss << "n["<< std::fixed << std::setprecision(0) << nterm_mass << "]";
    }
    map<int,double> pos_mass;
    if(modificationstr != "UNMODIFIED")    {
        vector<string> mod;
        split_string(modificationstr,mod,'|');

        for( const string& mod_i : mod)    {
            if(verbose)cout << mod_i << endl;
            if(mod_i.empty()) continue;
            else{
                vector<string> tokens;
                split_string(mod_i,tokens,'@');
                if(verbose)    cout << tokens[0] << tokens[1] << endl;
                pos_mass[stringTo<int>(tokens[1])]=stringTo<double>(tokens[0]);
            }
        }
    }

    for(int j = 0; j < peptideseq.size(); j ++)
    {
        if(verbose)cout << "dictionary has key: " << pos_mass.count(j+1) << endl;
        oss << peptideseq[j];
        if(pos_mass.count(j+1)==0)        {
            // not modified, do nothing
        }        else        {
            oss << "[" << std::fixed << std::setprecision(0) << pos_mass[j+1] << "]";
        }
    }
    // this is because I have never seen any c-term modification ...
    if(verbose)cout << "cterm modification is not considered" << endl;
    return oss.str();
}

SPsmAnnotation::SPsmAnnotation(int scan, double precursor, int chg, double rt) {
    initMembers();
    set(scan, precursor, chg, rt,-1,-1,"",-1);

}

void SPsmAnnotation::set(int scan, double precursor, int chg, double rt, int fileId, long total_idx,
                         string specfile,
                         int ms2counts) {
    // only update for valid values.
    if(scan >=0)     ms2_scan = scan;
    if(precursor>0) precursorMz = precursor;
    if(chg > 0) charge = chg;
    if(rt >0) retentiontimeinsec = rt;

    if(fileId >=0) fileid = fileId;
    if(total_idx >=0) idx = total_idx;
    if(not specfile.empty())mzxml_filename = specfile;
    if(ms2counts>=0) ms2idx = ms2counts;
}

string SPsmAnnotation::createInsertSql() const {
    ostringstream oss;
    oss << idx << "," << fileid << "," << ms2idx << ",'" << peptideseq
        << "'," << score << "," << ms2_scan << "," << cterm_mass << "," << nterm_mass
        << ",'" << modificationstr << "'," << precursorMz << "," << charge << "," << retentiontimeinsec
        << "," << pProb << "," << iProb << "," << isDecoy << "," << significance << ",'" << protein
        << "','" << m_collision_energy << "'," << rfscore;

    string sql = "INSERT INTO GROUNDTRUTH (ID,FILEID,MS2COUNTS,PEPTIDE,SCORE,SCAN, "
                 "CTERM,NTERM, MODIFICATION,PRECURSOR,CHARGE,RT,PEPTIDEPROPHETPROB, "
                 "IPROPHETPROB, ISDECOY, SIGNIFICANCE, PROTEIN, CE, RFSCORE) "
                 "VALUES (" + oss.str() + "); ";
    return sql;
}

string SPsmAnnotation::createUpdateSql() const {
    ostringstream oss;
    oss << "PEPTIDE='" << peptideseq << "'"
        << ",MODIFICATION='" << modificationstr << "'"
        << ",CHARGE=" << charge
        << ",SCORE=" << score
        << ",CTERM=" << cterm_mass
        << ",NTERM=" << nterm_mass
        << ",PEPTIDEPROPHETPROB=" << pProb
        << ",IPROPHETPROB=" << iProb
        << ",ISDECOY=" << isDecoy
        << ",SIGNIFICANCE=" << significance
        << ",PROTEIN='" << protein << "'"
        << ",CE='" << m_collision_energy << "'"
        << ",RFSCORE=" << rfscore
            ;

    string sql = "update GROUNDTRUTH set " + oss.str() + " where ID=" + to_string(idx) + ";";
    return sql;
}

string SPsmAnnotation::createChargeUpdateSql() const {
    // update charge only
    ostringstream oss;
    oss << "CHARGE=" << charge ;

    string sql = "update GROUNDTRUTH set " + oss.str() + " where ID=" + to_string(idx) + ";";
    return sql;
}

void SPsmAnnotation::toOstringStreamNoId(ostringstream &oss) const{
    oss << R"("peptide": ")" << peptideseq << "\", "
        << "\"rescued_peptide\": \"" << m_rescued_peptide << "\", "
        << " \"filename\": \"" << mzxml_filename << "\", "
        << " \"precursor\": " << precursorMz
        << R"(, "charge": )" << charge
        << R"(, "scan": )" << ms2_scan
        << R"(, "cterm": ")" << cterm_mass
        << R"(", "nterm": ")" << nterm_mass
        << R"(", "othermod": ")" << modificationstr
        << R"(", "rt": )" << retentiontimeinsec
        << R"(, "score": ")" << score
        << R"(", "pProb": )" << pProb
        << R"(, "iProb": )" << iProb
        << R"(, "rfscore": )" << rfscore
        << R"(, "isDecoy": )" << isDecoy
        << R"(, "protein": ")" << protein
        << R"(", "significance": )" << significance;
}

ostream &operator<<(ostream &out, const SPsmAnnotation &gt) {
    out << gt.idx << ","
            << gt.ms2idx << ","
            << gt.fileid << ","
            << gt.isDecoy << ","
            << gt.significance << ","
            << gt.protein << ","
            << gt.peptideseq << ","
            << gt.score << ","
            << gt.cterm_mass << ","
            << gt.nterm_mass << ","
            << gt.modificationstr << ","
            << gt.precursorMz << ","
        <<gt.charge << ","
        <<gt.ms2_scan << ","
        <<gt.mzxml_filename << ","
        <<gt.retentiontimeinsec << ","
        <<gt.pProb << ","
        <<gt.iProb << ","
        << gt.rfscore << ","
        << gt.peaknum << ",";
    return out;
}

void SPsmAnnotation::setCollisionEnergy(string collisionEnergy) {
    m_collision_energy = collisionEnergy;
}

// estimate peptide mass based on average AA mass, ~ 111 Da.
//          e.g. 10 amino acid will be 1110 Da.
// the precursor mass is actually precursor m/z
//  estimated mass / precursor mz  ~  charge.
// todo: change name of variable. precursor mass.
bool SPsmAnnotation::fixChargeStates() {

    bool ret = false;
    if(charge == 0 and peptideseq!="UNKNOWN"){
        int lenPep= peptideseq.length();
        double avarageAAmass = 111;
        double estimatedMass = avarageAAmass * lenPep;
        int estimatedCharge = round(estimatedMass / precursorMz);

        charge = estimatedCharge;
        ret = true;
    }
    return ret;

}


void getmodificationfrompeptidestring(string peptidestr, SPsmAnnotation &gtinfo) {
    gtinfo.nterm_mass = 0;
    gtinfo.cterm_mass = 0;
    map<int, double> pos_mass;
    // remove charge first
    int split_charge = peptidestr.find_first_of('/');
    if (split_charge != string::npos) {
        peptidestr = peptidestr.substr(0, split_charge);
    }
    string AAseq;
    if ('n' == peptidestr[0])     {
        int startpos = peptidestr.find_first_of('[');
        int endpos = peptidestr.find_first_of(']');
        if (startpos != string::npos and endpos != string::npos) {
            gtinfo.nterm_mass = atof(peptidestr.substr(startpos + 1, endpos).c_str());
            peptidestr = peptidestr.substr(endpos + 1); // remove the n[abc] head
        }
    }
    while (peptidestr.find_first_of('[') != string::npos)     {
        int startpos = peptidestr.find_first_of('[');
        int endpos = peptidestr.find_first_of(']');
        if (startpos != string::npos and endpos != string::npos) {
            double modmass = atof(peptidestr.substr(startpos + 1, endpos).c_str());
            AAseq += peptidestr.substr(0, startpos);
            double modpos = AAseq.length();
            pos_mass.insert(make_pair(modpos, modmass));
            peptidestr = peptidestr.substr(endpos + 1); // remove the n[abc] head
        }
    }
    AAseq += peptidestr;
    gtinfo.peptideseq = AAseq;
    gtinfo.modificationstr = CPosMass(pos_mass).toString();
}