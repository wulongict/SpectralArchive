//
// Created by wulong on 10/29/15.
//
#include <spdlog/spdlog.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ctime>
//#include "DatabaseManager.h"
#include "Util.h"
#include "PeakList.h"
using namespace std;

void trim_space_only(std::string &str) {
    str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
    str.erase(str.find_last_not_of(' ') + 1);
}

string getmodificationfrompeptidestring(string peptidestr, modification &mod) {
    mod.nterm_mass = 0;
    mod.cterm_mass = 0;
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
            mod.nterm_mass = atof(peptidestr.substr(startpos + 1, endpos).c_str());
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
            mod.pos_mass.insert(make_pair(modpos, modmass));
            peptidestr = peptidestr.substr(endpos + 1);
        }
    }
    AAseq += peptidestr;
    return AAseq;
}


CTable::CTable(const string& filename, char delimitor, bool has_header, int skipNum) {
    cout << "[Info] Loading table " << filename << endl;
    m_filename = filename;
    m_delim = delimitor;
    m_has_header = has_header;
    string line;
    if (not File::isExist(filename)) {
        throw ios_base::failure("File does not exist");
    }
    ifstream fin(filename.c_str(), ios::in);

    if (fin.is_open()) {
        CountProgress cps(100000, "loading " + filename);
        while (skipNum-- && getline(fin, line)) {
            cout << "skip line " << line << endl;
        }
        // after skip, start reading has_header
        if (m_has_header)         {
            getline(fin, line);
            trim_space_only(line);
            split_string(line, m_column_header, m_delim);
            // cout << endl << "Header: " << endl << line << endl;
        }
        while (getline(fin, line)) {
            cps.increase();
            vector<string> tokens;

            trim_space_only(line);
            split_string(line, tokens, m_delim);
            m_table.push_back(tokens);
        }
        fin.close();
    }
    // cout << "[Info] " << m_table.size() << " lines loaded." << endl;
    m_row = m_table.size();
    m_col = 0;
    if (m_row > 0) { m_col = m_table[0].size(); }
    // verify table on each row.
    for (int i = 0; i < m_row; i++) {
        vector<string> x = m_table[i];
        if (x.size() > m_col) {
            cout << "Warning: row i=" << i << " columns = " << x.size() << " greater than  " << m_col << endl;
        } else if (x.size() < m_col) {
            cout << "Warning: row i=" << i << " columns = " << x.size() << " less than  " << m_col << endl;
        }
    }

    if (m_column_header.size() != m_col && m_has_header) {
        cout << "Invalid has_header size: has_header size=" << m_column_header.size() << " and column number = "
             << m_col << endl;
        throw logic_error("invalid has_header");
    }

    if (not m_has_header) {
        cout << "Using place-holder header: ";
        for (int i = 0; i < m_col; i++) {
            m_column_header.push_back(string("C") + to_string(i));
            cout << "C" << i << ",";
        }
        cout << endl;
    }
    // cout << "[Info] Table loaded " << filename << endl << endl;
}

string CTable::getEntry(int row, int col) const {
    if (row >= 0 && col >= 0 && m_table.size() > row && m_table[row].size() > col) {
        return m_table[row][col];
    } else {
        return CTable_unavailableEntry;
    }
}

void CTable::setEntry(int row, int col, string value) {
    if (row >= 0 && row < m_row && col >= 0 && col < m_col) {
        m_table[row][col] = std::move(value);
    } else {
        cout << "Error: invalid row col pair (" << row << ", " << col << ")" << endl;
        throw logic_error("invlaid row or column value");
    }
}

// this function is not good, because it changed column number.
// Warning: after using this function, the column row is no longer aligned with other rows. as it has one more element.
void CTable::appendEntry(int row, const string& value) { // not a good one; the header is left out
    if (row >= 0 && row < m_row) {
        m_table[row].push_back(value);
        if (m_col != m_table[row].size()) {
            m_col = m_table[row].size();
        }

        if (m_col > m_column_header.size()) {
            m_column_header.push_back(string("C") + to_string(m_col - 1));
        }
    } else {
        cout << "Error: invalid row number " << row << endl;
        throw logic_error("CTable row number out of range");
    }
}

// column concatenate 
// The header2column might broken after this step, if there are two column with same name
void CTable::Join(CTable &other) {
    if (m_row != other.m_row) {
        cout << "Warning: The row number should be equal! " << endl;
        throw logic_error("Table with different row number can not be merged into one!");
    }
    spdlog::get("A")->info("Start to merge table ");
    int prev_header_size = m_column_header.size();

    for (int j = 0; j < other.m_col; j++) {
        for (int i = 0; i < m_row; i++) {
            appendEntry(i, other.getEntry(i, j));
        }
    }
    buildheader2column(true);

    if (other.hasHeader()) {
        for (int j = prev_header_size; j < m_col; j++) {
            m_column_header[j] = other.m_column_header[j - prev_header_size];
        }
    }
}

void CTable::saveAs(const string& filename, bool with_header, char delimtor) {
    ofstream fout(filename.c_str(), ios::out);
    if (fout.is_open()) {
        Progress ps(m_row, "export table " + filename);
        if (with_header) {
//            cout << "Exporting header.." << endl;
            for (int j = 0; j < m_col; j++) {
                fout << m_column_header[j];
                if (j != m_col - 1) fout << delimtor;
            }
            fout << endl;
//            cout << "Header done!" << endl;
        }
        for (int i = 0; i < m_row; i++) {
            ps.increase();
            for (int j = 0; j < m_col; j++) {
                fout << getEntry(i, j);
                if (j != m_col - 1) fout << delimtor;
            }
            fout << endl;
        }
        fout.close();
    } else {
        cout << "Error: Fail to open file " << filename << endl;
    }
}

void CTable::printRow(int i) {
    if (i < m_row && i >= 0) {
        for (const auto& x: m_table[i]) cout << x << ",";
        cout << endl;
    } else{
        cout << "[Error] invalid row num " << i << " valid range: 0 <= rownum <= " << m_row-1 << endl;
    }
}

void CTable::resizeRow(int k, int newsize) {
    if(m_table[k].size() > newsize)    {
        m_table[k].resize(newsize);
    } else  {
        cout << "new size is larger: " << newsize << " > " << m_table[k].size() << "! do nothing!" << endl;
    }
}

void CTable::build_table_index(int col) {
    if (not m_tableindex.empty())    {
        m_tableindex.clear();
        cout << "Table index exist!" << endl;
    }
    for(int i = 0; i < m_row; i ++)    {
        string key = m_table[i][col];
        if(m_tableindex.find(key)==m_tableindex.end())        {
            m_tableindex[key]=vector<int>({i});
        } else {
            m_tableindex[key].push_back(i);
        }
    }
}

void CTable::addRow(const vector<string>& row) {
    m_table.push_back(row);
    m_row ++;
    if (m_col == 0)    {
        m_col = row.size();
        for(int i = 0; i < m_col; i ++)        {
            m_column_header.push_back(string("C")+to_string(i)); // fake col number
        }
    } else if (m_col != row.size()) {
        throw logic_error("can not add row with different column size!");
    }
}

void CTable::setHeader(const vector<string> &header) {
    m_has_header = true;
    if (m_col == 0 or m_col == header.size())   {
        m_col = header.size();
        m_column_header = header;
    } else {
        throw logic_error("input header size is less than the column number!");
    }
    buildheader2column(true);
}

int CTable::getRowByKey(const string& key, int col) {
    if(m_tableindex.empty())   {
        build_table_index(col);
    }
    if(m_tableindex.find(key) == m_tableindex.end())    {
        return -1;
    }  else  {
        if(m_tableindex[key].size()>1){
            cout << "[Warning] multiple row found for key " << key << endl;
        }
        return m_tableindex[key][0];
    }
}

void CTable::print() {
    for(size_t i = 0; i < m_row; i ++)
    {
        printRow(i);
    }
}

CTable::CTable() {
    m_row = 0;
    m_col = 0;
    m_filename = "";
    m_delim = '\t';
    m_has_header = false;

}

void CTable::appendHeader(const vector<string> &Header) {
    for(auto & i : Header)     {
        m_column_header.push_back(i);
        m_col += 1;
    }
    buildheader2column(true);

}

string CTable::getHeaderByColumn(size_t column) {
    if(column >=m_col)        {
        throw logic_error("column number out of header size");
    }
    return m_column_header[column];
}

bool CTable::hasColumn(const string& colname)
{
    buildheader2column();
    return m_header2col.count(colname)>0;
}

int CTable::getColByHeader(const string &header) {
    buildheader2column();
    if(m_header2col.count(header)==0)
    {
        throw logic_error("header not found: " + header );
    }
    return m_header2col[header];
}

void CTable::buildheader2column(bool rebuild) {
    if(m_header2col.empty() or rebuild)
    {
        for(size_t i = 0; i < m_col; i ++)
        {
            m_header2col[m_column_header[i]]=i;
        }
        // cout << "column index build" << endl;
    }
}

unsigned long CTable::LengthofRow(size_t k) {
    return m_table[k].size();
}

bool CTable::hasHeader() const {
    return m_has_header;
}

namespace statistic {
    double calcGeneralizedMean(const vector<double>& v, double p) {
        double genrealizedMean = 0;
        for (double i : v) {
            genrealizedMean += pow(i, p);
        }
        genrealizedMean /= v.size();
        genrealizedMean = pow(genrealizedMean, 1 / p);
        return genrealizedMean;
    }

    double calcGeneralizdMean(vector<double> &v, vector<int> &nonzeros, double p) {
        double generalizedMean = 0;
        for (int nonzero : nonzeros) {
            generalizedMean += pow(fabs(v[nonzero]), p);
        }
        generalizedMean /= v.size();
        generalizedMean = pow(generalizedMean, 1 / p);
        return generalizedMean;
    }

    double getNonzeroMin(const vector<double> &v) {
        double nonzeromin = DBL_MAX;
        for (double i : v) {
            if (i > 0 && i < nonzeromin) nonzeromin = i;
        }
        return nonzeromin;
    }

    void checkAndFixVector(vector<double> &v, double &min, double &max) {
        if (min < 0) {
            cout << "Vector should be non-negative real number" << endl;
            cout << min << endl;
            throw "[Error] vector for calculation Harmonic Mean and Geometric Mean should be non negative real numbers.";
        } else if (min == 0 && min != max) {
            double nonzeromin = getNonzeroMin(v);
            double fakeValue = nonzeromin / v.size();
            for (double & i : v) {
                if (i == 0) i = fakeValue;
            }
        }
    }

    double calcHarmonicMean(vector<double> v) {
        // if there are some zeros, we should take care of Devide by zero errordouble min;
        double max, min;
        min = *(std::min_element(v.begin(), v.end()));
        max = *std::max_element(v.begin(), v.end());
        checkAndFixVector(v, min, max);

        double ret = 0;
        if (min == 0 && min == max) { ret = 0; }
        else if (min < 0) {
            cout << "Error" << endl;
            ret = -1;
        } else { ret = calcGeneralizedMean(v, -1); }
        return ret;
    }

    double calcGeometricMean(vector<double> v) {
        double max, min;
        min = *std::min_element(v.begin(), v.end());
        max = *std::max_element(v.begin(), v.end());
        checkAndFixVector(v, min, max);
        double ret = 0;
        for (double i : v) {
            ret += log(i);
        }
        ret /= v.size();
        return exp(ret);
    }

    double calcArithmeticMean(const vector<double>& v) {
        return calcGeneralizedMean(v, 1);
    }

    double calcQuadraticMean(const vector<double>& v) {
        return calcGeneralizedMean(v, 2);
    }

    double calcCubicMean(const vector<double>& v) {
        return calcGeneralizedMean(v, 3);
    }
}


std::vector<double> AnalysisMassDiff::CalculateMassDifference(std::vector<PeakList *> &vpl, double tolerance) {
    std::cout << "[Info] Calculate mass difference with tolerance " << tolerance << "Th" << endl;
    int Count = 0;
    int MaxDataPointsNum = 10000;
    vector<double> massdiff;
    for (int i = 0; i < vpl.size() - 1; ++i) {
        PeakList *currentpl = vpl[i];
        cout << i << " / " << vpl.size() << "\r" << flush;
        if (Count > MaxDataPointsNum) break;
        for (int j = 0; j < currentpl->getM_mzList().size(); ++j) {
            double mz = currentpl->getM_mzList()[j];
            PeakList tmp = vpl[i + 1]->getPeaksWithin(mz - tolerance, mz + tolerance);
            if (!tmp.getM_mzList().empty()) {
                massdiff.push_back(minDiff(mz, tmp.getM_mzList()));
                Count++;
            }
        }
    }
    cout << endl;
    double AvgDiff = statistic::calcmean(massdiff);
    double stdDiff = statistic::calcstd(massdiff);
    cout << "[Resu] Mass deviation:" << endl;
    cout << "[Resu] mean " << AvgDiff << endl;
    cout << "[Resu] std " << stdDiff << endl;
    return massdiff;
}

double AnalysisMassDiff::minDiff(double mz, const vector<double>& mzlist) {
    double minMassDiff = 100;
    for (double i : mzlist) {
        double massdiff = i - mz;
        if (fabs(minMassDiff) > fabs(massdiff)) {
            minMassDiff = massdiff;
        }
    }
    return minMassDiff;
}

long CMyMatrix::getM_col() const {
    return m_col;
}

long CMyMatrix::getM_row() const {
    return m_row;
}

CMyMatrix::CMyMatrix(const CMyMatrix &other) {
    m_row = other.getM_row();
    m_col = other.getM_col();
    m_Size = m_row * m_col;
    m_entries = new double[m_Size];
    if (m_entries == nullptr) {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i, j, other.get(i, j));
        }
    }
}

CMyMatrix::CMyMatrix(long r, long c) {
    m_row = r;
    m_col = c;
    m_Size = r * c;
    m_entries = new double[m_Size];
    if (m_entries == nullptr) {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    initialize();
}

CMyMatrix::CMyMatrix(long r, long c, const double *matrix) {
    m_row = r;
    m_col = c;
    m_Size = r * c;
    m_entries = new double[m_Size];
    if (m_entries == nullptr) {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i, j, matrix[i * m_col + j]);
        }
    }
}

CMyMatrix::~CMyMatrix() {

        delete[] m_entries;
}

void CMyMatrix::Print() const {
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            cout << get(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void CMyMatrix::mergeMatrix(vector<CMyMatrix> &vMatrix, const string& mergeMethod) {
    cout << "[Info] Calculate sum of Matrix" << endl;
    initialize();// This is just add all of them directly together, without any tune of the weight.
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            double val = 0;
            vector<double> v;
            for (auto & k : vMatrix) {
                v.push_back(k.get(i, j));
            }
            //----------------------a bunch of methods-------------------
            if (mergeMethod == "Geometric") {
                val = statistic::calcGeometricMean(v);
            } else if (mergeMethod == "Arithmetic") {
                val = statistic::calcArithmeticMean(v);
            } else if (mergeMethod == "Harmonic") {
                val = statistic::calcHarmonicMean(v);
            } else {
                cout << "[Error] incorrect merge method: " << mergeMethod << endl;
                throw "Invald merge method.";
            }
            set(i, j, val);
        }
    }
}

void CMyMatrix::HomonicMatrix(vector<CMyMatrix> &vMatrix) {
    cout << "[Info] Calculate sum of Matrix" << endl;
    initialize();// This is just add all of them directly together, without any tune of the weight.
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            double val = 0;
            for (auto & k : vMatrix) {
                double sim = k.get(i, j);
                if (fabs(sim) <= EPSILON) {
                    sim = 1e-3;
                }
                val += 1.0 / sim;
            }
            if (fabs(val) <= EPSILON)
                throw "Devide by Zero";
            set(i, j, 1.0 / val);
        }
    }
}

void CMyMatrix::gmOfMatrix(vector<CMyMatrix> &vMatrix) {
    cout << "[Info] Calculate point-wise product of Matrix" << endl;
    initialize();// This is just add all of them directly together, without any tune of the weight.
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            for (auto & k : vMatrix) {
                double val = get(i, j) * (0.1 + k.get(i, j));
                set(i, j, val);
            }
            double gm = pow(get(i, j), 1 / vMatrix.size());
            set(i, j, gm);
        }
    }
}

void CMyMatrix::outputAsText(const string& outputfilename) {
    cout << "[Info] Export a matrix as text file: " << outputfilename << endl;
    FILE *pfile = fopen(outputfilename.c_str(), "w");
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            fprintf(pfile, "%.4lf\t", m_entries[i * m_col + j]);
        }
        fprintf(pfile, "\n");

    }
    fclose(pfile);
}

void CMyMatrix::outputBinary(const string& outputfilename) {
    cout << "[Info] Export a matrix as binary file:" << outputfilename << endl;
    FILE *pfile = fopen(outputfilename.c_str(), "wb");
    fwrite(&m_row, sizeof(long), 1, pfile);
    fwrite(&m_col, sizeof(long), 1, pfile);
    fwrite(m_entries, sizeof(double), m_Size, pfile);
    fclose(pfile);
}

CMyMatrix::CMyMatrix(string filename) {
    FILE *pfile = fopen(filename.c_str(), "rb");
    if (pfile == nullptr) {
        cout << "can not open file" << endl;
        throw "Fail to open file";
    }

    int len = fread(&m_row, sizeof(long), 1, pfile);
    len = fread(&m_col, sizeof(long), 1, pfile);
    m_Size = m_row * m_col;
    m_entries = new double[m_Size];
    len = fread(m_entries, sizeof(double), m_Size, pfile);
    fclose(pfile);
}

double CMyMatrix::get(long i, long j) const {
    return m_entries[i * m_col + j];
}

void CMyMatrix::set(long i, long j, double value) {
    m_entries[i * m_col + j] = value;
}

void CMyMatrix::initialize(double val, long r, long c) {
    if (r == -1 || c == -1) {
        r = m_row;
        c = m_col;
    }
    if (m_row != r || m_col != c) {
        cout << "[Info] matrix resize to " << r << "-by-" << c << endl;
        if (m_entries != nullptr) {
            delete[] m_entries;
            m_entries = nullptr;
        }
        m_row = r;
        m_col = c;
        m_Size = r * c;
        m_entries = new double[m_Size];
        if (m_entries == nullptr) {
            cout << "[Info] Fail to create matrix" << endl;
            throw "Fail to create matrix";
        }
    }
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i, j, val);
        }
    }
}


void TestMatrix() {
    try {
        CMyMatrix m(2, 2);
        m.set(1, 1, 2.0);
        m.set(0, 0, 1.0);
        string filename = "testMatrix.binary";
        m.outputBinary(filename);
        m.Print();
        CMyMatrix n(filename + "x");
        n.Print();
        cout << "[True] Testing Matrix" << endl;
    }
    catch (char const *s) {
        cout << "[Error] " << s << endl;
        cout << "[False] Testing Matrix" << endl;
    }
}




SimpleTimer::SimpleTimer(bool verbose) {
    m_start = std::chrono::steady_clock::now();//clock();
    m_end = m_start;
    m_used = m_end - m_start;
    m_taskname = "";
    m_verbose = verbose;
//    cout << m_verbose << " timer verbosity: " << endl;
    if(m_verbose) {
        cout << "[Timer] " << flush;
        cout << "Timer is started..." << endl;
    }
}

SimpleTimer::SimpleTimer(string taskname, bool verbose) {
    m_start = std::chrono::steady_clock::now();
    m_end = m_start;
    m_used = m_end - m_start;
    m_taskname = taskname;
    m_verbose = verbose;
//    cout << m_verbose << " timer verbosity: " << taskname << "\t" << m_taskname << endl;
    if(m_verbose){
        cout << "[Timer] " << flush;
        cout << m_taskname << " Timer is started..." << endl;
    }

}
SimpleTimer::SimpleTimer(const char *x, bool verbose): SimpleTimer(string(x), verbose) {}

SimpleTimer::~SimpleTimer() {
    stop();
}

double SimpleTimer::secondsElapsed() {
    m_used = std::chrono::steady_clock::now() - m_start;
    return m_used.count();
}

double SimpleTimer::stop() {
    m_end = std::chrono::steady_clock::now();
    m_used = m_end - m_start;
    if(m_verbose){
        cout << "[Timer] " << flush;
        std::cout <<  m_taskname << " Time used: "
                  << std::fixed << setprecision(4)
                  << m_used.count() << "s." << endl;
    }

    return m_used.count();
}

double SimpleTimer::restart(const string& taskname) {
    double x = stop();
    // initialization of new clock
    m_start = std::chrono::steady_clock::now();
    m_end = m_start;
    m_used = m_end - m_start;
    if (!taskname.empty()) {
        m_taskname = taskname;
    }

    if(m_verbose){
        cout << "[Timer] " << flush;
        cout << m_taskname << " Timer is started..." << endl;
    }


    return x;
}




int getProperThreads(int threadnum) {
    int threads = std::thread::hardware_concurrency() - 1;
    if (threads == 0) {
        threads += 1;
    }
    if(threadnum > 0 and threadnum < threads){
        threads = threadnum;
    }
    return threads;
}

void get_FDR_CorrectNum(vector<double> &tProbs, vector<double> &dProbs, vector<tuple<double, double>> &FDR_CorrectNum) {
    sort(tProbs.begin(), tProbs.end(), [](const double &x, const double &y) -> bool { return x > y; });
    sort(dProbs.begin(), dProbs.end(), [](const double &x, const double &y) -> bool { return x > y; });

    int tcount = 0, dcount = 0;
    FDR_CorrectNum.push_back(make_tuple(0.0, 0.0));
    while (tcount < tProbs.size() and dcount < dProbs.size()) {
        if (tProbs[tcount] >= dProbs[dcount]) {
            tcount++;
        } else  {
            dcount++;
            // todo: decoy is defined as
            //  FDR = D/T x 100%, and we should use pi_0
            double FDR = dcount * 1.0 / tcount;
            FDR_CorrectNum.emplace_back(FDR, tcount * 1.0);
        }

    }
}

CMyMatrix::CMyMatrix() {
    m_row = 0;
    m_col = 0;
    m_Size = 0;
    m_entries = nullptr;
}

#include <mutex>

std::mutex progress_mtx;

void Progress::increase(int n ) {
    lock_guard<mutex> guard(progress_mtx);

    for(int k = 0; k < n; k ++)
    {
        m_task_num_finished++;
        while (100.0 * m_task_num_finished / m_task_num >= m_percentage + 1) {
            cout << "." << flush;
            m_percentage += 1;
            if (m_percentage % m_print_percentage_gap == 0)// print value very 10%
            {
                cout << m_percentage << "%:(" << m_task_num_finished << "/" << m_task_num << ")"<< flush;
            }
            if (m_percentage == 100)
                cout << endl;
        }
    }



}


Progress::Progress(long task_num) : m_print_percentage_gap(10) {
    m_task_num = task_num;
    m_task_num_finished = 0;
    m_percentage = 0;
    cout << "Progress:" << flush;

}

Progress::Progress(long task_num, const string& task_name) : m_print_percentage_gap(10) {
    m_task_num = task_num;
    m_task_num_finished = 0;
    m_percentage = 0;
    cout << "[" << task_name << "] Progress:" << flush;
}

void File::splitpath(const string& inputpath, string &path, string &file) {
    int found = inputpath.find_last_of('/');
    if (found != string::npos) {
        path = inputpath.substr(0, found);
        file = inputpath.substr(found + 1);
    } else {
        path = ".";
        file = inputpath;
    }

}

void File::parent(string inputpath, string &path, string &folder) {
    int found = inputpath.find_last_of('/');
    if (found == inputpath.length() - 1) {
        inputpath = inputpath.substr(0, inputpath.length() - 1);
        parent(inputpath, path, folder);

    } else if (found != string::npos) {
        path = inputpath.substr(0, found);
        folder = inputpath.substr(found + 1);
    } else {
        path = ".";
        folder = inputpath;
    }

}

void File::splitname(const string& filename, string &filenameprefix, string &ext) {
    int found = filename.find_last_of('.');

    if (found != string::npos) {
        filenameprefix = filename.substr(0, found);
        ext = filename.substr(found + 1);
    } else {
        filenameprefix = filename;
        ext = "";
    }
}

bool File::isExist(const string &curOutputfile,bool beQuiet) {
    bool ret = false;
    FILE *pfile = fopen(curOutputfile.c_str(), "r");
    if (pfile != nullptr) {
        fclose(pfile);
        ret = true;
    } else {
        if(not beQuiet)  cout << "[Error] File does not exist: \"" << curOutputfile << "\"" << endl;
    }
    return ret;
}

void File::getfilesize(const string &filename, long &filesize) {
    FILE *pfile = fopen(filename.c_str(), "r");
    if (pfile == nullptr) {
        printf("File %s does not exist!", filename.c_str());
        throw runtime_error("Fail to open file");
    }
    fseek(pfile, 0L, SEEK_END);
    filesize = ftell(pfile);
    rewind(pfile);
    fclose(pfile);
}

char *File::loadFileToBuffer(const string &filename) {
    long filesize =0;
    getfilesize(filename, filesize);
    FILE *pfile = fopen(filename.c_str(), "rb");
    char *buf = new char[filesize +1];
    if(buf == nullptr){
        cout << "Error: fail to allocate momery"  << endl;
    }
    long nSize = fread(buf, sizeof(char),filesize,pfile);
    buf[filesize]='\0';
    fclose(pfile);
    return buf;
}

// only support splitting with space
void split_string(string &line, vector<string> &tokens) {
    istringstream iss(line);

    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter(tokens));
//    for(auto x:tokens) cout << x << endl;
}

// support user defined delimiters.
void split_string(const string &line, vector<string> &tokens, char delim) {
    string token;
    istringstream iss(line);
    while (getline(iss, token, delim)) {
        tokens.push_back(token);
    }

}

string make_timestring() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << put_time(&tm, "%d-%m-%Y %H:%M:%S");
    auto str = oss.str();
    return str;

}

string readlinesfromfile(const string& filename) {
    cout << "loading readlines from file: " << filename << endl;
    std::ifstream t(filename);
    if(not t.good()){
        // file not exist!
        throw runtime_error("Fail to open file: " + filename);
    }
    std::string str;

    t.seekg(0, std::ios::end);
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(t)),
               std::istreambuf_iterator<char>());
    t.close();
    return str;
}

bool find_value(vector<long> &a, long val, int &idx) {
    idx = -1;
    bool found = false;
    auto it = std::find(a.begin(), a.end(),val);
    if(it != a.end())
    {
        found = true;
        idx = std::distance(a.begin(), it);
    }
    return found;
}

string getTimeStemp() {
    auto timeNow = std::chrono::system_clock::now();
//    // Some computation here
//    auto end = std::chrono::system_clock::now();

//    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t t = std::chrono::system_clock::to_time_t(timeNow);
    return string(std::ctime(&t));
//    std::cout << "finished computation at " << std::ctime(&t)
//              << "elapsed time: " << elapsed_seconds.count() << "s"
//              << std::endl;
}


CPath::CPath(string fullpath) {
    m_fullpath = std::move(fullpath);
    parse();
}

void CPath::print() {
    cout << "fullpath: " << m_fullpath << endl;
    cout << "folder: " << m_folder << endl;
    cout << "filename: " << m_filename << endl;
    cout << "filename no ext: " << m_filename_no_ext << endl;
    cout << "ext" << m_ext << endl;
}

CFolder CPath::getParentFolder() const {
    return CFolder(m_folder);
}

void CPath::parse() {
    File::splitpath(m_fullpath, m_folder, m_filename);
    File::splitname(m_filename, m_filename_no_ext, m_ext);
}

CFolder::CFolder() {
    m_folder = "";
}

CFolder::CFolder(string folder) {
    m_folder=std::move(folder);
}

CFolder::CFolder(const CFolder &cf) {
    m_folder = cf.m_folder;
}

CFolder& CFolder::operator=(const CFolder &cf) {
    m_folder = cf.m_folder;
    return *this;
}

CFolder CFolder::getParent() {
    string path, foldername;
    File::parent(m_folder,path,foldername);
    return CFolder(path);
}

string CFolder::getFolderName() {
    string path, foldername;
    File::parent(m_folder,path,foldername);
    return foldername;

}

void CFolder::print() {
    cout << "folder: " << m_folder << endl;
}

string CFolder::toStr() {
    return m_folder;
}






void CSqlGtTableRow::toOstringStream(ostringstream &oss) {
    oss << R"("id": ")" << m_results[IDX]
        << R"(", "peptide": ")" << m_results[PEPTIDE]
        << R"(", "precursor": )" << m_results[PRECURSOR]
        << R"(, "charge": )" << m_results[CHARGE]
        << R"(, "cterm": ")" << m_results[CTERM]
        << R"(", "nterm": ")" << m_results[NTERM]
        << R"(", "othermod": ")" <<m_results[MODIFICATION]
        << R"(", "rt": )" << m_results[RT]
        << R"(, "score": ")" << m_results[SCORE]
        << R"(", "pProb": )" << m_results[PEPTIDEPROPHETPROB]
        << R"(, "iProb": )" << m_results[IPROPHETPROB]
        << R"(, "isDecoy": )" << m_results[ISDECOY]
        << R"(, "protein": )" << m_results[PROTEIN]
        << R"(, "significance": )" << m_results[SIGNIFICANCE];
}

void CSqlGtTableRow::toOstringStreamNoId(ostringstream &oss) {
    oss << R"("peptide": ")" << m_results[PEPTIDE]
        << R"(", "precursor": )" << m_results[PRECURSOR]
        << R"(, "charge": )" << m_results[CHARGE]
        << R"(, "cterm": ")" << m_results[CTERM]
        << R"(", "nterm": ")" << m_results[NTERM]
        << R"(", "othermod": ")" <<m_results[MODIFICATION]
        << R"(", "rt": )" << m_results[RT]
        << R"(, "score": ")" << m_results[SCORE]
        << R"(", "pProb": )" << m_results[PEPTIDEPROPHETPROB]
        << R"(, "iProb": )" << m_results[IPROPHETPROB]
        << R"(, "isDecoy": )" << m_results[ISDECOY]
        << R"(, "protein": ")" << m_results[PROTEIN]
        << R"(", "significance": )" << m_results[SIGNIFICANCE];
}

string CSqlGtTableRow::getJsonNode(bool getfilename) {
    ostringstream oss;
    oss << "{"
        << R"("id": ")" << m_results[IDX]
        << R"(", "peptide": ")" << m_results[PEPTIDE]
            << R"(", "filename": ")" << m_results[FILENAME_APPENDED]
            << R"(", "scan": ")" << m_results[SCAN]
        << R"(", "precursor": )" << m_results[PRECURSOR]
        << R"(, "charge": )" << m_results[CHARGE]
        << R"(, "cterm": ")" << m_results[CTERM]
        << R"(", "nterm": ")" << m_results[NTERM]
        << R"(", "othermod": ")" <<m_results[MODIFICATION]
        << R"(", "rt": )" << m_results[RT]
        << R"(, "score": ")" << m_results[SCORE]
        << R"(", "pProb": )" << m_results[PEPTIDEPROPHETPROB]
        << R"(, "iProb": )" << m_results[IPROPHETPROB]
        << R"(, "isDecoy": )" << m_results[ISDECOY]
        << R"(, "protein": ")" << m_results[PROTEIN]
        << R"(", "significance": )" << m_results[SIGNIFICANCE]
        << "}" << endl;
    string jsonnode = oss.str();
    return jsonnode;
}

CHistogram::CHistogram(int N) {
    m_hist.assign(N,0);
}

void CHistogram::add_data(int d) {
    if(d<m_hist.size() and d >=0)
    {
        m_hist[d] ++;
    } else{
        m_hist.back() ++;
    }
}

void CHistogram::display() {
    cout << "-------------Histogram of Ranking ----------------\n//Ranking 0 is not found, \n//Ranking 1 is correct, \n"
            "//>1 scoring funtion to be improved, reranking top hits:\n"
            "--------------- BEGIN----------------------------" << endl;
    int sum = accumulate(m_hist.begin(),m_hist.end(),0);
    cout << "//Ranking\tCounts\tFrequency" << endl;
    for(int i = 0; i < m_hist.size(); i ++)
    {
        if(m_hist[i]>0)
        {
            cout << i << "\t" << m_hist[i] << "\t" << m_hist[i]*1.0/sum << endl;
        }
    }
    cout << "--------------Histogram of Ranking END---------------------" << endl;
}

void CKeyValuesParser::findValue(string &nextstr, string &value, const string& pairseparator) {
    if(nextstr[0]=='\"')        {
        int qmpos=nextstr.substr(1).find_first_of( '\"');
        if(qmpos!=string::npos)            {
            value = nextstr.substr(1,qmpos);
            nextstr=nextstr.substr(qmpos+2);
        }
    }
    else        {
        int npos = nextstr.find_first_of(pairseparator);
        if(npos!=string::npos)       {
            value = nextstr.substr(0,npos);
            nextstr=nextstr.substr(npos+1);
        }     else     {
            value = nextstr;
            nextstr = "";
        }
    }
}

CKeyValuesParser::CKeyValuesParser(string &keyvaluestr, const string& keyvalueseparator, const string& pairseparator, bool verbosity) {
//    m_vmsg = make_shared<CVerboseMessage>();
    if(verbosity)
    {
        m_vmsg=make_shared<CVerboseMessage>();

    } else{
        m_vmsg = make_shared<CDummyMessageCollector>();
    }
//    cout << "start parsing" << endl;
    string key, value;
    string localstr = keyvaluestr;
    int spos = localstr.find_first_not_of(pairseparator);
    int npos= localstr.find_first_of(keyvalueseparator);
    while(npos!=string::npos and spos!=string::npos)   {
        key = localstr.substr(spos,npos-spos);
        localstr=localstr.substr(npos+1);
        findValue(localstr, value,pairseparator);
        setKeyValue(key, value);

        spos= localstr.find_first_not_of(pairseparator);
        npos= localstr.find_first_of(keyvalueseparator);
    }
}

void CKeyValuesParser::setKeyValue(string &key, string &value) {
    m_data[key]=value;
//    m_vmsg->add("setting"+key+":"+value + "\n");
//    cout << key <<": " << value << endl;
}

void CKeyValuesParser::display() {
    cout << "key: value " << endl;
    for(auto & it : m_data)
    {
        cout <<  it.first  << ": " << it.second << endl;
    }
}

string CANSIConsole::getColorStr(const string& str, CANSIConsole::color x) {
    return m_colors_inside[x].toStr()+str+"\033[0m";
}

void CANSIConsole::reset() {
    if(active)
    {
        cout << "\033[0m" ;
        active = false;

    }

}

void CANSIConsole::set(CANSIConsole::color x) {
    ansi_colors &ac = m_colors_inside[x];
    if(not active)
    {
        cout << ac.toStr() ;
        active = true;

    }
}

CANSIConsole::~CANSIConsole() {
    if(active) reset();
}

CANSIConsole::CANSIConsole() {
    active = false;
    m_colors_inside.emplace_back(206,57); // default
    m_colors_inside.emplace_back(1,0); // RED
    m_colors_inside.emplace_back(2,0); // GREEN
    m_colors_inside.emplace_back(4,0); // BLUE
    m_colors_inside.emplace_back(3,0); // YELLOW

}

File::CFile::CFile(const string& inputPath) {
    m_fullpathname = inputPath;
    splitpath(inputPath, path, filename);
    splitname(filename, basename, ext);
}

bool File::CFile::isFileExist(bool beQuiet) const {
    return isExist(m_fullpathname, beQuiet);
}

double CountFrequency::getEntropy(bool normalized, bool verbosity) {
    double sum = 0;
    for(auto &x: m_counts){
        sum += x.second;
        if(verbosity) cout << "sum = " << x.second << endl;
    }
    double entropy = 0;
    if(sum >0){
        for(auto &x: m_counts){
            double p = x.second/sum;
            entropy -= p * log(p);
            if(verbosity)cout << "p: " << p << "\t" << -p*log(p) << " to " << entropy << endl;
        }
        if(normalized and m_counts.size()>1){
            entropy /= log(1.0*m_counts.size());
            if(verbosity) cout << "normalization: " << entropy << " by "<< log(1.0*m_counts.size()) << endl;
        }
    }
    if(verbosity) cout << "entropy" << endl;
    return entropy;
}

void CountFrequency::add_data(int d) {
    m_sum += d;
    if(m_counts.find(d)!=m_counts.end()){
        m_counts[d]++;
    }else{
        m_counts[d] = 1;
    }
}

int CountFrequency::getFreq(int key) {
    if(m_counts.find(key)!=m_counts.end())
        return m_counts[key];
    else return 0;
}

CountFrequency::CountFrequency() {
    m_sum = 0;
}

void CountFrequency::print() {
    cout << "cluster:" ;
    for(auto &item: m_counts){
        cout << "\t" << item.first ;
    }
    cout << endl << "counts";
    for(auto &item: m_counts){
        cout << "\t" << item.second ;
    }
    cout << endl;

}

int CountFrequency::getSampleSum() const {
    return m_sum;
}
