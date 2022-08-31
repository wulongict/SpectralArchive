//
// Created by wulong on 10/29/15.
//

#ifndef PROJECT_UTIL_H
#define PROJECT_UTIL_H


#include <iostream>
#include <thread>
#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>
#include <iterator>

using namespace std;
const double proton_mass = 1.007276;
const double EPSILON = 1e-8;


/**
   *  @brief Shuffle the elements of a sequence using a uniform random
   *         number generator.
   *  @ingroup mutating_algorithms
   *  @param  __first   A forward iterator.
   *  @param  __last    A forward iterator.
   *  @param  __g       A UniformRandomNumberGenerator (26.5.1.3).
   *  @return  Nothing.
   *
   *  Reorders the elements in the range @p [__first,__last) using @p __g to
   *  provide random numbers.
  */
  template<typename _RandomAccessIterator,
	   typename _UniformRandomNumberGenerator>
    void
    gcc5shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last,
	    _UniformRandomNumberGenerator&& __g)
    {
      // concept requirements
      __glibcxx_function_requires(_Mutable_RandomAccessIteratorConcept<
	    _RandomAccessIterator>)
      __glibcxx_requires_valid_range(__first, __last);

      if (__first == __last)
	return;

      typedef typename iterator_traits<_RandomAccessIterator>::difference_type
	_DistanceType;

      typedef typename std::make_unsigned<_DistanceType>::type __ud_type;
      typedef typename std::uniform_int_distribution<__ud_type> __distr_type;
      typedef typename __distr_type::param_type __p_type;
      __distr_type __d;

      for (_RandomAccessIterator __i = __first + 1; __i != __last; ++__i)
	std::iter_swap(__i, __first + __d(__g, __p_type(0, __i - __first)));
    }

class ICVerboseMsg
{
public:
    virtual ~ICVerboseMsg()= default;
    virtual void add(string msg)=0;
    virtual string getMessage()=0;
};

class CDummyMessageCollector: public ICVerboseMsg
{
public:
    CDummyMessageCollector()= default;
    void add(string msg) override{}
    string getMessage() override{return "";}
};

class CVerboseMessage: public ICVerboseMsg
{
    string m_verbose_message;
public:
    CVerboseMessage()= default;
    void add(string msg) override
    {
        m_verbose_message += msg;
    }
    string getMessage() override
    {
        return m_verbose_message;
    }
    ~CVerboseMessage() override{
        cout << m_verbose_message << flush;
    }
};



class CKeyValuesParser {
    map<string, string> m_data;
    shared_ptr<ICVerboseMsg> m_vmsg;
public:
    void findValue(string &nextstr, string &value, const string& pairseparator);
    explicit CKeyValuesParser(string &keyvaluestr, const string& keyvalueseparator="=", const string& pairseparator=" ",bool verbosity=false);
    void setKeyValue(string &key, string &value);
    void display();
    string getvalue(const string& key) {
        if(m_data.find(key)!=m_data.end())  return m_data.at(key);
        else{
            m_vmsg->add("Error: key "+ key + " does not exist!\n");
            return "";
        }
    }

};

bool find_value(vector<long> &a, long val, int &idx);

class CHistogram {
    vector<int> m_hist;
public:
    explicit CHistogram(int N);
    void add_data(int d);
    void display();
};

class CountFrequency{
    // tested!
    map<int, int> m_counts;
    int m_sum;
public:
    CountFrequency();
    void add_data(int d);
    int getFreq(int key);
    int getSampleSum() const;
    double getEntropy(bool normalized=false,bool verbosity=false);
    void print();
    int totalInsertion(){
        int sum = 0;
        for(auto &x:m_counts) sum += x.second;
        return sum;
    }
};


class CANSIConsole {
    struct ansi_colors {
        int front;
        int back;
        ansi_colors(int f, int b)
        {
            front = f;
            back = b;
        }
        string toStr() const
        {
            return string("\033[38;5;") + to_string(front) +";48;5;" + to_string(back) + "m";
        }
    };
    vector<ansi_colors> m_colors_inside;
    bool active;
public:
    enum color {DEFAULT, RED, GREEN, BLUE, YELLOW};
    CANSIConsole();
    string getColorStr(const string& str, color x);
    void reset();

    void set(color x=DEFAULT);
    ~CANSIConsole();

};

class CInstantColorConsole
{
    CANSIConsole a;
public:
    CInstantColorConsole()=default;

    CInstantColorConsole(string &x, bool newline,CANSIConsole::color c=CANSIConsole::DEFAULT)    {
        print(x,newline,c);

    }
    void print(const string &x,bool newline=true, CANSIConsole::color c=CANSIConsole::DEFAULT)    {
        a.set(c);
        cout << x ;
        a.reset();
        if(newline) cout << endl;
    }

};

void split_string(string &line, vector<string> &tokens);

void split_string(const string &line, vector<string> &tokens, char delim);

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6){

    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


template<class T>
std::string FormatWithCommas(T value){
    std::stringstream ss;
    ss.imbue(std::locale("en_US.UTF-8"));
    ss << std::fixed << value;
    return ss.str();
}

string getTimeStemp();

void get_FDR_CorrectNum(vector<double> & tProbs, vector<double> & dProbs, vector<tuple<double, double>> & FDR_CorrectNum);
void trim_space_only(std::string &str);

string make_timestring();

//
//template<typename T>
//void exportTable(vector<T> &data, const string& filename, char delimitor);

template<typename T>
void exportTable(vector<T> &data, const string &filename, char delimitor) {
    ofstream fout(filename.c_str(), ios::out);
    for (int i = 0; i < data.size(); i++) {
        fout << data[i];
        if (i != data.size() - 1) fout << delimitor;

        fout << endl;
    }
    cout << "saved to file: " << filename << endl;

}
template<typename T>
void exportTable(vector<vector<T>> &data, const string& filename, char delimitor) {
    ofstream fout(filename.c_str(), ios::out);
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[j].size(); j++) {
            fout << data[i][j];
            if (j != data[j].size() - 1) fout << delimitor;
        }
        fout << endl;
    }
    fout.close();
}



template<typename T>
string to_string(const string& start, string delimiter, const T& t){
    ostringstream oss;
    oss << t << delimiter ;
    return start + oss.str();
}


template<typename T>
string to_string_with_oss(ostringstream &oss, const string& delimiter, const T& t)
{
    oss << t ;
    return oss.str();

}

template<typename T, typename ... Ts>
string to_string_with_oss(ostringstream &oss, string delimiter, const T& t, const Ts & ... others)
{
//    cout << std::is_function<decltype(t)>::value << endl;

    oss << t ;
    if(strcmp(typeid(T).name(),typeid(decltype(std::fixed)).name())==0
    or strcmp(typeid(T).name(),typeid(decltype(std::setprecision(9))).name())==0)
    {
        // find std::fixed or setprecision
    } else{
        oss << delimiter;
    }
    return to_string_with_oss(oss, delimiter,others...);
}



template<typename T, typename ... Ts>
string to_string(string start, string delimiter, const T &t, const Ts & ... others){
    ostringstream oss;
//    oss  << t << delimiter ;
    return  start+to_string_with_oss(  oss, delimiter, t,others...);
}





template<typename T>
void print(const string& delimiter, const T& t){
    cout << t;
}

template<typename T, typename ... Ts>
void print(string delimiter, const T &t, const Ts & ... others)
{
    cout << delimiter << t;
    print(delimiter,others...);
}


template< typename T>
void print_on_fixed_width(int width,const T& t){
    cout << setw(width) << setfill(' ') << t ;
}

template< typename T, typename... Ts>
void print_on_fixed_width(int width,const T &t, const Ts&... ts){
    cout << setw(width) << setfill(' ') << t ;
    print_on_fixed_width(width,  ts ...);
}

string readlinesfromfile(const string& filename);



class CPosMass
{
    map<int, double> &m_pos_mass;
public:
    explicit CPosMass(map<int, double> &pos_mass): m_pos_mass(pos_mass)
    {}
    string toString()    {
        string pos_mass_mapstr;
        for(auto & m_pos_mas : this->m_pos_mass)        {
            pos_mass_mapstr += to_string_with_precision(m_pos_mas.second,4) + "@"+to_string(m_pos_mas.first)+"|";
        }
        if(pos_mass_mapstr.empty())        {
            pos_mass_mapstr = "UNMODIFIED";
        }
        return pos_mass_mapstr;
    }
};


struct modification
{
    double cterm_mass;
    double nterm_mass;
    map<int, double> pos_mass;

    string getModStr()    {
        return CPosMass(pos_mass).toString();
    }
};

template<typename T>
T stringTo(string &s){
//    cout << s << endl;
    istringstream iss(s);
    T ret;
    iss >> ret;
//    cout << "type : " << typeid(T).name() << " value " << ret << endl;
    return ret;
}

class ICastStrings{
public:
    virtual ~ICastStrings()= default;
    virtual vector<string>& getResults() = 0;
    template<typename T>
    T to(std::size_t col)   {

        return stringTo<T>(getResults()[col]);
    }
    double toDouble(size_t col)    {
        return to<double>(col);
    }
    float toFloat(size_t col)    {
        return to<float>(col);
    }
    long toLong(size_t col)    {
        return to<long>(col);
    }
    int toInt(size_t col)    {
        return to<int>(col);
    }
    string toString(size_t col)    {
        return getResults()[col];
    }

    string operator[](size_t col)    {
        return getResults()[col];
    }
};

class CSqlSpecfileTableRow: public ICastStrings
{
    vector<string> &m_results;
public:
    enum HEADER {FILE_ID, FILENAME,START,END};
    CSqlSpecfileTableRow(vector<string> &results):m_results(results){}
    vector<string> & getResults() override {return m_results;}
};


// How this works!??
// All of the columns
// ID          FILEID      MS2COUNTS   PEPTIDE      SCORE       SCAN        CTERM       NTERM       MODIFICATION  PRECURSOR   CHARGE      RT
// PEPTIDEPROPHETPROB  IPROPHETPROB  ISDECOY     SIGNIFICANCE  PROTEIN               CE          ALTERPEPTIDE  RFSCORE

// ID  FILEID MS2COUNTS  PEPTIDE SCORE  SCAN CTERM  NTERM MODIFICATION  PRECURSOR CHARGE  RT
// PEPTIDEPROPHETPROB  IPROPHETPROB RFSCORE  ISDECOY SIGNIFICANCE  PROTEIN CE  ALTERPEPTIDE NEIGHBOR
class CSqlGtTableRow : public ICastStrings
{
public:
    vector<string> &m_results;
    enum HEADER { IDX,  FILEID, MS2COUNTS, PEPTIDE,
        SCORE, SCAN, CTERM, NTERM,MODIFICATION, PRECURSOR,
        CHARGE, RT, PEPTIDEPROPHETPROB, IPROPHETPROB,
        ISDECOY, SIGNIFICANCE, PROTEIN,CE,ALTERPEPTIDE,RFSCORE,FILENAME_APPENDED};
    explicit CSqlGtTableRow(vector<string> &results):m_results(results) {
//        cout << "This is it : row num: " << results.size() << " -- and our last items: " << PROTEIN << endl;
    }
    vector<string> & getResults() override {return m_results;}

    // todo: to be used later
    void toOstringStream(ostringstream &oss);

    void toOstringStreamNoId(ostringstream &oss);

    string getJsonNode(bool getfilename=false);
};



string getmodificationfrompeptidestring(string peptidestr, modification &mod);


const string CTable_unavailableEntry = "N/A";

class CTable {
    // Every table got its column; if not, put C0, C1, C2,.... as header
    //    const static string unavailableEntry;
    string m_filename;
    char m_delim;
    bool m_has_header;
    vector<vector<string>> m_table; // content
    vector<string> m_column_header; // this is special row.
    std::map<string,vector<int>> m_tableindex;
    map<string, int> m_header2col;
public:
    size_t m_row;
    size_t m_col;
    CTable();
    CTable(const string& filename, char delimitor, bool has_header, int skipNum = 0);
    void buildheader2column(bool rebuild=false);
    void build_table_index(int col);

    bool hasHeader() const;
    void appendHeader(const vector<string>& Header);
    string getHeaderByColumn(size_t column);
    int getColByHeader(const string& header);
    void setHeader(const vector<string> &header);

    unsigned long LengthofRow(size_t k);
    void resizeRow(int k, int newsize);

    int getRowByKey(const string& key, int col);    // key and key column
    string getEntry(int row, int col) const ;

    void addRow(const vector<string>& row);
    void setEntry(int row, int col, string value);
    void appendEntry(int row, const string& value);

    void Join(CTable &other);

    void saveAs(const string& filename, bool with_header, char delimtor);

    void print();
    void printRow(int i);
};

//template <typename T>
//void CTable::exportTable(vector<vector<T>> &data, string filename, char delimitor) ;

class Progress{
    long m_task_num;
    long m_task_num_finished;
    int m_percentage;
    const int m_print_percentage_gap;
public:
    explicit Progress(long task_num);
    Progress(long task_num, const string& task_name);
    void increase(int n=1);// todo: thread safe
    ~Progress()= default;
};

class CountProgress{
    int m_counts;
    int m_task_notify;
    int m_taskbatch;
    bool m_running;
    void init(int batchsize)    {
        m_taskbatch = batchsize;
        m_counts = 0;
        m_task_notify = m_taskbatch;
        m_running = true;
    }
public:
    CountProgress(int batchsize){
        init(batchsize);
        cout << "Progress: " << flush;
    }

    CountProgress(int batchsize, const string& task_name){
        init(batchsize);
        cout << "[" << task_name << "] Progress: " << flush;

    }
    void increase(int n=1)
    {
        m_counts += n;
        if(m_counts>=m_task_notify) {
            cout << "..." << FormatWithCommas(m_counts) << flush;
            m_task_notify += m_taskbatch;
        }
    }
    void stop(){
        if(m_running)
        cout << "..." << FormatWithCommas(m_counts) << "-DONE" << endl;
        m_running = false;
    }
    ~CountProgress(){stop();}

};


namespace File{

    void splitpath(const string& inputpath, string &path, string &file);
    void parent(string inputpath, string &path, string &folder);
    void splitname(const string& filename, string &filenameprefix, string &ext);
    template<typename T>
    void saveas(vector<T> data, const string& filename, bool verbose)
    {
        ofstream fout(filename, ios::out);
        copy(data.begin(),data.end(), ostream_iterator<T>(fout,"\n"));
        if(verbose)cout << "Data saved as : " << filename << endl;
    }
    bool isExist(const string &curOutputfile, bool beQuiet=false);
    struct CFile{
            string m_fullpathname;
        string path; // path do not contain "/".
        string filename;
        string ext;
        string basename;
        explicit CFile(const string& inputPath);
        bool isFileExist(bool beQuiet=false) const;
    };

    void getfilesize(const string &filename, long &filesize);
    char * loadFileToBuffer(const string & filename);


}


namespace statistic {
    template<typename T>
    T calcmean(vector<T> v) {
        T s = 0;
        for (int i = 0; i < v.size(); ++i) {
            s += v[i];
        }
        return s / v.size();
    }

    template <class T> T calcstd(vector<T> v) {
        T avg = calcmean(v);
        T sum_of_square = 0;
        for (int i = 0; i < v.size(); ++i) {
            sum_of_square += (v[i] - avg) * (v[i] - avg);
        }
        T std = 0;
        if(not v.empty())
        {
            std = sqrt(sum_of_square / v.size());
        }
        return std;
    }

    template <class T> T calcpow4Std(vector<T> v) {
        T avg = calcmean(v);
        T sum_of_pow4 = 0;
        for (int i = 0; i < v.size(); ++i) {
            sum_of_pow4 += (v[i] - avg) * (v[i] - avg)*(v[i] - avg) * (v[i] - avg);
        }
        T std = 0;
        if(not v.empty())
        {
            std = sqrt(sqrt(sum_of_pow4 / v.size()));
        }
        return std;
    }

    template <class T> T calcabsavg(vector<T> v) {
        T avg = calcmean(v);
        T sum_of_absavg = 0;
        for (int i = 0; i < v.size(); ++i) {
            sum_of_absavg += fabs(v[i] - avg);//
            // * (v[i] - avg)*(v[i] - avg) * (v[i] - avg);
        }
        T std = 0;
        if(not v.empty())
        {
            std = (sum_of_absavg / v.size());
        }
        return std;
    }

    template<typename T>
    T BHfdrThreshold(vector<T> &p_values, T fdrThreshold, bool useQvalue, bool verbose) {
        if(p_values.empty()) return 0;
        T threshold = 0, pvalueThresholdQvalue=0;
        sort(p_values.begin(), p_values.end());
        vector<T> FDR(p_values.size(), 0);
        for (int i = 0; i < p_values.size(); i++) {
            FDR[i] = p_values[i] * p_values.size() / (i + 1);
        }
        // get threshold
        int num_ssm_identified_with_fdr = 0;
        for (int j = 0; j < FDR.size(); j++) {
            if (FDR[j] < fdrThreshold) {
                threshold = p_values[j];
            } else {
                num_ssm_identified_with_fdr = j;
                break;
            }
        }
        // FDR will not be empty!! use size()-1
        int num_ssm_identified_with_qvalue = 0;
        for (int j = FDR.size()-1; j >=0 ; j--) {
            if (FDR[j] <= fdrThreshold) {
                pvalueThresholdQvalue = p_values[j];
                num_ssm_identified_with_qvalue = j ;

                break;
            } else {
                continue;
            }
        }
        if(verbose){
            cout << "[FDR] threshold of FDR  < 1% is : " << std::scientific << threshold
                 << " identified ssm " << num_ssm_identified_with_fdr << " total " << p_values.size()<< endl;
            cout << "[FDR] threshold of qvalue < 1% is : "<< std::scientific << pvalueThresholdQvalue
                 << "identified ssm " << num_ssm_identified_with_qvalue << " total " << p_values.size()<< endl;
        }

        return useQvalue? pvalueThresholdQvalue: threshold;
    }



    double calcGeneralizedMean(const vector<double>& v, double p);

    double getNonzeroMin(const vector<double> &v);

    void checkAndFixVector(vector<double> &v, double &min, double &max);

    double calcHarmonicMean(vector<double> v);

    double calcGeometricMean(vector<double> v);

    double calcArithmeticMean(const vector<double>& v);

    double calcQuadraticMean(const vector<double>& v);

    double calcCubicMean(const vector<double>& v);

}


class SimpleTimer {
private:
    std::chrono::time_point<std::chrono::steady_clock> m_start;
//    int m_start;
    std::chrono::time_point<std::chrono::steady_clock> m_end;
//    int m_end;
//    double m_used;
    std::chrono::duration<double> m_used;
    string m_taskname;
    bool m_verbose;
public:


    // avoid implicit conversion from const char * to bool.
    explicit SimpleTimer(char * x, bool verbose=true);
    explicit  SimpleTimer(bool verbose=true);
    explicit SimpleTimer(string taskname, bool verbose=true);

    double stop();
    double restart(const string& taskname="");
    double secondsElapsed();

    ~SimpleTimer();

};

int getProperThreads(int threadnum=-1);


void TestMatrix();




class CMyMatrix {
private:
    long m_row;
    long m_col;
public:
    long getM_row() const;

    long getM_col() const;

private:
    long m_Size;
    double *m_entries;
public:
    CMyMatrix();
    CMyMatrix(const CMyMatrix &other);
    CMyMatrix(long r, long c);

    CMyMatrix(long r, long c, const double *matrix);
    ~CMyMatrix();
    void Print() const;

    void initialize(double val = 0, long r = -1, long c = -1);
    void set(long i, long j, double value);
    double get(long i, long j) const;

    explicit CMyMatrix(string filename);

    void outputBinary(const string& outputfilename);
    void outputAsText(const string& outputfilename);
    void mergeMatrix(vector<CMyMatrix> &vMatrix, const string& mergeMethod);
    void HomonicMatrix(vector<CMyMatrix> &vMatrix);
    void gmOfMatrix(vector<CMyMatrix> &vMatrix);
};


// do we still need this?
class CFolder{
    string m_folder;
public:
    CFolder();
    explicit CFolder(string folder);
    CFolder(const CFolder &cf);
    CFolder& operator=(const CFolder &cf);
    ~CFolder()= default;

    CFolder getParent();
    string getFolderName();
    string toStr();
    void print();
};

class CPath{
    string m_fullpath;

public:
    string m_folder;
    string m_filename;
    string m_filename_no_ext;
    string m_ext;
public:
    explicit CPath(string fullpath);
    void print();

    CFolder getParentFolder() const;

private:
    void parse();
};


#endif //PROJECT_UTIL_H
