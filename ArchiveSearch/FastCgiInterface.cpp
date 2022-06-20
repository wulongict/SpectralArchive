//
// Created by wulong on 11/14/18.
//

#include <iostream>
#include <spdlog/spdlog.h>
#include <fcgio.h>
#include "const_params.h"
#include "../librarymsms/ICMzFile.h"
#include "../librarymsms/Util.h"
#include "FastCgiInterface.h"
#include "CSpectralArchive.h"
#include "../librarymsms/ProteomicsDataTypes.h"

using namespace std;



// Maximum bytes
const unsigned long STDIN_MAX = 1000000;

/**
 * Note this is not thread safe due to the static allocation of the
 * content_buffer.
 */
string get_request_content(const FCGX_Request &request) {
    char *content_length_str = FCGX_GetParam(HTTP::CONTENT_LENGTH, request.envp);
    unsigned long content_length = STDIN_MAX;

    if (content_length_str) {
        content_length = strtol(content_length_str, &content_length_str, 10);
        if (*content_length_str) {
            cerr << "Can't Parse 'CONTENT_LENGTH='"
                 << FCGX_GetParam(HTTP::CONTENT_LENGTH, request.envp)
                 << "'. Consuming stdin up to " << STDIN_MAX << endl;
        }

        if (content_length > STDIN_MAX) {
            content_length = STDIN_MAX;
        }
    } else {
        // Do not read from stdin if CONTENT_LENGTH is missing
        content_length = 0;
    }

    char *content_buffer = new char[content_length];
    cin.read(content_buffer, content_length);
    content_length = cin.gcount();

    // Chew up any remaining stdin - this shouldn't be necessary
    // but is because mod_fastcgi doesn't handle it correctly.

    // ignore() doesn't set the eof bit in some versions of glibc++
    // so use gcount() instead of eof()...
    do cin.ignore(1024); while (cin.gcount() == 1024);

    string content(content_buffer, content_length);
    delete[] content_buffer;
    return content;
}


void CFastCGIServer::startFastCGIServer() {
    cout << "FCGX_IsCGI" << FCGX_IsCGI() << endl;
    int isOK = FCGX_Init();
    if (isOK != 0){
        cout << "FCGX_Init() fail " << isOK << endl;
    } else{
        cout << "FCGX_Init() ready " << endl;
    }
    FCGX_InitRequest(&m_request, 0, 0);


    // The loop of the server starts from here.
    // How can we make multiple thread call ?
    while (FCGX_Accept_r(&m_request) == 0) {
        char **env = m_request.envp;
        if (verbose) {
            while (*(++env)) {
                puts(*env);
                // cout << *env << endl;
            }
        }

        string uristr = getEnv(HTTP::REQUEST_URI, m_request);
        string methodstr = getEnv(HTTP::REQUEST_METHOD, m_request);

//        cout << "Stop here ??? URI=" << uristr << " Method=" << methodstr << endl;
        spdlog::get("A")->info("{} {} {}:{}", getEnv(HTTP::REQUEST_METHOD, m_request),
                               getEnv(HTTP::HTTP_REFERER, m_request),
                // Sometimes,Google Chrome does not send HTTP_REFERER ... so we may have empty string
                               getEnv(HTTP::REMOTE_ADDR, m_request),
                               getEnv(HTTP::REMOTE_PORT, m_request));

        string content = getContent();

        if (methodstr == "POST") {
            spdlog::get("A")->info("POST content: {}", content);
//            cout << "POST Content: " << content << endl;
            if (uristr.find("/id/") != string::npos) {
                this->searchQueryId(content);
            } else if (uristr.find("/peptideseq") != string::npos) {
                this->searchPeptideSequence(content);
            } else if (uristr.find("/identification") != string::npos) {
                // it is time to use this
                this->identification(content);
            } else if (uristr.find("/remark") != string::npos) {
                this->addToRemark(content);
//                this->identification(content);
            } else {
                cout << "URI not found: " << uristr << endl;
            }
            if (content.find("STOPIT") != string::npos) {
                cout << "Find stop in the content" << endl;
                break;
            }
        } else if (methodstr == "GET") {
            spdlog::get("A")->info("Get URI: '{}'", uristr);

            if (uristr.find("/id/") != string::npos) {
                // Long: update html_template_str every time... [to be optimized in the future]
                string htmlfile = "/data/wulong/bitbucket/lorikeet/faiss_ivf256_pq16_gpu.index.json.html";
                string html_template_str = readlinesfromfile(htmlfile);
                this->getPageWithId(uristr, html_template_str);
            } else if (uristr.find("/identification") != string::npos) {
                // it is time to use this
                this->identification(content);
            } else if (uristr.find("/remark/") != string::npos) {
                this->getRemark(uristr);
                // add get remark code
            } else if (uristr.find("/spectrum?id") != string::npos) {
                this->getPeakListWithId(uristr);
            }
        } else {
            cout << "[Error] Get Request other than: GET/POST! ==> " << methodstr << endl;
        }
    }
}


CFastCGIServer::CFastCGIServer(int maxconnections, int port, CSpectralArchive &archive, int topn, int id, shared_ptr<CSocketServerSummary> socketSummary)
        : m_archive(archive), m_request() {
    if(socketSummary == nullptr){
        m_summary = make_shared<CSocketServerSummary>();
    }else{
        m_summary = socketSummary;
    }


    verbose = false;
    m_cin_streambuf = cin.rdbuf();
    m_cout_streambuf = cout.rdbuf();
    m_cerr_streambuf = cerr.rdbuf();
    m_id = id;

//    startFastCGIServer();
}
// More options can be added,
// such as Cache-Control: max-age=3600
//
void CFastCGIServer::response(string &message, string contentType) {
//    cout << "MESSAGE is : " << message << endl;
    fcgi_streambuf cout_fcgi_streambuf(m_request.out);
    fcgi_streambuf cerr_fcgi_streambuf(m_request.err);

    cout.rdbuf(&cout_fcgi_streambuf);
    cerr.rdbuf(&cerr_fcgi_streambuf);
    cout << "Content-type: " << contentType << "\r\n"
            << "Cache-Control: max-age=3600" << "\r\n"
         << "\r\n"
         << message;

    // Note: the fcgi_streambuf destructor will auto flush
    // restore stdio streambufs

    cout.rdbuf(m_cout_streambuf);
    cerr.rdbuf(m_cerr_streambuf);
    if (verbose) cout << "Json string send to browser: size = " << message.size() / 1024 << "KB" << endl;
}

string CFastCGIServer::getContent() {
    fcgi_streambuf cin_fcgi_streambuf(m_request.in);
    cin.rdbuf(&cin_fcgi_streambuf);

    string content = get_request_content(m_request);
    cin.rdbuf(m_cin_streambuf);
    return content;
}

struct contentParser {
    int calcEdge;
    int nprobe;
    long queryindex;
    int topN;
    vector<uint16_t> mzspec;
    int visualization;
    string remarks;

    contentParser() {
        visualization = 0; // NO
        calcEdge = 0; // NO

        nprobe = 256;
        topN = 10;
        queryindex = -1;  // NO
        remarks = "NONE";
    }


    bool update_int_through_string(string x, int &val, string key) {
        int pos = x.find(key + "=");
        if (pos == string::npos) {
            return false;
        } else {
            string info = x.substr(pos + key.length() + 1);
            int pos = info.find_first_not_of("0123456789");
            val = atoi(info.substr(0, pos).c_str());
            return true;
        }
    }

    bool update_remarks_through_string(string x, string &remarks, string endwith = "$") {
        string key = "SpecRemarks=";
        int pos = x.find(key);
        if (pos == string::npos) {
            return false;
        } else {
            string info = x.substr(pos + key.length());
            int pos = info.find_first_of(endwith);
            if (pos == string::npos) {
                return false;
            }
            remarks = info.substr(0, pos);
        }
        return true;
    }

    bool update_long_through_string(string x, long &val, string key) {
        int pos = x.find(key + "=");
        if (pos == string::npos) {
            return false;
        } else {
            string info = x.substr(pos + key.length() + 1);
            int pos = info.find_first_not_of("0123456789");
            val = atol(info.substr(0, pos).c_str());
            return true;
        }
    }

    bool update_spec_from_string(string &x, vector<uint16_t> &queryspec) {
        int pos = x.find("SPEC=");
        if (pos == string::npos) {
            return false;
        } else {
            string info = x.substr(pos + 5); // key.length() = 4 SPEC +1 = 5
            int pos = info.find_first_not_of("0123456789,");
            string values = info.substr(0, pos);
            vector<string> tokens;
            split_string(values, tokens, ',');

            const int MAX_PEAK_ALLOWED = 50;
            if (tokens.size() > MAX_PEAK_ALLOWED) {
                cout << "[Warning] More than 50 peaks in given spectrum!!" << endl;
            }
            queryspec.assign(MAX_PEAK_ALLOWED, 0);
            for (int i = 0; i < tokens.size() and i < MAX_PEAK_ALLOWED; i++) {
                int val = atoi(tokens[i].c_str());
                if (val >= 0 and val <= UINT16_MAX) {
                    queryspec[i] = val;
                } else {
                    cout << "Invalid value: " << val << ". Replaced with 0 and break now" << endl;
                    break;
                }
            }
            return true;
        }
    }

    void parse(string &content) {
        bool a = update_int_through_string(content, topN, "TOPN");
        bool b = update_long_through_string(content, queryindex, "QUERYID");
        bool c = update_int_through_string(content, calcEdge, "EDGE");
        bool d = update_int_through_string(content, nprobe, "NPROBE");
        bool f = update_spec_from_string(content, mzspec);
        bool g = update_int_through_string(content, visualization, "VISUALIZATION");
        bool h = update_remarks_through_string(content, remarks, "$");
//        return a or b or c or d or f or g;
    }
};

void CFastCGIServer::addToRemark(string &content) {
    if (content.size() > 0) {
        contentParser ps;
        ps.parse(content);
        string message = "DONE!";
        if (ps.queryindex >= -1 and ps.queryindex < m_archive.size()) {
            spdlog::get("A")->info("add remark: {}", ps.remarks);
            m_archive.addRemark(ps.queryindex, ps.remarks);
//            m_archive.searchQuery(ps.queryindex, message, ps.topN, ps.calcEdge, ps.nprobe, ps.mzspec, ps.visualization==1);
            response(message, "text/html");
        } else {
            cout << "Undefined request " << endl;
        }
    }
}

void CFastCGIServer::searchQueryId(string &content) {
    if (content.size() > 0) {
        SimpleTimer st;
        m_summary->num_of_search_done++;

        contentParser ps;
        ps.parse(content);

        string message;

        // Note: query id can be -1
        // When first query id is -1, we search with the spectrum in query pointer!
        if (ps.queryindex >= -1 and ps.queryindex < m_archive.size()
            and ps.topN >= 1 and ps.topN <= MAX_TOPN_ALLOWED) {  // now -1 is meaningful...
            spdlog::get("A")->info("start searching: queryindex={}, topn={} edge={}", ps.queryindex, ps.topN,
                                   ps.calcEdge);
            m_archive.searchQuery(ps.queryindex, message, ps.topN, ps.calcEdge, ps.nprobe, ps.mzspec,
                                  ps.visualization == 1);
            response(message, "text/html");
            m_summary->update(st.stop());
            cout << "\n    Server Id    " << m_id << endl;
            m_summary->print();
        } else {
            cout << "In construction --- to do ----" << endl;
        }
    }
}

// this is good implementation compared to another one. in socketserver
void CFastCGIServer::getPageWithId(string &uristr, string &html_template_str) {
    vector<string> param;
//    cout << __FUNCTION__ << " " << uristr << endl;

    split_string(uristr, param, '/');
//    cout << "-params " << endl;
//    for (auto each: param) cout << each << endl;
    if (param.size() > 2) {
        int queryindex = atoi(param[2].c_str());
        if(queryindex<0){
            // using magic number here
            queryindex = 16189;
        }

        string tmpstr = html_template_str;
        tmpstr.replace(tmpstr.find("16189"), string("16189").size(), to_string(queryindex));
        tmpstr.replace(tmpstr.find("7,849,350"), string("7,849,350").size(),
                       FormatWithCommas(m_archive.size()) + "(on " + m_archive.getPlatform() + ")");

        tmpstr.replace(tmpstr.find("August 14, 2018"), string("August 14, 2018").size(), __DATE__);

        response(tmpstr, "text/html");
    }
}

void CFastCGIServer::getRemark(string &uristr) {
    vector<string> param;
    //    cout << __FUNCTION__ << " " << uristr << endl;

    split_string(uristr, param, '/');
    //    cout << "-params " << endl;
    //    for (auto each: param) cout << each << endl;
    if (param.size() > 2) {
        int queryindex = atoi(param[2].c_str());
        if(queryindex<0){
            // using magic number here
            queryindex = 16189;
        }
        string remarks;
        m_archive.getRemark(queryindex, remarks);
        string tmpstr = "{'id': " + to_string(queryindex)+", 'remark': '" + remarks+  "'}";
        response(tmpstr, "text/json");
    }
}

//#include "randomSpecParser.h"
void CFastCGIServer::getPeakListWithId(string &uristr) {

    vector<string> param;
    split_string(uristr, param, '=');
    string newuristr=uristr.substr(uristr.find_first_of("?")+1);
//    int queryindex = atoi(param[1].c_str());
    CKeyValuesParser kvp(newuristr,"=","&");
    string chrome = kvp.getvalue("chrom");
    string updateSVG = kvp.getvalue("updateSVG");
    string rtTol = kvp.getvalue("rtTol");
    string mzTol = kvp.getvalue("mzTol");
    string mzOffset = kvp.getvalue("mzOffset");

    double rt_tol = 200, mz_tol = 0.03, mz_offset = 0.0;
    int charge = 1;
    if(rtTol != "")
    {
        rt_tol = stringTo<double>(rtTol);
    }
    if(mzTol != "")
    {
        mz_tol = stringTo<double>(mzTol);
    }
    if(mzOffset != "")
    {
        mz_offset = stringTo<double>(mzOffset);
    }
    string tmp = kvp.getvalue("charge");
    if(kvp.getvalue("charge") != ""){
        charge = stringTo<int>(tmp);
    }

    int queryindex = atoi(kvp.getvalue("id").c_str());

    cout <<"URI str: " << uristr << "\tqueryindex: " << queryindex << endl;
    if(queryindex < 0){
        queryindex = 0;
        cout << "-1 invalid query id changed to magic id: " << queryindex << endl;

    }

    // part one: mz string
    CMzSpec qspec = m_archive.getMzSpec(queryindex);
//    qspec.display();
    string messageMz = qspec.getpeakmsg();

    // part two: real spectra
    vector<double>   mz, intensity;
    bool getChromatogram = chrome == "1";
    shared_ptr<Chromatogram> chr = nullptr;
    map<int, shared_ptr<Chromatogram>> ChromatogramOfIsotopes = {
            {-1,nullptr},
            {0,nullptr},
            {1,nullptr},
            {2,nullptr},
            {3,nullptr}
    };
    bool updateChromSVG= updateSVG == "1";
    string mzrtmap="[]";
    m_archive.getRawSpec(queryindex, mz, intensity);
    m_archive.getChromatograms(queryindex, getChromatogram, updateChromSVG, mz_tol, charge, ChromatogramOfIsotopes,mzrtmap);

    cout << "get real peaks: size" << mz.size() << "\t" << intensity.size() << endl;
    if(mz.size()>0)
    {
    cout << "first peak: " << mz[0] << "\t" << intensity[0] << endl;

    }
    ostringstream realmz, realinten;
    for(int i = 0; i < mz.size(); i ++)    {
        if(i!=0){
            realmz<< " ";
            realinten << " ";
        }
        realmz << mz[i] ;
        realinten << intensity[i];
    }
    string realpeaks="{\"mz\": \""+ realmz.str()
            +"\", \n\"intensity\": \""
            + realinten.str() + "\"}";

    string message = "{\"mzs\": \""+messageMz+"\"," + "\n\"realpeaks\": " + realpeaks ;
    if(chr!= nullptr){
        message += ",\n\"chrom\": "+ chr->toJsonStr(rt_tol);
    }
    for(auto it=ChromatogramOfIsotopes.begin(); it != ChromatogramOfIsotopes.end(); it ++)
    {
        if(it->second != nullptr){
            string chrom_id = to_string(it->first);
            if(it->first == -1){
                chrom_id = "x";
            }
            message += ",\n\"chrom_" + chrom_id+"\": " + it->second->toJsonStr(rt_tol);
        }
    }
    if(not mzrtmap.empty()){
        message += ", \n\"mzrtmap\": " + mzrtmap;
    }
    message += "}";
    qspec.printMZStr();

    response(message, "text/json");
}


void CFastCGIServer::identification(string &content) {
    int found = content.find_last_not_of("\r\n\t ");
    if (found != string::npos) {
        cout << "content " << content << "\t" << found << endl;
        content = content.substr(0, found + 1);
    }
    // do something
    cout << "do identification " << content << endl;
    // identify the following spectra
    long queryidx = 8893380;
    int topN = 30;


    string jsonstring;
    vector<uint16_t> query;
    m_archive.searchQuery(queryidx, jsonstring, 30, 1, 20, query, false);
}

void CFastCGIServer::searchPeptideSequence(string &content) {
    int found = content.find_last_not_of("\r\n\t ");
    if (found != string::npos) {
        content = content.substr(0, found + 1);
    }
    // get peptide
    found = content.find("Peptide_Sequence=");
    string peptide;
    if (found != string::npos) {
        peptide = content.substr(found + string("Peptide_Sequence=").length());
    }

    string jsonmessage;

    if (peptide.empty()) {
        CKeyValuesParser kvp(content, "=", ";", true);
        jsonmessage = m_archive.searchFileName(kvp.getvalue("FILENAME"), kvp.getvalue("STARTSCAN"),
                                 kvp.getvalue("ENDSCAN"));

    } else {
        jsonmessage = m_archive.searchPeptide(peptide);
    }

    jsonmessage = "{\"nodes\": \n[" + jsonmessage + "]\n}";

    response(jsonmessage, "text/json");

}

string CFastCGIServer::getEnv(string key, FCGX_Request &request) {
    const char *envstr = FCGX_GetParam(key.c_str(), request.envp);
    if (envstr == nullptr) return key + ": N/A";
    return string(envstr);
}
