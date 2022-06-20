//
// Created by wulong on 11/14/18.
//

#ifndef MYTOOL_FASTCGIINTERFACE_H
#define MYTOOL_FASTCGIINTERFACE_H

#include <fcgio.h>
#include <string>
#include <fstream>
using namespace std;
class CSocketServerSummary;
class CSpectralArchive;



namespace HTTP
{
    const char GATEWAY_INTERFACE[] = "GATEWAY_INTERFACE";
    const char SERVER_SOFTWARE[] = "SERVER_SOFTWARE";
    const char QUERY_STRING[] = "QUERY_STRING";
    const char REQUEST_METHOD[] = "REQUEST_METHOD";
    const char CONTENT_TYPE[] = "CONTENT_TYPE";
    const char CONTENT_LENGTH[] = "CONTENT_LENGTH";
    const char SCRIPT_FILENAME[] = "SCRIPT_FILENAME";
    const char SCRIPT_NAME[] = "SCRIPT_NAME";
    const char REQUEST_URI[] = "REQUEST_URI";
    const char DOCUMENT_URI[] = "DOCUMENT_URI";
    const char DOCUMENT_ROOT[] = "DOCUMENT_ROOT";
    const char SERVER_PROTOCOL[] = "SERVER_PROTOCOL";
    const char REMOTE_ADDR[] = "REMOTE_ADDR";
    const char REMOTE_PORT[] = "REMOTE_PORT";
    const char SERVER_ADDR[] = "SERVER_ADDR";
    const char SERVER_PORT[] = "SERVER_PORT";
    const char SERVER_NAME[] = "SERVER_NAME";
    const char HTTP_HOST[] = "HTTP_HOST";
    const char HTTP_CONNECTION[] = "HTTP_CONNECTION";
    const char HTTP_CONTENT_LENGTH[] = "HTTP_CONTENT_LENGTH";
    const char HTTP_PRAGMA[] = "HTTP_PRAGMA";
    const char HTTP_CACHE_CONTROL[] = "HTTP_CACHE_CONTROL";
    const char HTTP_ORIGIN[] = "HTTP_ORIGIN";
    const char HTTP_USER_AGENT[] = "HTTP_USER_AGENT";
    const char HTTP_CONTENT_TYPE[] = "HTTP_CONTENT_TYPE";
    const char HTTP_ACCEPT[] = "HTTP_ACCEPT";
    const char HTTP_REFERER[] = "HTTP_REFERER";
    const char HTTP_ACCEPT_ENCODING[] = "HTTP_ACCEPT_ENCODING";
    const char HTTP_ACCEPT_LANGUAGE[] = "HTTP_ACCEPT_LANGUAGE";
    const char HTTP_COOKIE[] = "HTTP_COOKIE";

}


class CFastCGIServer {
    const int m_PeakNumPerSpec = 50; // todo: wrong
    streambuf *m_cin_streambuf;
    streambuf *m_cout_streambuf;
    streambuf *m_cerr_streambuf;

    FCGX_Request m_request;

    CSpectralArchive & m_archive;
    shared_ptr<CSocketServerSummary> m_summary;
    bool verbose;
    int m_id;
public:
    // 1. return -Werror=return-type
    // 2. reference initialization as member list.
    CFastCGIServer(int maxconnections, int port, CSpectralArchive &archive, int topn, int id, shared_ptr<CSocketServerSummary> socketSummary);
    void response(string &message, string contentType);
    void startFastCGIServer();
    string getContent();
    void searchPeptideSequence(string &content);
    void searchQueryId(string &content);
    void addToRemark(string &content);
    void getRemark(string &uristr);
    void getPageWithId(string &uristr, string& html_template_str);
    void getPeakListWithId(string &uristr);
    void identification(string &content);
    string getEnv(string key, FCGX_Request &request);
};

#endif //MYTOOL_FASTCGIINTERFACE_H
