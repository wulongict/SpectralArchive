//
// Created by wulong on 6/25/18.
//

#ifndef MYTOOL_SOCKETSERVER_H
#define MYTOOL_SOCKETSERVER_H



// to solve : error: 'inet_addr' 'sockaddr_in' was not declared in this scope
#include <netinet/in.h>

#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;
class CSpectralArchive;
class CSocketServerSummary;
// struct sockaddr_in;

bool update_queryindex_from_post_message(string x, int &topN, int &queryindex, int &calcEdge, int &nprobe,
                                         vector<uint16_t> &mzspec, int &visualization);

class CSocketServer {
    int forClientSockfd;
    int sockfd;
    struct sockaddr_in serverInfo;
    struct sockaddr_in clientInfo;
    int addrlen;
    int m_topN;
    const int PeakNumPerSpec;
    CSpectralArchive &m_archive;
    shared_ptr<CSocketServerSummary> m_summary;
public:
    CSocketServer(int maxconnections, int port, CSpectralArchive &archive, int topn);
    void startService();
    ~CSocketServer(){  }
};

#endif //MYTOOL_SOCKETSERVER_H
