//
// Created by wulong on 6/25/18.
//

#include <sys/socket.h>
#include <arpa/inet.h>
#include "CSocketServer.h"
#include "../librarymsms/ProteomicsDataTypes.h"
#include "const_params.h"
#include "../External/spdlog-1.x/include/spdlog/spdlog.h"
#include "CSpectralArchive.h"   // archive
#include "../librarymsms/PeakList.h"
#include "../librarymsms/Util.h"

using namespace std;

int myrandom(int i) { return std::rand() % i; }

namespace PSMTable {
    enum resTableheader {
        index,
        filename,
        fileid,
        ms2_id,
        peptide,
        score,
        scan,
        cterm_mass,
        nterm_mass,
        modificationstr,
        precursor,
        charge,
        retentiontime,
        peptideprophetprob,
        iprophetprob
    };
}

namespace distance {
    double L2dist(float *x, float *y, int len) {
        double d = 0;
        for (int i = 0; i < len; i++) {
            d += (*(x + i) - *(y + i)) * (*(x + i) - *(y + i));
        }
        return sqrt(d);
    }

    double cosine(float *x, float *y, int len) {
        double cosine = 0;
        vector<float> zero(len, 0.0f);
        double xnorm = L2dist(x, zero.data(), len);
        double ynorm = L2dist(y, zero.data(), len);
        for (int i = 0; i < len; i++) {
            cosine += (*(x + i)) * (*(y + i));
        }
        if (fabs(xnorm) < EPSILON or fabs(ynorm) < EPSILON) {
            cosine = 0; // does this OK? maybe not!
        } else {
            cosine = cosine / xnorm / ynorm;
        }
        return cosine;

    }

    double test() {
        vector<float> x = {1, 2}, y = {3, 4};
        assert(fabs(cosine(x.data(), y.data(), 2)) < EPSILON);
        assert(fabs(L2dist(x.data(), y.data(), 2) - sqrt(8.0)) < EPSILON);
        return 0.0;

    }
}


bool update_int_through_string(string x, int &val, string key) {

    int pos = x.find(key + "=");
    if (pos == string::npos) {
        cout << "Key " << key << " NOT found in string: " << x << endl;
        return false;
    } else {
        string info = x.substr(pos + key.length() + 1);
        int pos = info.find_first_not_of("0123456789");
        val = atoi(info.substr(0, pos).c_str());

        cout << "Key: " << key << " is  udpated to " << val << " via string: " << x << endl;
        return true;
    }
}

bool update_spec_from_string(string &x, vector<uint16_t> &queryspec) {
    int pos = x.find("SPEC=");
    if (pos == string::npos) {
        cout << "key SPEC does not exist!" << endl;
        return false;
    } else {
        string info = x.substr(pos + 5); // key.length() = 4 SPEC +1 = 5
        int pos = info.find_first_not_of("0123456789,");
        string values = info.substr(0, pos);
        vector<string> tokens;
        split_string(values, tokens, ',');
        // must be 50 peaks.
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

bool update_queryindex_from_post_message(string x, int &topN, int &queryindex, int &calcEdge, int &nprobe,
                                         vector<uint16_t> &mzspec, int &visualization) {
    cout << "Looking for information: TOPN, QUERYID, and EDGE" << endl;
    bool a = update_int_through_string(x, topN, "TOPN");
    bool b = update_int_through_string(x, queryindex, "QUERYID");
    bool c = update_int_through_string(x, calcEdge, "EDGE");
    bool d = update_int_through_string(x, nprobe, "NPROBE");
    bool f = update_spec_from_string(x, mzspec);
    bool g = update_int_through_string(x, visualization, "VISUALIZATION");
    return a or b or c or d or f or g;

}


int safe_bind(int &sockfd, sockaddr_in &serverInfo) {
    int bind_status = bind(sockfd, (struct sockaddr *) &serverInfo, sizeof(serverInfo));
    if (bind_status == 0) {
        spdlog::get("A")->info("Successfully binding with socket!");
    } else if (bind_status == -1) {
        spdlog::get("A")->info("Fail to binding with socket!");
        // try for another ten times
        int n = 10;
        while (n--) {
            spdlog::get("A")->info("Trying to binding socket handle with socket address!");
            bind_status = bind(sockfd, (struct sockaddr *) &serverInfo, sizeof(serverInfo));
            if (bind_status == 0) {
                spdlog::get("A")->info("Successfully bining with socket!");
                break;
            } else if (bind_status == -1) {
                spdlog::get("A")->info("Fail to binding with socket! Try again after 10 seconds");

                sleep(10);
            }
        }
    } else {
        cout << "Unknow binding status: " << bind_status << endl;
    }
    return bind_status;
}

void set_timeout_for_send_recv(int &forClientSockfd, const int seconds) {
    struct timeval timeout = {seconds, 0};//3s
    int ret = setsockopt(forClientSockfd, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout));
    if (ret == 0) {
        spdlog::get("A")->info("Send time out: 3 seconds");
    } else {
        spdlog::get("A")->error("Fail to set timeout for Send: 3 seconds");
    }
    ret = setsockopt(forClientSockfd, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout));
    if (ret == 0) {
        spdlog::get("A")->info("Recieve time out: 3 seconds");
    } else {
        spdlog::get("A")->error("Fail to set timeout for Recieve: 3 seconds");
    }
}

class HTTPResponse {
    string messageheader;
    string message;
public:

    HTTPResponse(string messagebody, string responsetype) {
        messageheader = ""
                        "HTTP/1.1 200 OK\n"
                        //                  "Date: Sun, 18 Oct 2009 08:56:53 GMT\n"
                        //                  "Server: Apache/2.2.14 (Win32)\n"
                        //                  "Last-Modified: Sat, 20 Nov 2004 07:16:26 GMT\n"
                        //                  "ETag: \"10000000565a5-2c-3e94b66c2e680\"\n"
                        //                  "Accept-Ranges: bytes\n"
                        "Content-Length: %CONTENT_LEN%\n"
                        "Connection: close\n"
                        "Access-Control-Allow-Origin: *\n"
                        "Access-Control-Allow-Methods: POST, GET, OPTIONS\n"
                        "Access-Control-Allow-Headers: X-PINGOTHER, Content-Type\n"
                        "Access-Control-Max-Age: 15400\n"
                        "Content-Type: text/%RESPONSE_TYPE%\nCache-Control: max-age=36000\n"
                        //                  "X-Pad: avoid browser bug\n"
                        "\n";

        messageheader.replace(messageheader.find("%CONTENT_LEN%"), strlen("%CONTENT_LEN%"),
                              to_string(messagebody.size()));
        messageheader.replace(messageheader.find("%RESPONSE_TYPE%"), strlen("%RESPONSE_TYPE%"), responsetype);
        message += messageheader + messagebody;
    }

    string getMessage() { return message; }

};

class ReplyFileContent {
    enum file_response_table {
        QUERY_HEAD, RESPONSE_FILE
    };
    CTable *m_fileresponse;
public:
    ReplyFileContent(string configfile) {

        m_fileresponse = new CTable(configfile, ',', false, 0);
        m_fileresponse->build_table_index(ReplyFileContent::QUERY_HEAD);
        for (int i = 0; i < m_fileresponse->m_row; i++) {
            m_fileresponse->printRow(i);

        }
    }

    bool reply(string &message, string &responsetype, string queryhead) {
        spdlog::get("A")->info("Requesting \"{}\"", queryhead);
        // looking for root page
        int row = m_fileresponse->getRowByKey(queryhead, ReplyFileContent::QUERY_HEAD);
        if (row == -1) {
            message = "Not Found!";
            return false;
        }

        string replyfile = m_fileresponse->getEntry(row, ReplyFileContent::RESPONSE_FILE);
//        string htmlfile = "/data/wulong/data/honeybee/run.js";
        string filenameprefix, ext;
        File::splitname(replyfile, filenameprefix, ext);
        if (ext == "html") {
            responsetype = "html";
        } else if (ext == "css") {
            responsetype = "css";
        }
        //

        spdlog::get("A")->info("replywith  {} responsetype {}", replyfile, ext);


        string tmpstr = readlinesfromfile(replyfile);
        message = tmpstr;
        return true;
    }
};


CSocketServer::CSocketServer(int maxconnections, int port, CSpectralArchive &archive, int topn)
        : PeakNumPerSpec(50), m_archive(archive) {
    m_summary = make_shared<CSocketServerSummary>();
    spdlog::get("A")->debug("Getting to start the server!--1");

    m_topN = topn;
    sockfd = 0;
    forClientSockfd = 0;
    // Domain:
    // AF_INET: IPV4;
    // AF_INET6: IPV6;
    // Type:
    // SOCK_STREAM: TCP;
    // SOCK_DGRAM: UDP

    sockfd = socket(AF_INET, SOCK_STREAM, 0);

    if (sockfd == -1) {
        spdlog::get("A")->error("Fail to create a socket...");
        throw runtime_error("Fail to create the socket");
    }

    //socket的連線
//    struct sockaddr_in serverInfo,clientInfo;
    addrlen = sizeof(clientInfo);
    bzero(&serverInfo, sizeof(serverInfo));

    serverInfo.sin_family = PF_INET;
    serverInfo.sin_addr.s_addr = INADDR_ANY;

    serverInfo.sin_port = htons(port);

    // safe bind: socket and address
    if (safe_bind(sockfd, serverInfo) != 0) {
        throw logic_error("Fail to bind socket with internet address!");
    }
    if (maxconnections <= 0) {
        maxconnections = 10;
    }

    int listen_status = listen(sockfd, maxconnections);
    spdlog::get("A")->info("Listening on port {} with status {}; Max Connection allowed: {}",
                           port, listen_status, maxconnections);

    startService();


}

void CSocketServer::startService() {
    ReplyFileContent rfc("configfile");
    SimpleTimer st("server");
    while (true) {
        const int message_max_len = 25600;  // 25 KB information recieved!
        char inputBuffer[message_max_len] = {};

        m_summary->update(st.stop());
        m_summary->print();

        forClientSockfd = accept(sockfd, (struct sockaddr *) &clientInfo, &addrlen);
        m_summary->num_of_search_done++;
        st.restart("server");

        if (forClientSockfd == -1) {
            cout << "Fail to open socket for the connection!" << endl;
            spdlog::get("A")->error("Fail to open socket for connection!");
            continue;
        }
        // get timeout for send and recv
        set_timeout_for_send_recv(forClientSockfd, 9);

        int recv_bytes = recv(forClientSockfd, inputBuffer, sizeof(inputBuffer), 0);

        if (recv_bytes == 0 or recv_bytes == -1) {
            spdlog::get("A")->info("Client closed connection. Status={}", recv_bytes);
            close(forClientSockfd);
            continue;
        } else {
            spdlog::get("A")->info("bytes recieved: {}", recv_bytes);
        }
        inputBuffer[message_max_len - 1] = '\0';

        string x = inputBuffer;
        char *address_ip = inet_ntoa(clientInfo.sin_addr);
        string address_ipstr;
        spdlog::get("A")->info("IP={} Request=\"{}\"", address_ip, x);


        string message = "Not found";
        string responsetype = "html";
        spdlog::get("A")->info("Looking for POST in query information.  ");

        if (x.size() > 0 and x.find("POST") != string::npos) {
            int queryindex = -1;
            m_topN = 10;
            // update queryindex
            int calcEdge = 1;
            int nprobe = 256;
            vector<uint16_t> queryspec;
            int visualization = 0;
            // to be replaced!
            update_queryindex_from_post_message(x, m_topN, queryindex, calcEdge, nprobe,
                                                queryspec, visualization);

            if (queryindex >= 0 and queryindex < m_archive.size()
                and m_topN >= 1 and m_topN <= MAX_TOPN_ALLOWED) {
                spdlog::get("A")->info("start searching: queryindex={}, topn={} edge={}", queryindex, m_topN, calcEdge);
                // get json string response! --> message
                m_archive.searchQuery(queryindex, message, m_topN, calcEdge, nprobe, queryspec, visualization == 1);
//                searchQuery(queryindex, message, m_archive, m_topN);
            } else {

                cout << "Inconstruction!!!" << endl;

            }
        } else if (x.find("GET / HTTP") != string::npos) {
            string htmlfile = "/data/wulong/data/honeybee/faiss_ivf256_pq16_gpu.index.json.html";
            spdlog::get("A")->info("Get new query, reply with entry page: {}", htmlfile);
            cout << "Found GET / root " << endl;
            string tmpstr = readlinesfromfile(htmlfile);
            message = tmpstr;
        } else if (x.find("GET /id/") != string::npos) {

            int found = x.find_first_of("\r\n");
            if (found != string::npos) {
                string queryhead = x.substr(0, found);
                cout << "Looking for queryhead: " << queryhead << endl;
                vector<string> infos, param;
                split_string(queryhead, infos);
                for (auto eachinfo: infos) cout << eachinfo << endl;
                cout << "\n--" << flush << endl;
                string queryparam = infos[1];
                split_string(queryparam, param, '/');
                for (auto eachqueryparam: param) cout << eachqueryparam << endl;
                cout << "\n--" << flush << endl;
                if (param.size() > 2) {
                    cout << param[2] << "---> " << endl;
                    int queryindex = atoi(param[2].c_str());


                    string htmlfile = "/data/wulong/data/honeybee/faiss_ivf256_pq16_gpu.index.json.html";
                    string tmpstr = readlinesfromfile(htmlfile);
                    tmpstr.replace(tmpstr.find("1312081"), string("1312081").size(), to_string(queryindex));
                    tmpstr.replace(tmpstr.find("7,849,350"), string("7,849,350").size(),
                                   FormatWithCommas(m_archive.size()) + "(on " + m_archive.getPlatform() + ")"); // todo

                    tmpstr.replace(tmpstr.find("August 14, 2018"), string("August 14, 2018").size(), __DATE__);

                    message = tmpstr;
                } else {

                    cout << "Got empty " << endl;
                    message = "404 Not Found! id should not be empty";
                }


            }

        } else if (x.find("GET /spectrum?id") != string::npos) {

            int found = x.find_first_of("\r\n");
            if (found != string::npos) {
                string queryhead = x.substr(0, found);
                cout << "Looking for queryhead: " << queryhead << endl;
                vector<string> infos, param;
                split_string(queryhead, infos);
                string queryparam = infos[1];
                split_string(queryparam, param, '=');
                int queryindex = atoi(param[1].c_str());

                CMzSpec qspec = m_archive.getMzSpec(queryindex);
                qspec.display();

                message = qspec.getpeakmsg();
                qspec.printMZStr();
            }

        } else if (x.find("GET /run.js HTTP") != string::npos) {

            // looking for root page
            string htmlfile = "/data/wulong/data/honeybee/run.js";
            spdlog::get("A")->info("Requesting run.js reply with entry page: {}", htmlfile);


            string tmpstr = readlinesfromfile(htmlfile);
            message = tmpstr;


        } else if (x.size() > 0 and x.find("GET") != string::npos)// reply all the files
        {
            int found = x.find_first_of("\r\n");
            if (found != string::npos) {
                string queryhead = x.substr(0, found);
                cout << "Looking for queryhead" << queryhead << endl;
                rfc.reply(message, responsetype, queryhead);

            }
        }
        // todo: compress the message to save the internet bandwith

        message = HTTPResponse(message, responsetype).getMessage();

        spdlog::get("A")->info("IP={} Bytes to be sent:  {}", address_ip, message.size());
        int sent_bytes = send(forClientSockfd, message.c_str(), message.size(), 0);

        if (sent_bytes != message.size()) {
            spdlog::get("A")->error("Incomplete message sent: {} < {}", sent_bytes, message.size());
        } else {
            cout << "[Info] current query process is done. The server is waiting for new query...< < < = = =" << endl;
        }

        int status_exit = close(forClientSockfd);
        spdlog::get("A")->info("IP={} Socket closed with exitcode: {}", address_ip, status_exit);

        if (x.size() > 0 and x.find("CLOSE_SOCKET_NOW") != string::npos) {
            spdlog::get("A")->info("Socket will be closed now");
            message = "* [Notice] Server in maintainance\n";

            break;
        }


    }
}




