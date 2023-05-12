//
// Created by wulong on 2/1/18.
//

#ifndef MYTOOL_DEBUGMODE_H
#define MYTOOL_DEBUGMODE_H


#include <string>
// define some global vaiables
#include <atomic>

extern std::atomic<bool> g_quit_flag;


class CDebugMode {

private:
    bool m_debug;

    static CDebugMode * pDebug;
    CDebugMode(){
        m_debug = false;
        rank=0;
        outputNullPep=false;
        useScanCharge=false;
    }
public:
    int rank;
    bool outputNullPep;
    bool useScanCharge;
    static CDebugMode * callDebug();
    static void releasePtr();
    void setDebug(bool debug)
    {
        m_debug = debug;
    }
    bool getMdebug()
    {
        return m_debug;
    }
};

void initlog(std::string logfile, std::string logname);



#endif //MYTOOL_DEBUGMODE_H
