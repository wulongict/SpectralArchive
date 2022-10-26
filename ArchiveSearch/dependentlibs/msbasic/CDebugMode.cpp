//
// Created by wulong on 2/1/18.
//

#include "CDebugMode.h"
#include "spdlog/sinks/sink.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"
#include "spdlog/sinks/ansicolor_sink.h"

using namespace std;
using namespace spdlog;

CDebugMode * CDebugMode::pDebug = nullptr;

CDebugMode *CDebugMode::callDebug() {
    if(CDebugMode::pDebug == nullptr)
    {
        CDebugMode::pDebug = new CDebugMode();

    }
    return CDebugMode::pDebug;

}

void CDebugMode::releasePtr() {
    if(pDebug!=nullptr) {
        delete pDebug;
        pDebug = nullptr;
    }
}

#include <iostream>

void initlog(std::string logfile, std::string logname) {
    auto filelog = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile.c_str(),
                                                                          3*1024*1024*1024, 1);
    auto stdsink = std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(filelog);
    sinks.push_back(stdsink);
   if( spdlog::get(logname)==nullptr){
       auto combined_logger = std::make_shared<spdlog::logger>(logname, std::begin(sinks), std::end(sinks));
    if(combined_logger == nullptr){
        std::cout << "Logger not exist! " << std::endl;
    }
    spdlog::register_logger(combined_logger);
   }else{
       std::cout << "logger already exist" << std::endl;
   }
    

}
