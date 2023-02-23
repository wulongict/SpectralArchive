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

// two more loggers added, namely file and stdout.
void initlog(std::string logfile, std::string logname)
{
    
    // write log to file <logfile> and registered as  <logname>
    auto filesink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile.c_str(),
                                                                           1024 * 1024 * 1024, 1);
    
    auto stdsink = std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(filesink);
    sinks.push_back(stdsink);
    if (spdlog::get(logname) == nullptr)
    {
        auto combined_logger = std::make_shared<spdlog::logger>(logname, std::begin(sinks), std::end(sinks));
        auto filelogger = std::make_shared<spdlog::logger>("file", filesink);
        auto stdlogger = std::make_shared<spdlog::logger>("stdout", stdsink);
        if (combined_logger == nullptr or filelogger==nullptr or stdlogger == nullptr)
        {
            std::cerr << "Fail to create loggers. " << std::endl;
        }
        spdlog::register_logger(combined_logger);
        spdlog::register_logger(filelogger);
        spdlog::register_logger(stdlogger);
    }
    else
    {
        std::cout << "logger already exist" << std::endl;
    }
}
