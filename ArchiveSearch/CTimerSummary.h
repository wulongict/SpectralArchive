//
// Created by wulong on 8/18/22.
//

#ifndef MYTOOL_CTIMERSUMMARY_H
#define MYTOOL_CTIMERSUMMARY_H

#include <chrono>
#include <string>
#include <map>
#include <ostream>

class CTimeSummary{
    struct timerTuple{
        std::chrono::time_point<std::chrono::steady_clock> m_start;
        std::chrono::duration<double> m_used;
    };
    static CTimeSummary * g_time_summary;
    std::map<std::string, timerTuple> m_record;
    CTimeSummary();
public:
    static CTimeSummary * getInstance();
    void startTimer(std::string name);
    void pauseTimer(std::string name);
    double timeUsed(std::string name);
    void print(std::ostream &os, long num_tasks);


};

#endif //MYTOOL_CTIMERSUMMARY_H
