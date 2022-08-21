//
// Created by wulong on 8/18/22.
//

#include "CTimerSummary.h"
#include <iostream>

CTimeSummary * CTimeSummary::g_time_summary = nullptr;

CTimeSummary::CTimeSummary() {}

CTimeSummary *CTimeSummary::getInstance() {
    if(g_time_summary == nullptr){
        g_time_summary = new CTimeSummary;
    }
    return g_time_summary;
}

void CTimeSummary::startTimer(std::string name) {
    if(m_record.count(name)==0){
        m_record[name] = timerTuple();
        m_record[name].m_used = std::chrono::duration<double>(0);
    }
    m_record[name].m_start = std::chrono::steady_clock::now();
}

void CTimeSummary::pauseTimer(std::string name) {
    if(m_record.count(name)==0){
        return;
    }
    m_record[name].m_used += (std::chrono::steady_clock::now()-m_record[name].m_start);

}

double CTimeSummary::timeUsed(std::string name) {
    if(m_record.count(name)==0) {
        return 0.0;
    }
    return m_record[name].m_used.count();
}

void CTimeSummary::print(std::ostream &os, long num_tasks) {
    os << "name_of_timer\ttime_used(s)\tavgtime_used(ms)" << std::endl;
    for(auto it = m_record.begin(); it != m_record.end(); it ++){
        double elapsed = timeUsed(it->first);
        os << it->first << "\t" << elapsed << "\t" << elapsed/num_tasks*1000 << std::endl;

    }
}
