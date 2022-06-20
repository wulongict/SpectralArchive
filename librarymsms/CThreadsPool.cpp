// Created by wulong on 4/24/19.
//

#include "CThreadsPool.h"
#include <vector>
#include <mutex>
#include <iostream>
#include "Util.h"
#include <thread>
#include <numeric>
#include <functional> // bind is here!
using namespace std;

CTaskPool::CTaskPool(vector<ICThreadTask *> &_vtask, bool verbose, bool progress ) :vTask(_vtask){
    nextjob = 0;
    m_verbose = verbose;
    m_progress = progress;
}

void CTaskPool::doTasks() {
    int taskID = getNextTaskID();
    while(taskID!=-1)  {
        run(taskID);
        taskID=getNextTaskID();
    }
}

// param: threadnum -1, will be reset to half of the cores.
void CTaskPool::start(int threadnum, string taskname) {
    if(threadnum == -1 or threadnum > getProperThreads())    {
        threadnum = getProperThreads()/2+1;
    }

    if(m_progress)cout << "["<< taskname<< " threads: " << threadnum<<  "] " << flush;
    vector<thread> vThreads;
    for(int j = 0; j < threadnum; j ++)    {
        vThreads.push_back(std::thread(std::bind(&CTaskPool::doTasks, this)));
    }
    for(int j = 0; j < threadnum; j ++)    {
        vThreads[j].join();
        if(m_verbose  )cout << "thread " << j  << " exit" << endl;
    }
}

void CTaskPool::run(int taskID) {
    if(taskID!=-1 and taskID<vTask.size())
    {
        vTask[taskID]->run();
    }
}

int CTaskPool::getNextTaskID() {
    int nextTask = -1;
    std::lock_guard<std::mutex> guard(nextjob_lock);
    if(nextjob < vTask.size()) {
        if(m_verbose)cout << "Assigning task " << nextjob << endl;
        int percentage = (nextjob+1)*100/vTask.size(), prepercentage=(nextjob)*100/vTask.size();
        if(prepercentage<percentage and m_progress){
            cout << "." <<flush;
            if(percentage%10==0) cout << percentage <<"%";
            if(percentage==100) cout << endl;
        }
        nextjob++;
        nextTask = nextjob-1;
    }
    else {
        if(m_verbose)cout << "No task left!" << endl;
    }
    return nextTask;
}

void CTaskPool::issueTasks(vector<ICThreadTask *> &_vtask, bool verbose, bool progress, int threadnum, string taskname) {
    if(verbose) cout << "Thread num is : " << threadnum << endl;
    CTaskPool ctp(_vtask,verbose, progress);
    ctp.start(threadnum,taskname);
}

XTask::XTask(const XTask &other) :v(other.v),sum(other.sum),f(other.f){}

void XTask::run() {
    sum = 0;
    for(int i = 0; i < v.size(); i ++)    {
        sum += f(v[i]);
        cout << sum << "--" << endl;
    }
}

void threadpool_example() {
    int N = 10;
    int MAXPower = 10;
    vector<double> sumA(MAXPower, 0);
    vector<double> a(N,0);
    iota(a.begin(), a.end(),0);
    cout << a[0] << "  " << a[1] << endl;
    vector<XTask> vTask;  // task pool
    vector<CFptr> vFptr;
    for(int i = 0; i < MAXPower; i ++)    {
        vFptr.emplace_back(CFptr(i));
    }
    for(int i = 0; i < MAXPower; i ++)    {
        double &sum = sumA[i];
        vTask.emplace_back(XTask(a,sum,vFptr[i]));
    }
    vector<ICThreadTask *> vtaskptrs;
    for(int i = 0; i < vTask.size(); i ++)    {
        vtaskptrs.push_back(&vTask[i]);
    }
    CTaskPool ctp(vtaskptrs);
    ctp.start(4, "example task");
    // tasks done

    for(int i = 0; i < MAXPower; i ++)    {
        cout << sumA[i] << endl;
    }
}
