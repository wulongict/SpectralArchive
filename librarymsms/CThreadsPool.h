//
// Created by wulong on 4/24/19.
//

#ifndef MYTOOL_CTHREADSPOOL_H
#define MYTOOL_CTHREADSPOOL_H

#include <vector>
#include <mutex>
#include <iostream>
#include <cmath>
#include <memory>
using namespace std;

class ICThreadTask{
public:
    virtual ~ICThreadTask(){}
    virtual void run() = 0;
};

class CTaskPool{
    bool m_verbose;
    vector<ICThreadTask*> &vTask;
    int nextjob;
    mutex nextjob_lock;
    bool m_progress;
public:
    CTaskPool(vector<ICThreadTask*> &_vtask, bool verbose=true,bool progress=true);
    void start(int threadnum, string taskname);
    static void issueTasks(vector<ICThreadTask*> &_vtask, bool verbose,bool progress, int threadnum, string taskname );
private:
    void doTasks();
    int getNextTaskID();
    void run(int taskID);
};

// Here is a example to use the thread pool tool

class CFptr{
    double p;
public:
    CFptr(double _p):p(_p){}
    double operator()(double x){
        cout << "this is fp with power: " << p << endl;return pow(x,p);
    }
};


class XTask: public ICThreadTask{
    vector<double> &v;
    double &sum;

    CFptr &f;
public:
    XTask(vector<double> &_v, double &_sum, CFptr &_f):v(_v),sum(_sum),f(_f){}
    XTask(const XTask &other);
    void run();
};


void threadpool_example();

#endif //MYTOOL_CTHREADSPOOL_H
