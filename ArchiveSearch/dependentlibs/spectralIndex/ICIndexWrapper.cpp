//
// Created by wulong on 2/18/19.
//


#include "CMyIndex.h"
#include "CFaissIndex.h"
#include "CKMeans.h"

void ICIndexWrapper::read() {
    read(getfilename());
}

void ICIndexWrapper::write() {
    write(getfilename());
}

void ICIndexWrapper::setfilename(string filename) {
    m_filename = filename;
}

string ICIndexWrapper::getfilename() {
    return m_filename;
}

ICIndexWrapper::~ICIndexWrapper() {}

shared_ptr<ICIndexWrapper> IndexFactory(bool myIndex, shared_ptr<CPQParam> option) {
    if (myIndex) {
//            CPQParam option;
        return make_shared<CMyIndex>(option);
    } else {
        return make_shared<CFaissIndexWrapper>();
    }
}
