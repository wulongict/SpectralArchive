//
// Created by wulong on 6/19/19.
//



#include "spdlog/spdlog.h"

#include <faiss/AutoTune.h>
#include "CFaissIndex.h"
#include <faiss/index_io.h>
#include "GpuResourceObj.h"
#include <iostream>
#include <faiss/AutoTune.h>
// #include <faiss/gpu/GpuAutoTune.h>
#if FAISS_VERSION_MAJOR == 1 and FAISS_VERSION_MINOR >=6
// this header is available in faiss 1.6.4
#include <faiss/index_factory.h>
#endif
// GPU CUDA

faiss::Index *CFaissIndexWrapper::getPtr() {
    return m_index;
}

CFaissIndexWrapper::CFaissIndexWrapper() {
    m_index = nullptr;
}

void CFaissIndexWrapper::read(string filename) {
    m_isCPU = true;
//    spdlog::get("A")->info("Loading index from disk... ");
    setPtr(faiss::read_index(filename.c_str()));
}

void CFaissIndexWrapper::display() {
#if FAISS_VERSION_MAJOR == 1 and FAISS_VERSION_MINOR <6
    getPtr()->display();
#endif
}

void CFaissIndexWrapper::write(string filename) {
    faiss::write_index(getPtr(), filename.c_str());
}

long CFaissIndexWrapper::total() {
    return getPtr()->ntotal;
}

void CFaissIndexWrapper::search(int num, float *queries, int ret_num, vector<float> &dist,
                                vector<long> &ind) {
//    cout << "Searching with pointer for " << num << " queries " << ret_num << endl;
    getPtr()->search(num, queries, ret_num, dist.data(), ind.data());
//    cout << "Query done!" << endl;
}

void CFaissIndexWrapper::setnProb(int nprobe) {
    if (m_isCPU) {
        //cout << "Setting nprobe as " << nprobe << endl;
        faiss::ParameterSpace().set_index_parameter(getPtr(), "nprobe", nprobe);

    } else {
        if (nprobe > 64) {
            nprobe = 64;
            cout << "Using GPU: #Buckets should be less than 64: nprobe= " << nprobe << endl;
        }
#ifdef __CUDA__
        faiss::gpu::GpuParameterSpace().set_index_parameter(getPtr(), "nprobe", nprobe);
        // faiss::GpuParameterSpace().set_index_parameter(getPtr(), "nprobe", nprobe);
        
#endif
    }
}

bool CFaissIndexWrapper::empty() {
    return getPtr() == nullptr;
}

int CFaissIndexWrapper::dim() {
    return getPtr()->d;
}

bool CFaissIndexWrapper::istrained() {
    return getPtr()->is_trained;
}

void CFaissIndexWrapper::setPtr(faiss::Index *newptr) {
//    std::cout << "pointer: " << m_index << " --> " << newptr << std::endl;
    m_index = newptr;
}

void CFaissIndexWrapper::backup() {
    m_tmp = m_index;
}

void CFaissIndexWrapper::toCPU() {
    delete m_index;
    m_index = m_tmp;
}

void CFaissIndexWrapper::createEmptyIndex(int dim, string indexstr) {
    cout << "Creating empty index ..." << endl;
    setPtr(faiss::index_factory(dim, indexstr.c_str()));
}

void CFaissIndexWrapper::toGPU() {
#ifdef __CUDA__
    cout << "[Info] Moving index from CPU to GPU" << endl;
    backup();
    GpuResourceObj *gpuRes = GpuResourceObj::getInstance();
    setPtr(gpuRes->cpu_to_gpu(getPtr()));
    m_isCPU = false;
#endif
}

void CFaissIndexWrapper::train(int batchSize, float *vecBatch) {
    m_index->train(batchSize, vecBatch);
}

CFaissIndexWrapper::~CFaissIndexWrapper() {
    if (not empty()) {
        delete getPtr();
        setPtr(nullptr);
    }
}

void CFaissIndexWrapper::add(long newspecnum, float *vec) {
    getPtr()->add(newspecnum, vec);
}