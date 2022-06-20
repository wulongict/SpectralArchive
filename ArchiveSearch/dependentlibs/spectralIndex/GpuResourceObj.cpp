//
// Created by wulong on 6/19/19.
//
#ifdef __CUDA__

#include <iostream>
#include "GpuResourceObj.h"

//---------------------
#include <faiss/gpu/GpuIndexIVFPQ.h>
#include <faiss/gpu/StandardGpuResources.h>

// #include <faiss/gpu/GpuAutoTune.h>
#include <faiss/gpu/GpuCloner.h>
#include <faiss/index_io.h>

using namespace std;

GpuResourceObj *GpuResourceObj::getInstance() {
    if (m_gpuptr == nullptr) {
        m_gpuptr = new GpuResourceObj();
        long memSize = 4 * 1024L * 1024L * 1024L;
//        memSize = 0;
        m_gpuptr->m_res.setTempMemory(memSize); // 4GB
    }
    return m_gpuptr;
}

GpuResourceObj *GpuResourceObj::m_gpuptr = nullptr;

faiss::Index *GpuResourceObj::cpu_to_gpu(faiss::Index *cpuIndexPtr) {
    
    return faiss::gpu::index_cpu_to_gpu(&m_res, 0, cpuIndexPtr, &m_gco);
}

faiss::Index *GpuResourceObj::gpu_to_cpu(faiss::Index *gpuIndexPtr) {
    return faiss::gpu::index_gpu_to_cpu(gpuIndexPtr);
}
#else
#pragma message("__FILE__: __CUDA__: XXX NO CUDA")
#warning NO CUDA

#endif
