//
// Created by wulong on 6/19/19.
//

#ifndef MYTOOL_GPURESOURCEOBJ_H
#define MYTOOL_GPURESOURCEOBJ_H
#ifdef __CUDA__
#pragma message("__FILE__: __CUDA__: XXX using CUDA")
//#warning "C Preprocessor got here!"

#include <faiss/gpu/GpuAutoTune.h>
#include <faiss/gpu/StandardGpuResources.h> // standard GPU Resources
#include <faiss/gpu/GpuCloner.h>

// --------------------------------------------
#include <faiss/gpu/GpuIndexIVFPQ.h>
#include <faiss/index_io.h>


class GpuResourceObj
{
    static GpuResourceObj * m_gpuptr ;
    faiss::gpu::StandardGpuResources m_res;
    faiss::gpu::GpuClonerOptions m_gco;
public:
    static  GpuResourceObj * getInstance();
    faiss::Index * cpu_to_gpu(faiss::Index * cpuIndexPtr);
    faiss::Index * gpu_to_cpu(faiss::Index *gpuIndexPtr);
private:
    GpuResourceObj() {}

};
#else
#pragma message("__FILE__: __CUDA__: XXX NO CUDA")
#warning NO CUDA
#endif


#endif //MYTOOL_GPURESOURCEOBJ_H
