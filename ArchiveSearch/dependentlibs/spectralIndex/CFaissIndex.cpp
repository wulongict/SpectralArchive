//
// Created by wulong on 6/19/19.
//



#include "spdlog/spdlog.h"

#include <faiss/AutoTune.h>

#include "CFaissIndex.h"
#include <faiss/index_io.h>

#include <iostream>
#include <faiss/AutoTune.h>

#if FAISS_VERSION_MAJOR == 1 and FAISS_VERSION_MINOR >=6
// this header is available in faiss 1.6.4
#include <faiss/index_factory.h>
#include <faiss/IndexIVFPQ.h>
#include <fstream>

#endif

#ifdef __CUDA__
// #pragma message("__FILE__: __CUDA__: XXX using CUDA")
//#warning "C Preprocessor got here!"

#include <faiss/gpu/GpuAutoTune.h>
#include <faiss/gpu/StandardGpuResources.h> // standard GPU Resources
#include <faiss/gpu/GpuCloner.h>

// --------------------------------------------
#include <faiss/gpu/GpuIndexIVFPQ.h>
#include <faiss/index_io.h>

#else
// #pragma message("__FILE__: __CUDA__: XXX NO CUDA")
#warning NO CUDA
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
    setPtr(faiss::read_index(filename.c_str()));
}

// write the bucket size to a file, measure the imbalance as entropy.
void CFaissIndexWrapper::display() {
    int number_of_lists = ((faiss::IndexIVFPQ*)m_index)->nlist;
    vector<long> size_of_list(number_of_lists,0);


    double sum = 0;
    for(int i = 0; i < number_of_lists; i ++){
        long sizelist = ((faiss::IndexIVFPQ*)m_index)->get_list_size(i);
        size_of_list[i] = sizelist;
        sum += sizelist;
    }
    double entropy = 0;
    for(int i = 0; i < number_of_lists; i ++){
        if(size_of_list[i]>0){
            double p = size_of_list[i]/sum;
            entropy -= p* log2(p); 
        }
    }
    int sacodesize = ((faiss::IndexIVFPQ*)m_index)->sa_code_size(); 

    string outputfile = getfilename() + "_nlist.txt";
    ofstream fout(outputfile.c_str(), ios::out);
    fout << "Entropy\t" << entropy << endl;
    fout << "code size\t" << sacodesize << endl;
    for(int i = 0; i < number_of_lists; i ++){
        fout << i << "\t" << size_of_list[i] << endl;
    }
    fout.close();

    

#if FAISS_VERSION_MAJOR == 1 and FAISS_VERSION_MINOR <6
    getPtr()->display();
#endif
}

// save index to file.
void CFaissIndexWrapper::write(string filename) {
    faiss::write_index(getPtr(), filename.c_str());
}

// get number of vectors in the index
long CFaissIndexWrapper::total() {
    return getPtr()->ntotal;
}

void CFaissIndexWrapper::search(int num, float *queries, int ret_num, vector<float> &dist,
                                vector<long> &ind) {
    
    faiss::Index * ptr = getPtr();
    ptr->search(num, queries, ret_num, dist.data(), ind.data());
}

void CFaissIndexWrapper::setnProb(int nprobe) {
    if (m_isCPU) {
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
    m_index = newptr;
    m_index->verbose = true;
}



void CFaissIndexWrapper::toCPU() {
    if(not m_isCPU){
        #ifdef __CUDA__
        cout << "[Info] Moving index from GPU to CPU" << endl;
        m_index = faiss::gpu::index_gpu_to_cpu(m_index);
        m_isCPU = true;
        #endif
    }

}

// get number of gpus with Faiss library 
int getNumGPUs(){
    int num_gpus = 0;
#ifdef __CUDA__
    num_gpus = faiss::gpu::getNumDevices();
#endif
    return num_gpus;
}


void CFaissIndexWrapper::createEmptyIndex(int dim, string indexstr) {
    setPtr(faiss::index_factory(dim, indexstr.c_str()));
}

#include <map>
void CFaissIndexWrapper::toMultipleGPUs(vector<int> &gpu_idx) {

    if(m_gpu_idx.empty()){
        m_gpu_idx = gpu_idx;
        
    }
    if (not m_isCPU){
        // already using GPU. no change.
    }
    else{

#ifdef __CUDA__
        cout << "[Info] Moving index from CPU to MultipleGPU" << endl;
        int ngpus = faiss::gpu::getNumDevices();
        printf("Number of GPUs found: %d\n", ngpus);
        int num_valid_gpu_id = 0;
        std::map<int, bool> id_map;
        for (auto x : m_gpu_idx)
        {
            if (x < ngpus)
            {
                num_valid_gpu_id += 1;
                id_map[x] = true;
            }
        }
        if (num_valid_gpu_id == 0)
        {
            cout << "[Info] all the gpu ids are invalid, will try first GPU. " << endl;
            m_gpu_idx.erase(gpu_idx.begin(), gpu_idx.end());
            m_gpu_idx.push_back(0);
            id_map[0] = true;
        }
        std::vector<faiss::gpu::GpuResourcesProvider *> res;
        std::vector<int> devs;
        for (int i = 0; i < ngpus; i++)
        {
            if (id_map.count(i))
            {
                // cout << "adding gpu " << i << endl;
                res.push_back(new faiss::gpu::StandardGpuResources);
                devs.push_back(i);
            }
        }
        m_index = faiss::gpu::index_cpu_to_gpu_multiple(res, devs, m_index);


        m_isCPU = false;
#endif
    }
}

void CFaissIndexWrapper::toGPU() {
    if(m_gpu_idx.empty()){
        m_gpu_idx.push_back(0);
        }
    toMultipleGPUs(m_gpu_idx);
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

void CFaissIndexWrapper::removeIds(vector<long> &idx) {
    faiss::IDSelectorRange sel(idx.front(), idx.back()+1);
    cout << "removing index of size " << idx.size() << endl;
    bool usingGPU = not m_isCPU;
    if(usingGPU){
        toCPU();
    }
    int n = m_index->remove_ids(sel);
    if(usingGPU){
        toGPU();
    }
    cout << n << " vectors removed from index" << endl;

}
