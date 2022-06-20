//
// Created by wulong on 8/4/18.
//

#ifndef MYTOOL_LIBSEARCHENGINE_H
#define MYTOOL_LIBSEARCHENGINE_H

#include <string>
using namespace std;



//
//// maybe used someday
//class LibSearchEngine {
//    faiss::Index *m_index;
//    int m_dim = 4096;
//    // annotation of each spectrum
//
//public:
//    LibSearchEngine(string libraryname)
//    {
//        m_dim = 4096;
//        long psmnum = 0;
//        DataFile df(libraryname);
//        float *vec = df.toFloatVector(m_dim, psmnum, false, false, 0.05);
//        m_index = faiss::index_factory(m_dim, "IVF256,PQ64");
//        m_index->train(psmnum, vec); // train and added
//        m_index->add(psmnum, vec);
//        delete [] vec;
//    }
//    void QueryTopN(int topn, float *query, vector<float> &dist, vector<faiss::Index::idx_t > &ind)
//    {
//        // dim is defined
//
//        m_index->search(1,query,topn,dist.data(),ind.data());
//    }
//    void MultiQueryTopN(int numQuery, int topn, float *query, vector<float> &dist, vector<faiss::Index::idx_t > &ind)
//    {
//        m_index->search(numQuery,query,topn, dist.data(), ind.data());
//    }
//    void AppendDataToIndex(string filename)
//    {
//        DataFile df(filename);
//        long psmnum = 0;
//        float *vec = df.toFloatVector(m_dim, psmnum, false, false, 0.05);
//        m_index->add(psmnum,vec);
//        delete [] vec;
//    }
//    // batch addition is safer
//    void AppendDataToIndexBatch(string filename, long batchsize=1000)
//    {
//        DataFile df(filename);
//        long psmnum = 0;
//        float *vec = df.toFloatVector(m_dim, psmnum, false, false, 0.05);
//
//        for(long i = 0; i*batchsize < psmnum; i ++)
//        {
//            long nsize = batchsize;
//            if (i*batchsize + batchsize > psmnum)
//            {
//                nsize = psmnum - i*batchsize;
//            }
//
//            m_index->add(nsize, vec+i*batchsize);
//            cout << "batch i=" << i+1 << "\t" << nsize << " vectors added to index" << endl;
//        }
//        m_index->add(psmnum,vec);
//        delete [] vec;
//    }
//
//    void SaveIndex(string indexfilename){
//        faiss::write_index(m_index, indexfilename.c_str());
//        cout << "index saved into file: " << indexfilename << endl;
//    }
//
//
//};


#endif //MYTOOL_LIBSEARCHENGINE_H
