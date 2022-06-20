//
// Created by wulong on 6/27/19.
//

#ifndef MYTOOL_CONCRETEQUERIES_H
#define MYTOOL_CONCRETEQUERIES_H

#include "ICQuery.h"

class ICMzFile;
class ICQuery;

class CMzQuery : public ICQuery{
    int m_dim;
    ICMzFile &m_Mzs;
    vector<long> &m_idx;
    float *m_queries;
    bool m_useflankingbins;
public:
    CMzQuery(ICMzFile &Mzs, int dim, vector<long> &idx, bool useflankingbins);
    ~CMzQuery() override;
    void get(vector<float> &queries) override;
    float * get() override;
    long dim() override;
    uint16_t * getPtrUint16(int i) override;
    int size() override;
    long getQueryIndex(int i) override;
};

class CSpecQuery: public ICQuery
{
    int m_dim;
    float *m_queries;
    bool m_useflankingbins;
    vector<uint16_t> &m_queryspec;
public:
    CSpecQuery(vector<uint16_t> &query, int dim, bool useflankingbins);
    ~CSpecQuery() override;
    float *get() override;
    void get(vector<float> &queries) override;
    long dim() override;
    int size() override;
    uint16_t *getPtrUint16(int i ) override;
    long getQueryIndex(int i) override;
};


// factory method
ICQuery * createICQuery(vector<uint16_t> *query,vector<long> *idx1R, bool useflankingbins, int dim, ICMzFile &csa);


#endif //MYTOOL_CONCRETEQUERIES_H
