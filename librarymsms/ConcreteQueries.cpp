//
// Created by wulong on 6/27/19.
//

#include "ConcreteQueries.h"
#include "ICMzFile.h"
#include "PeakList.h"
#include "ProteomicsDataTypes.h"


CMzQuery::CMzQuery(ICMzFile &Mzs, int dim, vector<long> &idx, bool useflankingbins) : m_dim(dim), m_Mzs(Mzs), m_idx(idx)
{
    m_queries = nullptr;
    m_useflankingbins = useflankingbins;
}

CMzQuery::~CMzQuery()
{
    if (m_queries)
    {
        delete [] m_queries;
        m_queries = nullptr;

    }
}

float * CMzQuery::get()
{
    //cout << "get pointer first" << endl;
    if (m_queries == nullptr)
    {
        m_queries = m_Mzs.buildNormalizedQueries(m_dim, m_idx, m_useflankingbins, m_idx.size());
    }
    return m_queries;
}


void CMzQuery::get(vector<float> &queries)
{
    float *p = get();
    if (get() != nullptr)
    {
        queries.assign(m_queries, m_queries + dim()*size());
    }
}

long CMzQuery::dim()
{
    return m_dim;
}

int CMzQuery::size()
{
    return m_idx.size();
}

uint16_t *CMzQuery::getPtrUint16(int i) {
    return m_Mzs.getMzSpec(m_idx[i]).getPeakPtr();
}

long CMzQuery::getQueryIndex(int i) {
    return m_idx[i];
}


CSpecQuery::CSpecQuery(vector<uint16_t> &query, int dim, bool useflankingbins) : m_dim{dim}, m_useflankingbins{useflankingbins}, m_queryspec{query}
{
    m_queries = nullptr;
}

CSpecQuery::~CSpecQuery() {
    if (m_queries)
    {
        delete [] m_queries;
        m_queries = nullptr;

    }
}

float *CSpecQuery::get() {
    vector<float> v;
    int peaknum = m_queryspec.size();
    const int PEAKNUM_PER_SEPC = 50;

    if(m_queries == nullptr)
    {
        // build m_queries
        vector<double> mz(peaknum, 0), intensity(peaknum, 0);
        const double MAX_PEAK_NUM = UINT16_MAX;
        const double MAX_PEAK_MZ = 2000;
        for (int i = 0; i < m_queryspec.size(); i++) {
            if (m_queryspec[i] == 0) {
                mz.resize(i);
                intensity.resize(i);
                break;
            }
            mz[i] = m_queryspec[i];
            mz[i] /= MAX_PEAK_NUM;
            mz[i] *= MAX_PEAK_MZ;
            intensity[i] = PEAKNUM_PER_SEPC - i;
        }

        BinningPeakList bpl(mz, intensity, m_useflankingbins);

        getFloatVecPaddingZeros(&bpl, dim(), false, v);
        m_queries = new float[dim()]; // todo
        std::copy(v.begin(),v.end(), m_queries);
    }
    return m_queries;
}

void CSpecQuery::get(vector<float> &queries) {
    float *p = get();
    if(get() != nullptr)
    {
        queries.assign(m_queries,m_queries + dim() * size());
    }
}

long CSpecQuery::dim() {return m_dim;}

int CSpecQuery::size() {return 1;}

uint16_t *CSpecQuery::getPtrUint16(int i) {
    return m_queryspec.data();
}

long CSpecQuery::getQueryIndex(int i) {
    cout << "[Info] Spectrum out of archive got query id -1" << endl;
//    cout << "Error: SpecQuery Object always return zero..." << endl;
    return -1;
}

// Note: query id can be -1
// When first query id is -1, we search with the spectrum in query pointer!
ICQuery *createICQuery(vector<uint16_t> *query, vector<long> *idx1R, bool useflankingbins, int dim, ICMzFile &csa) {
    if (idx1R != nullptr and idx1R->size() == 1 and idx1R->at(0) == -1 and query != nullptr) {
        return new CSpecQuery(*query, dim, useflankingbins);
    } else if (idx1R != nullptr) {
        return new CMzQuery(csa, dim, *idx1R, useflankingbins);
    } else {
        cout << "Invalid parameter for initialize ICQuery object" << endl;
        cout << "try (-1, query, nullptr)   or (queryindex, nullptr, idx1R)" << endl;
        return nullptr;
    }
}
