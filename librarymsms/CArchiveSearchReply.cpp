//
// Created by wulong on 4/13/19.
//

#include "CArchiveSearchReply.h"
#include "CAnnotationDB.h"
#include "Util.h"
#include "CPSMAnnotation.h"
#include <iostream>
#include "ICQuery.h"
#include "ProteomicsDataTypes.h"
using namespace std;

// this one support the approxmated distance.
CAnnSpectra::CAnnSpectra(vector<long> &ind, vector<float> &dist,
                         vector<float> &appdist, long queryidx) : m_dim{0},  m_has_spec{false} {
    m_queryidx = queryidx;
    for (int i = 0; i < ind.size(); i++) {
        m_anns.emplace_back(SAnnSpectrum(ind[i], dist[i], appdist[i]));
        m_queueneighbors.emplace( SAnnSpectrum(ind[i], dist[i], appdist[i]) );
    }
}

// approxmated distances are all zeros.
CAnnSpectra::CAnnSpectra(vector<long> &neighborIdx, vector<float> &dist, long queryidx,vector<int> &dpscore) : m_dim{0}, m_has_spec{false}{
    m_queryidx = queryidx;
    if(neighborIdx.size() != dist.size() )    {
        cout << "Error: the size of index and distance should be the same" << endl;
        return;
    }
    else if (neighborIdx.empty())
    {
        //cout << "Warning: the size of neighborIdx is zero " << endl;
        return;
    }
    if(dpscore.size()!=neighborIdx.size()){
        cout << "Error: size of index and dp score are different" << endl;
        dpscore.assign(neighborIdx.size(),0);
    }

    for (int i = 0; i < neighborIdx.size(); i++) {
        m_anns.emplace_back(SAnnSpectrum(neighborIdx[i], dist[i], 0,dpscore[i]));
        m_queueneighbors.emplace(SAnnSpectrum(neighborIdx[i], dist[i], 0,dpscore[i]));
    }

}

// all of the approximate distance are zero values.
CAnnSpectra::CAnnSpectra(int dim, vector<long> &ind, vector<float> &dist, vector<float> &specs, long queryidx) : m_dim{
        dim}, m_has_spec{true} {
    m_queryidx = queryidx;
    for (int i = 0; i < ind.size(); i++) {
        m_anns.emplace_back(SAnnSpectrum(ind[i], dist[i], 0, specs.begin() + i * m_dim,
                specs.begin() + (i + 1) * m_dim)   );
        m_queueneighbors.emplace(SAnnSpectrum(ind[i],dist[i],0, specs.begin() + i * m_dim,
                specs.begin() + (i + 1) * m_dim)  );
    }
}

CAnnSpectra &CAnnSpectra::keepWithMinDP(double mindp) {
    sort(m_anns.begin(), m_anns.end(), [](const SAnnSpectrum &x, const SAnnSpectrum &y){return x.dist < y.dist;});

    int L = m_anns.size() ;
    double max_dist = sqrt(2-2*mindp);
    m_anns.erase(remove_if(m_anns.begin(), m_anns.end(),[&](const SAnnSpectrum &x){return x.dist>max_dist;}), m_anns.end());
    return *this;
}
CAnnSpectra &CAnnSpectra::keepTopN(int n, bool verbose) {
    sort(m_anns.begin(), m_anns.end(), [](const SAnnSpectrum &x, const SAnnSpectrum &y) {
        return x.dist < y.dist;
    });

    if (size() <= n) {
        if(verbose)cout << "size of query result: " << size() << " <= " << n << endl;
    } else {
        m_anns.resize(n);
    }
    // save anns to some where.

    return *this;
}

CAnnSpectra &CAnnSpectra::print(int i) {
    if(i<size())  m_anns[i].print();
    return *this;
}

vector<float> CAnnSpectra::ConcatnateSpec() {
    if(size() ==0) return vector<float>();
    vector<float> v(m_dim * m_anns.size(), 0);
    if (m_has_spec) {
        for (int i = 0; i < m_anns.size(); i++) {
            copy(m_anns[i].spec.begin(), m_anns[i].spec.end(), v.begin() + i * m_dim);
        }
    } else {
        cout << "Do not has spec " << endl;
    }
    return v;
}

vector<long> CAnnSpectra::ConcatenateIndex() {
    if(size() == 0) return vector<long>();
    vector<long> I(size(), 0);
    for (int i = 0; i < size(); i++) {
        I[i] = m_anns[i].idx;
//        if(m_queryidx == I[i])
//            cout << "-------------------------------Find itself ---------------" << endl;
    }
    return I;
}

int CAnnSpectra::size() const {
    return m_anns.size();
}

void CAnnSpectra::print() {
    int width=10;
    print_on_fixed_width(width,"Query","Neighbor","appDist","Dist","\n");

    for(int i = 0; i< size(); i ++)    {
        print_on_fixed_width(width, m_queryidx, m_anns[i].idx, m_anns[i].appdist, m_anns[i].dist, "\n");
    }
    print_on_fixed_width(width,"Query","Neighbor","appDist","Dist","\n");
    int cnt = size();
    while(cnt)    {
        cnt --;
        print_on_fixed_width(width,m_queryidx, m_queueneighbors.top().idx, m_queueneighbors.top().appdist, m_queueneighbors.top().dist,"\n");
        m_queueneighbors.pop();
    }
}

void CAnnSpectra::toLinkStr(string &links) {
    bool empty= true;
    ostringstream oss;
    for(auto it = m_anns.begin(); it!=m_anns.end(); it ++)      {
        if(it->idx == m_queryidx)
            continue;

        if(not empty)  {  oss << ",\n"; }

        tolink(oss, it);
        if(empty) empty = false;
    }
    if(not links.empty()) links += ",\n";
    links += oss.str();
}

vector<SAnnSpectrum>::iterator CAnnSpectra::find(long query_index) {
    auto it = find_if(m_anns.begin(), m_anns.end(),
                      [&query_index](const SAnnSpectrum &x) -> bool { return x.idx == query_index; });
    return it;
}

bool CAnnSpectra::isExist(long query_index) {
    return find(query_index)!=m_anns.end();
}

bool CAnnSpectra::empty() const {
    return m_anns.empty();
}

void CAnnSpectra::tolink(ostringstream &oss, vector<SAnnSpectrum>::iterator &it) const {
    oss << "{"
        << R"("source": ")" << m_queryidx << R"(",)"
        << R"("target": ")" << it->idx << R"(",)"
        << R"("value": )" << it->appdist << ","
        << R"("pvalue": )" << it->pvaluePartial << ","
        << R"("dpscore": )" << it->dotprod << ","
        << R"("realdist": )" << it->dist
        << "}";
}

string CAnnSpectra::createUpdateNeighborSQL() {
    ostringstream oss;
    oss << " NEIGHBOR = '";
    long queryId = getQueryIdx();
    for (int j = 0; j < size(); j ++)   {
        if(j>0) oss << ";";
        oss << m_anns[j].dist << "@" << m_anns[j].idx;
    }
    oss << "'";
    string sql = "update GROUNDTRUTH set " + oss.str() + " where ID=" + to_string(queryId) + ";";
    return sql;
}


void SAnnSpectrum::print() {
    if (!spec.empty()) {
        cout << "idx\tdist\tappdist\tpvalue(all)\tpvalue(partial)\tspec[0]" << endl
             << idx << "\t" << dist << "\t" << appdist << "\t"<<  pvalueAll << "\t" << pvaluePartial << "\t" << spec[0]
             << endl;

    } else {
        cout << "idx\tdist\tappdist\tpvalue(all)\tpvalue(partial)" << endl
             << idx << "\t" << dist << "\t" << appdist<< "\t"<<  pvalueAll << "\t" << pvaluePartial
                << endl;
    }
}

void SAnnSpectrum::setpvaluePartial(float pvaluepartial) {
    pvaluePartial = pvaluepartial;
    pvalues.push_back(pvaluePartial);

}

float SAnnSpectrum::getStdOfPvalue() const {
    using namespace statistic;
    return calcstd<float>(pvalues);
}

SAnnSpectrum::SAnnSpectrum() {
    idx = 0;
    dist = 0;
    appdist = 0;
    pvalueAll = 0;
    pvaluePartial = 10000;
    neighborPvalueAll = false;
    dotprod = 0;
}

SAnnSpectrum::SAnnSpectrum(long idx_, float dist_, float appdist_):SAnnSpectrum() {
    idx=idx_;
    dist = dist_;
    appdist=appdist_;
}

SAnnSpectrum::SAnnSpectrum(long idx_, float dist_, float appdist_, int dotprod_):SAnnSpectrum(idx_,dist_,appdist_) {
    dotprod = dotprod_;
}

SAnnSpectrum::SAnnSpectrum(long idx_, float dist_, float appdist_, vector<float>::iterator start,
                           vector<float>::iterator end) : SAnnSpectrum(idx_,dist_,appdist_){
    spec.assign(start, end);
}

CArchiveSearchReply::CArchiveSearchReply(vector <CAnnSpectra*> &vqr, long queryIndex) : m_vqr(vqr){
    m_query_index = queryIndex;
}

void CArchiveSearchReply::toNodeStr(string &myString, CAnnotationDB &m_AnnotationDB) {
    set<long> nodeIdx;
    getUniqueId(nodeIdx, false);

    vector<CNodeGtInfo> vNodes;
    createNodeList(nodeIdx,m_AnnotationDB, vNodes);

    getNodeStr(vNodes,myString);
}

void CArchiveSearchReply::createNodeList(set<long> &nodeIdx, CAnnotationDB &m_AnnotationDB, vector<CNodeGtInfo> &vNodes) {
    for(auto eachnode : nodeIdx) {
        SPsmAnnotation gtinfo;
        m_AnnotationDB.retrieveGtinfo(eachnode, gtinfo);
        m_AnnotationDB.fixChargeState(gtinfo);

        vNodes.emplace_back(CNodeGtInfo(gtinfo));
    }

    set<string> peptides={"UNKNOWN"};
    map<string, int> peptide_group;
    peptide_group["UNKNOWN"]=1;
    for(auto node : vNodes)    {
        if(peptides.count(node.getPepSeq())>0) {

        }
        else  {
            peptides.insert(node.getPepSeq());
            int t = peptide_group.size()+1;
            peptide_group[node.getPepSeq()] = t;
        }
    }

    for(auto & vNode : vNodes)    {
        vNode.setGroup(peptide_group[vNode.getPepSeq()]);
    }
}

void CArchiveSearchReply::toLinkStr(string &links) {
    double dist_max = 1.5;
    for(auto each : m_vqr)        {
        each->toLinkStr(links);
    }
}

void CArchiveSearchReply::getUniqueId(set<long> &nodeIdx, bool verbose) {
    nodeIdx.insert(m_query_index);
    for(auto eachQuery : m_vqr)  {
        if(m_query_index == eachQuery->getQueryIdx() and verbose)            {
            cout << "find query_ID:" << eachQuery->getQueryIdx() << endl;
        }
        nodeIdx.insert(eachQuery->getQueryIdx());
        for(const auto& eachneighbor : eachQuery->m_anns)  {
            nodeIdx.insert(eachneighbor.idx);
        }
    }
    // a set is supposed to be unique,
}

void CArchiveSearchReply::getNodeStr(vector<CNodeGtInfo> &vNodes, string &myNodes) {
    for(int i = 0; i < vNodes.size(); i ++)  {
        if(i>0) myNodes += ",\n";
        myNodes += "{" + vNodes[i].toNodeStr()+"}";
    }
//    cout << "debug:get node str" << endl;
}

void CArchiveSearchReply::tojsonstring(string &jsonstring, CAnnotationDB &m_AnnotationDB) {
    string nodes, links;
    toNodeStr(nodes, m_AnnotationDB);
    toLinkStr(links);
    stringstream oss;
    oss << "{\n\"nodes\": [\n" << nodes << "],\n" << endl;
    oss << "\"links\": [\n" << links << "]\n}\n" << endl;
    jsonstring = oss.str();
}

bool CompareEntryDist::operator()(const SAnnSpectrum &a, const SAnnSpectrum &b) {return a.dist > b.dist;}

void CArxivSearchResult::copyto(vector<CAnnSpectra *> &vqrs, int offset) {
    copy(m_queryresults.begin(), m_queryresults.end(), vqrs.begin() + offset);
}

// move results to vqrs.
void CArxivSearchResult::moveTo(vector<CAnnSpectra *> &vqrs, int offset) {
    copyto(vqrs, offset);
    for (int k = 0; k < size(); k++) {
        set(k, nullptr, false);
    }
}

void CArxivSearchResult::release(int i) {
    if(i<0 or i>size()) return;
    set(i, nullptr);
}

vector<long> CArxivSearchResult::concatenateIndex(int i) {
    return get(i)->ConcatenateIndex();
}

void CArxivSearchResult::keepTopN(int topN, bool verbose) {
    for(int i = 0; i < size(); i ++)  {
        keepTopN(i, topN, verbose);
    }
}

void CArxivSearchResult::keepTopN(int i, int topN, bool verbose) {
    get(i)->keepTopN(topN, verbose);
}

void CArxivSearchResult::push_back(CAnnSpectra *p) {
    m_queryresults.push_back(p);
}

void CArxivSearchResult::set(int i, CAnnSpectra *p, bool releasefirst) {
    if(releasefirst and get(i)!=nullptr) delete get(i);
    m_queryresults[i] = p;
}

CArxivSearchResult::~CArxivSearchResult() {
    for(int i = 0; i < m_queryresults.size(); i ++)    {
        release(i);
    }
}

CAnnSpectra *CArxivSearchResult::get(int i) {
    return m_queryresults.at(i);
}

int CArxivSearchResult::size() {return m_queryresults.size();}

void CArxivSearchResult::init(ICQuery &query, bool verbose, vector<vector<long>> &allRetIdx,
                              vector<vector<float>> &accDist,vector<vector<int>> &dpscores) {
    int querySize = query.size();
    int zeroNeighborCounts = 0;
    for (int i = 0; i < querySize; i++) {
        if (verbose) query.getMzSpec(i).display();

        // we are calling the second one. but only the third on support approximate distance.
        CAnnSpectra *p = new CAnnSpectra(allRetIdx[i], accDist[i], query.getQueryIndex(i),dpscores[i]);
        push_back(p);
        if(p!= nullptr and p->empty())
        {
            zeroNeighborCounts ++;
        }
    }
    if(verbose){
    cout << "Zero neighbor query counts : " << zeroNeighborCounts << "out of total " << querySize<< endl;

    }
}

void CArxivSearchResult::keepWithMinDP(double mindp) {
    for(int i = 0; i < size(); i ++)  {
        keepWithMinDP(i, mindp);
    }
}

void CArxivSearchResult::keepWithMinDP(int i, double mindp) {
    get(i)->keepWithMinDP(mindp);
}

CNodeGtInfo::CNodeGtInfo(const SPsmAnnotation& gtinfo) {
    m_gtinfo = make_shared<SPsmAnnotation>(gtinfo);
    m_group = -1;
}

string CNodeGtInfo::toNodeStr() {
    ostringstream oss;
    oss << R"("id": ")" << m_gtinfo->idx
        << R"(","group":)" << m_group << ",";
    m_gtinfo->toOstringStreamNoId(oss);
    string currentnode = oss.str();
    return currentnode;
}

string CNodeGtInfo::getPepSeq() {return m_gtinfo->peptideseq;}

void CNodeGtInfo::setGroup(int group) {m_group=group;}
