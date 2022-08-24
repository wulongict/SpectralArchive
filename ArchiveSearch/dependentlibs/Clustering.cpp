//
// Created by wulong on 9/2/20.
//

#include "Clustering.h"
#include <iostream>
#include "../../librarymsms/CAnnotationDB.h"
#include "../../librarymsms/Util.h"
#include "../../librarymsms/CPSMAnnotation.h"
#include "../../librarymsms/CMzFileReader.h"

using namespace std;

ClusterNode::ClusterNode() : isVisited(false), id(-1), clusterId(-1), nodedepth(0), peptide(""), isSignificant(false) {}

ClusterNode::ClusterNode(long id) : id(id), isVisited(false), clusterId(-1), nodedepth(0), peptide(""),
                                    isSignificant(false) {}

std::string ClusterNode::str() {
    return to_string("node:\t", "\t", id, clusterId, nodedepth, isSignificant, isVisited);
}

void Cluster::addNode(ClusterNode *nodePtr) {
    nodePtr->isVisited = true;
    m_nodeList.push_back(nodePtr);
    nodePtr->clusterId = cluster_id;
}

string Cluster::str(int topN) {
    ostringstream oss;
    oss << "size: " << m_nodeList.size() << " \nnodes: ";
    for (int i = 0; i < topN and i < m_nodeList.size(); i++) {
        oss << m_nodeList[i]->id << " ";
    }
    oss << "..." << endl;
    return oss.str();
}

void Cluster::print(int topN) {
    cout << "size: " << m_nodeList.size() << " \nnodes: ";
    for (int i = 0; i < topN and i < m_nodeList.size(); i++) {
        cout << m_nodeList[i]->id << " ";
    }
    cout << "..." << endl;
}

ostream &operator<<(ostream &os, const Cluster &node) {
    // Cluster ID: neighborSize, neighbor1, neighbor2, neighbor3 ...
    // ClusterSize: node1, node2, node3...
    os << "M\t" << node.cluster_id << "\t" << node.m_neighborClusters.size();
    for (auto &item: node.m_neighborClusters) {
        os << "\t" << item->cluster_id;
    }
    os << endl;

    os << "D\t" << node.cluster_id << "\t" << node.m_pepCounts.size() << "\t" << node.getDarkNodesNum()
       << "\t" << node.size() << "\t" << node.getDarkNodesNum() * 1.0 / node.size() << endl;

    os << "E\t" << node.cluster_id << "\t" << node.m_entropy << "\t" << node.m_pepCounts.size();
    for (auto &item: node.m_pepCounts) {
        os << "\t" << item.first << "\t" << item.second;
    }
    os << endl;

    os << "N\t" << node.cluster_id << "\t" << node.m_nodeList.size();
    for (int j = 0; j < node.m_nodeList.size(); j++) {
        os << "\t" << node.m_nodeList[j]->id;
    }
    return os;
}

int Cluster::size() const {
    return m_nodeList.size();
}

void Cluster::addNeighborClusters(std::unordered_set<Cluster *> neighborClusters) {
    for (auto &each: neighborClusters) {
        if (each->cluster_id != cluster_id) {
            m_neighborClusters.insert(each);
        }
    }
//    m_neighborClusters.insert(neighborClusters.begin(), neighborClusters.end());
}

int Cluster::getDarkNodesNum() const {
    int count = 0;
    for (auto &x: m_nodeList) {
        if (not x->isSignificant) count++;
    }
    return count;
}

void Cluster::mergeInto(Cluster *target) {
    for (auto &eachnode : m_nodeList) {
        target->addNode(eachnode);
    }
//    target->addNode()
//    target->m_nodeList.insert(target->m_nodeList.end(),m_nodeList.begin(), m_nodeList.end());
    // transfer all the node to a different cluster now.

    m_nodeList.resize(0);
    target->addNeighborClusters(m_neighborClusters);
    m_neighborClusters.clear();
}


Clusters::Clusters(long nodeNum, shared_ptr<CMzFileReader> csa) {
    m_csa = csa;
    m_nodeNum = nodeNum;
    m_nodes.assign(nodeNum, ClusterNode());
    for (long i = 0; i < nodeNum; i++) {
        m_nodes[i].id = i;
    }
    int minPeakNum = 6;
    m_nodes.erase(std::remove_if(m_nodes.begin(), m_nodes.end(),
                                 [&](const ClusterNode &x) {
                                     return csa->getPeakNum(x.id) <= minPeakNum;
                                 }), m_nodes.end());
    cout << "Nodes reduced from " << m_nodeNum << " to " << m_nodes.size() << " with min peak num: " << minPeakNum
         << endl;
    m_nodeNum = m_nodes.size();
    for (long i = 0; i < getNodeNum(); i++) {
        m_nodeId2NodePtr[m_nodes[i].id] = &m_nodes[i];
    }
}

// A tricky question is the distance is not symmetric, solved by merging
// merging case
// I. The 90% of the identified nodes are from same group
// II. There is only one neighbor cluster.

// todo:
// two things: mergePureDominantCluster; always merge;
void
Clusters::makeClusters(CAnnotationDB &mdb, double maxDist, long start, long end, bool mergeNode, int maxSearchDepth,
                       bool usePepWithPTM, bool alwaysMerge) {
    if (end < 0 or end > getNodeNum()) end = getNodeNum();
    if (start < 0 or start > end) start = 0;
    Progress ps(end - start, "task clustering");
    for (long i = start; i < getNodeNum() and i < end; i++) {
        if (isNodeVisited(i)) {
            continue;
        }
        setNodeDepth(i,0);
        try {
            unordered_set<ClusterNode *> neighborNodes = findNeighboringNodes(&m_nodes[i], mdb, maxDist, usePepWithPTM);
            // calculate entropy of neighbor nodes.
            CountFrequency clusterIdFreq;
            for (auto &eachNode: neighborNodes) {
                if (eachNode->isVisited) clusterIdFreq.add_data(eachNode->clusterId);// what if their cluster id is -1
            }

            unordered_set<Cluster *> neighborClusters = findNeighboringClusters(neighborNodes);

            // The distribution of cluster id in current set of neighbors. calculate entropy of neighbor nodes.
            removeVisited(neighborNodes);

            Cluster *newCluster = nullptr;
            // merge node and always merge can exist at the same time.
            if (mergeNode and not alwaysMerge) {
                double majorityProp = 0.9;
                // find dominant cluster!
                auto it = find_if(neighborClusters.begin(), neighborClusters.end(),
                                  [&](const Cluster *x) {
                                      return clusterIdFreq.getFreq(x->cluster_id) >
                                             majorityProp * clusterIdFreq.totalInsertion();
                                  });
//            more than one neighbor and there is a dominant one
                if (neighborClusters.size() > 1 and it != neighborClusters.end()) {
                    newCluster = *it;
                } else if (neighborClusters.size() == 1) {
                    // exact one, and it must be dominant one.
                    newCluster = *neighborClusters.begin();
                } else {
                    // no dominant one, create a new one.
                    newCluster = createCluster();
                    newCluster->addNeighborClusters(neighborClusters);
                }
            } else if (alwaysMerge) {
                // no neighbor, create new one
                // something else
                if (neighborClusters.empty()) {
                    newCluster = createCluster();
//                newCluster->addNeighborClusters(neighborClusters);
                } else if (neighborClusters.size() == 1) {
                    // merge to the first one, and delete the others.
                    newCluster = *neighborClusters.begin();
                } else {
                    // more than one neighbor, merge it to the min id cluster.
                    // one neighbor, merge to it
                    newCluster = mergeClusters(neighborClusters); // more than one cluster
                }
            } else {
                newCluster = createCluster();
                newCluster->addNeighborClusters(neighborClusters);
            }

            newCluster->addNode(&m_nodes[i]);
            // not merge, and not always merge. always create new cluster

            for (auto &eachitem: neighborNodes) {
                eachitem->nodedepth = m_nodes[i].nodedepth + 1;
            }
            int addCounts = BFS(neighborNodes, *newCluster, mdb, maxDist, maxSearchDepth, usePepWithPTM, alwaysMerge);
            ps.increase(addCounts + 1);
//            }
//            }
        } catch (exception &e) {
            cout << "i=" << i << endl << e.what() << endl;
        } catch (...) {
            cout << "i=" << i << endl;
        }
    }
    cout << "Clustering done" << endl;
    entropyOfEachCluster();
    getSummary();
}


#include <algorithm>
//// if you have two type of object
//// A1, A2, A3, A4,
//// B1, B2, B3, B4
////
//// and you would like to conversion between each pair of them, how can you do this!
//// you need 4 by 4 = 16 convesion functions.
//class ScoreConversionTools{
//public:
//    double m_score;
//    ScoreConversionTools(double score): m_score(score){
//
//    }
//    ScoreConversionTools(): ScoreConversionTools(0){
//
//    }
//
//    template<typename T>
//    T as(){
//        return T(m_score);
//    }
//    struct DPscore{
//        double score;
//        DPscore(double m_score){
//            score = m_score;
//        }
//        double toDist(){
//            //
//            // d^2 = x^2 + y^2 - 2xy = 2-2dp
//            // d^2 = 2(1-dp)
//            // d = sqrt(2.0) * sqrt(1-dp);
//            const double SQUARE_ROOT_OF_2 = sqrt(2.0);
//            return SQUARE_ROOT_OF_2 * sqrt(1 - score);
//        }
//    };
//};

double dp2dist(double dp) {
    //
    // d^2 = x^2 + y^2 - 2xy = 2-2dp
    // d^2 = 2(1-dp)
    // d = sqrt(2.0) * sqrt(1-dp);
    const double SQUARE_ROOT_OF_2 = sqrt(2.0);
    return SQUARE_ROOT_OF_2 * sqrt(1 - dp);
}

// get normalized dp
double topN_dp_int_2_float(int score) {
    int maxScore = 0;
    int N = 50;
    for (int i = 1; i <= N; i++) {
        maxScore += i * i;
    }
    return score * 1.0 / maxScore;

}

class CRawDotProduct {
    int m_raw_dp;
    const int m_maxScore;
public:
    struct normalizer {
        int maxScore;

        normalizer(int topN) {
            maxScore = 0;
            for (int i = 1; i <= topN; i++) {
                maxScore += i * i;
            }
        }
    };

    CRawDotProduct(int score, int topN = 50) : m_maxScore(normalizer(topN).maxScore) {
        m_raw_dp = score;
    }

    double normalize() {
        return m_raw_dp * 1.0 / m_maxScore;
    }

    double distance() {
        //
        // d^2 = x^2 + y^2 - 2xy = 2-2dp
        // d^2 = 2(1-dp)
        // d = sqrt(2.0) * sqrt(1-dp);
        const double SQUARE_ROOT_OF_2 = sqrt(2.0);
        return SQUARE_ROOT_OF_2 * sqrt(1 - normalize());

    }

};


std::unordered_set<ClusterNode *> Clusters::findNeighboringNodes(ClusterNode *node, CAnnotationDB &mdb, double maxDist, bool usePepWithPTM) {
    static SimpleTimer st("find neighbor nodes");

    static int called_times = 0;
    long node_id = node->id;
    SPsmAnnotation gtinfo;
    mdb.retrieveGtinfo(node_id, gtinfo);

    node->peptide = usePepWithPTM ? gtinfo.getModifiedPeptide(false) : gtinfo.peptideseq; // updated
    // wait, we should use seq with mod
    // not found gtinfo with 0.8
//    node->peptide=gtinfo.getModifiedPeptide(false);
    node->isSignificant = gtinfo.isSig();
    SPsmAnnotation::IdxDistPair idxDistPairs;
    bool userealdist = true;
    if (gtinfo.m_neighbors == "NULL") {
        // we could update the gtinfo to get neighbors.
        double elpased_init = st.secondsElapsed();
        vector<long> indexlist(m_csa->getSpecNum(), 0);
        iota(indexlist.begin(), indexlist.end(), 0);
        vector<int> scores;

        // to check the tolerance tomorrow.
        bool normalize_dp = true;
        m_csa->scorePartiallyWithVecForm(50, 1, 32, normalize_dp, indexlist, m_csa->getSpecBy(node->id), scores);
        int num_neighbors = 3072;
        double elpased_scoring = st.secondsElapsed();
        partial_sort(indexlist.begin(), indexlist.begin() + num_neighbors, indexlist.end(),
                     [&](const long &x, const long &y) { return scores[x] > scores[y]; });


        double elpsed_sorting = st.secondsElapsed();
        for (int i = 0; i < scores.size() and i < num_neighbors; i++) {
            idxDistPairs.m_data.push_back(
                    SPsmAnnotation::IdxDistPair::idxDist{CRawDotProduct(scores[indexlist[i]]).distance(),
                                                         indexlist[i]});
        }
        double elapsed_conversion = st.secondsElapsed();
        cout << "node:\t" << node->id << "\tcalled:\t" << called_times++ << "\ttime:\t" << st.secondsElapsed()
             << "\tper call:\t" << st.secondsElapsed() / called_times << "\tscore\t" << elpased_scoring - elpased_init
             << "\tsorting\t" << elpsed_sorting - elpased_scoring << "\tconversion\t"
             << elapsed_conversion - elpsed_sorting << endl;


    } else {
        idxDistPairs = gtinfo.neighborStrToIdxDistPair();
//        auto idxDistPairs = gtinfo.neighborStrToIdxDistPair();
    }
    // remove if visited
    auto it = std::remove_if(idxDistPairs.m_data.begin(), idxDistPairs.m_data.end(),
                             [&](SPsmAnnotation::IdxDistPair::idxDist &z) -> bool {
                                 return z.dist > maxDist;
                             }
    );
    idxDistPairs.m_data.erase(it, idxDistPairs.m_data.end());
    unordered_set<ClusterNode *> ret;
    // this is where the error happens
    for (auto &z: idxDistPairs.m_data) {
        if (z.idx != node->id) ret.insert(m_nodeId2NodePtr[z.idx]);
    }
    // remove if out of distance.
    return ret;
}

int
Clusters::BFS(std::unordered_set<ClusterNode *> &neighborNodes, Cluster &cluster, CAnnotationDB &mdb, double maxDist,
              int maxSearchDepth, bool usePepWithPTM, bool alwaysMerge) {
    int addCounts = 0;
    while (not neighborNodes.empty()) {
        // list of neighbor not empty
        ClusterNode *nodePtr = *neighborNodes.begin();
        neighborNodes.erase(neighborNodes.begin());
        if (nodePtr->nodedepth > maxSearchDepth) {
//            cout << "current node depth " << nodePtr->nodedepth << "  is greater than threshold " << maxSearchDepth << endl;
            break;
        }
// breadth first search!
// looking for neighbors neighbor within maxSearchDepth
//  alwaysMergeMode: absorb any cluster encountered.
//  otherwise: only merge new node that are not visited.

        cluster.addNode(nodePtr);
        addCounts++;
        unordered_set<ClusterNode *> neighborNodes_2 = findNeighboringNodes(nodePtr, mdb, maxDist, usePepWithPTM);
        unordered_set<Cluster *> neighborClusters = findNeighboringClusters(neighborNodes_2);
        // if always merge, we should do something...
        if (alwaysMerge and neighborClusters.size() > 0) {
//            mergeClusters(neighborClusters)
            for (auto &eachcluster: neighborClusters) {
                if (eachcluster->cluster_id == cluster.cluster_id) {
                    // same cluster
                    continue;
                }
//                cout << "merge: " << eachcluster->cluster_id  << " into " << cluster.cluster_id << endl;
                eachcluster->mergeInto(&cluster);
                m_emptyClusterIds.push_back(eachcluster->cluster_id);
            }
        } else if (neighborClusters.size() > 0) {
            cluster.addNeighborClusters(neighborClusters);
        }
        removeVisited(neighborNodes_2);
        if (neighborNodes_2.size() > 0 and nodePtr->nodedepth < maxSearchDepth) {
            for (auto &eachitem: neighborNodes_2) {
                eachitem->nodedepth = nodePtr->nodedepth + 1;
            }
            neighborNodes.insert(neighborNodes_2.begin(), neighborNodes_2.end());
        }
    }
    return addCounts;
}

// do not rely on address of vector!!!
// this step is very slow....
// clusters...
Cluster *Clusters::createCluster() {
    if (not m_emptyClusterIds.empty()) {
        int id = *m_emptyClusterIds.begin();
        try {
            m_emptyClusterIds.erase(m_emptyClusterIds.begin());
            if (id == 11478)
                cout << "reusing node: " << id << " from " << getNumCluster() << " ptr=" << m_clusterPtrs.at(id)
                     << endl;
        } catch (...) {
            cout << "Error: unknow in create Cluster" << endl;
        }
        return m_clusterPtrs[id];
    } else {
        Cluster *p = new Cluster();
        if (p == nullptr) cout << "Error ptr" << p << endl;
        m_clusterPtrs.push_back(p);
//        int size = getNumCluster();
        p->cluster_id = getNumCluster() - 1;
        return p;
    }


}

void Clusters::print() {
    cout << "number of clusters : " << getNonEmptyCluster() << endl;
    for (int i = 0; i < 3 and i < getNumCluster(); i++) {
        if (m_clusterPtrs[i]->size() > 0) {
            // skip empty clusters
            cout << "cluster " << i << ": \n";
            m_clusterPtrs[i]->print();

        }
    }
    cout << "..." << endl;

}

void Clusters::saveAs(string filename) {
    ofstream fout(filename.c_str(), ios::out);
    // output summary also .
    vector<string> linesOfSummary;
    split_string(m_summary, linesOfSummary, '\n');
    for (auto &eachline : linesOfSummary) {
        fout << "S\t" << eachline << endl;
    }
    for (auto &item : m_clusterPtrs) {
        if (item->size() > 0)  // plot only when its not empty
        {
            fout << *item << endl;
        }
    }
    fout.close();
}

// todo: allow merging cause problems...
// to be continued...
// fixed problems: minpeak num 6
Clusters::~Clusters() {
    cout << "start releasing nodes" << endl;
//    getSummary();
    for (auto &x: m_clusterPtrs) {
        // release the ram
        if (x != nullptr) {
            delete x;
            x = nullptr;
        }
    }
}

std::unordered_set<Cluster *> Clusters::findNeighboringClusters(unordered_set<ClusterNode *> &neighborNodes) {
    std::unordered_set<Cluster *> neighborClusters;
    for (auto &node: neighborNodes) {
        if (node->isVisited) {
            neighborClusters.insert(m_clusterPtrs.at(node->clusterId));
        }
    }
    return neighborClusters;
}

void Clusters::removeVisited(unordered_set<ClusterNode *> &neighborNodes) {
    unordered_set<ClusterNode *> tmp;
// get set of cluster ids
    for (auto it = neighborNodes.begin(); it != neighborNodes.end(); it++) {
        if (not(*it)->isVisited) {
            tmp.insert(*it);
        }
    }
    neighborNodes.swap(tmp);
}

string Clusters::getPeptideEntropy() {
    // create pep int map
// remove visited nodes.
    // we have m_pepcounts for the total results.
    map<string, int> pep2intmap;
    map<string, int> pep2countMap;
    for (auto &x: m_clusterPtrs) {
        if (x->size() == 0) continue; // skip empty clusters
        if (x->m_pepCounts.size() == 1 and pep2intmap.find(x->m_pepCounts.begin()->first) == pep2intmap.end()) {
            // single peptide. new peptide
            string peptide = x->m_pepCounts.begin()->first;
            int counts = pep2intmap.size();
            pep2intmap[peptide] = counts;

        }
        if (x->m_pepCounts.size() == 1) {
            if (pep2countMap.find(x->m_pepCounts.begin()->first) == pep2countMap.end()) {
                pep2countMap[x->m_pepCounts.begin()->first] = 1;
            } else {
                pep2countMap[x->m_pepCounts.begin()->first]++;
            }
        }
    }

    // peptide A: cluster i, cluster j, cluster k
    //
    CountFrequency entropyOfPepID, entropyOfPepCounts;
    for (auto &x: m_clusterPtrs) {
        if (x->size() == 0) continue; // skip empty clusters
        if (x->m_pepCounts.size() == 1) {
            string peptide = x->m_pepCounts.begin()->first;
            int pep_int_id = pep2intmap[peptide];
            entropyOfPepID.add_data(pep_int_id);
        }
    }

    for (auto &x: pep2countMap) {
        entropyOfPepCounts.add_data(x.second);
    }
    // frequency of frequency
//    CountFrequency freqOfPepFreq;
//    for(auto &x: )
    // print the summary of the entropyOfPepID


    ostringstream oss;
    oss << "Entropy of peptide count distribution (the samller the better, low frequency dominant): "
        << entropyOfPepCounts.getEntropy(true)
        << " with " << pep2countMap.size() << " peptide counts" << endl;

//    for(auto &x: )

    oss << "Entropy of peptide distribution: " << entropyOfPepID.getEntropy(true)
        << " with " << entropyOfPepID.totalInsertion() << " clusters" << endl;
    oss << "peptide\tfrequency\tfraction" << endl;
    int total_insertions = entropyOfPepID.totalInsertion();
    for (auto &item: pep2intmap) {
        string peptide = item.first;
        int freq = entropyOfPepID.getFreq(item.second);
        oss << peptide << "\t" << freq << "\t" << freq * 1.0 / total_insertions << endl;
    }
    return oss.str();

}

void Clusters::getSummary() {
    auto &maxNodeCluster = *max_element(m_clusterPtrs.begin(), m_clusterPtrs.end(),
                                        [](const Cluster *x, const Cluster *y) {
                                            return x->size() < y->size();
                                        });
    int num_mixture_clusters = count_if(m_clusterPtrs.begin(), m_clusterPtrs.end(),
                                        [](const Cluster *x) {
                                            return not x->m_neighborClusters.empty();
                                        }
    );

    auto &maxMixCluster = *max_element(m_clusterPtrs.begin(), m_clusterPtrs.end(),
                                       [](const Cluster *x, const Cluster *y) {
                                           return x->m_neighborClusters.size() < y->m_neighborClusters.size();
                                       });
    ostringstream oss;
    oss << "Summary begin" << endl;
    oss << "Number of clusters: " << getNonEmptyCluster()
        << "\taverage node number: " << m_nodes.size() * 1.0 / getNonEmptyCluster() << endl;

    CountFrequency cf;

    for (auto &x: m_clusterPtrs) {
        if (x->size() > 0) // only for non empty clusters
            cf.add_data(x->size());
    }
    oss << "== Entropy of clusters ==: " << cf.getEntropy(true) << endl;
    CountFrequency cf_entropy_type;
    map<int, CountFrequency> type2SizeFrequency = {{-1, CountFrequency()},
                                                   {0,  CountFrequency()},
                                                   {1,  CountFrequency()}};
    map<int, unordered_set<string>> type2PepSet = {{-1, unordered_set<string>()},
                                                   {0,  unordered_set<string>()},
                                                   {1,  unordered_set<string>()}};
    for (auto &x: m_clusterPtrs) {
        if (x->size() == 0) continue; // skip empty clusters

        if (x->m_entropy < 0) {
            cf_entropy_type.add_data(-1);
            type2SizeFrequency[-1].add_data(x->size());
        } else if (x->m_entropy > 0) {
            cf_entropy_type.add_data(1);
            type2SizeFrequency[1].add_data(x->size());
            for (auto &eachPepCount: x->m_pepCounts) {

                type2PepSet[1].insert(eachPepCount.first);
            }
        } else {
            cf_entropy_type.add_data(0);
            type2SizeFrequency[0].add_data(x->size());
            type2PepSet[0].insert(x->m_pepCounts.begin()->first);
        }
    }
    oss << "Entropy of each cluster summarized (type code: -1: no pep; 0: 1 pep; 1: more than 1 pep)" << endl;
    map<int, string> typeCode2Name = {{0,  "pure(one_pep_ID)"},
                                      {-1, "????(no_pep_ID)"},
                                      {1,  "mixed(multi_peps)"}};
    oss << "typename\ttype\tcounts\tfrequency\t#nodes\tnodeFreq\tPepCounts" << endl;
    for (int i = -1; i <= 1; i++) {
        int counts = cf_entropy_type.getFreq(i);
        double frequency = counts * 1.0 / cf_entropy_type.totalInsertion();
        double nodefreq = type2SizeFrequency[i].getSampleSum() * 1.0 / cf.getSampleSum();//m_nodes.size();
        oss << typeCode2Name[i] << "\t" << i << "\t" << cf_entropy_type.getFreq(i) << "\t" << frequency
            << "\t" << type2SizeFrequency[i].getSampleSum() << "\t" << nodefreq << "\t" << type2PepSet[i].size()
            << endl;
    }
    oss << endl;

    vector<string> fourTypesOfClusters = {"singletons", "doubletons", "tripletons", "quadrupletons"};
    int numNodesRemaining = m_nodes.size();
    int numClusterRemaining = getNonEmptyCluster();
    for (int i = 0; i < fourTypesOfClusters.size(); i++) {
        int nodeSize = i + 1;
        int numClusterWithNodeSize = cf.getFreq(nodeSize);
        numNodesRemaining -= nodeSize * numClusterWithNodeSize;
        numClusterRemaining -= numClusterWithNodeSize;
        oss << "Number of " << fourTypesOfClusters[i] << ": " << numClusterWithNodeSize << "\tFraction: "
            << numClusterWithNodeSize * 1.0 / getNonEmptyCluster()
            << "\taverage node number of remaining clusters (size>" << nodeSize << "): "
            << numNodesRemaining * 1.0 / numClusterRemaining << endl;
    }

    oss << "summary of clusters of different size" << endl;
    numNodesRemaining = m_nodes.size();
    numClusterRemaining = getNonEmptyCluster();
    oss << "Node\t" << cf.getSampleSum() << endl;
    oss << "Clusters\t" << getNonEmptyCluster() << endl;
    int accumulateNodes = 0;
    int accumulateClusters = 0;
    oss << endl << "Size\t#clusters\tFraction\tacc#Clusters\tacc#Nodes\tAvgSizeRemaining" << endl;
    for (int i = 0; i < maxNodeCluster->size(); i++) {
        int nodeSize = i + 1;
        int numClusterWithNodeSize = cf.getFreq(nodeSize);
        if (numClusterWithNodeSize > 0) {
            accumulateNodes += nodeSize * numClusterWithNodeSize;
            numNodesRemaining = cf.getSampleSum() - accumulateNodes;
            accumulateClusters += numClusterWithNodeSize;
            numClusterRemaining = getNonEmptyCluster() - accumulateClusters;
            oss << nodeSize << "\t" << numClusterWithNodeSize
                << "\t" << numClusterWithNodeSize * 1.0 / getNonEmptyCluster()
                << "\t" << accumulateClusters
                << "\t" << accumulateNodes
                << "\t" << numNodesRemaining * 1.0 / numClusterRemaining << endl;

        }
    }
    oss << "end of summary for clusters" << endl;

    oss << "Max cluster: \nCluster ID "
        << maxNodeCluster->cluster_id << "\tsize: " << maxNodeCluster->size() << endl;
    oss << "Num of mixclusters (neighbor-size>0): "
        << num_mixture_clusters << endl;
    oss << "Cluster with maximum number of neighbors: Number of neighbor: "
        << maxMixCluster->m_neighborClusters.size() << "\tClusterId: "
        << maxMixCluster->cluster_id << "\tClusterSize: " << maxMixCluster->size() << endl;
//    maxMixCluster->print();
    oss << maxMixCluster->str() << endl;
    // Other important calculations
    // 1. How many dark matter illuninated in the pure clusters?
    // 2. How about the peptide being separated by multiple class? Entropy
    // 3.
    oss << getPeptideEntropy() << endl;
    oss << "Summay end" << endl;
    m_summary = oss.str();
    cout << m_summary << endl;
}

void Clusters::entropyOfEachCluster() {
    // by this step, the cluster has been created.
    for (auto &eachCluster: m_clusterPtrs) {
        if (eachCluster->size() == 0) continue;
//        if(eachCluster == nullptr) continue;
        map<string, int> pep2intMap;
        CountFrequency cf;

        for (auto &eachNode : eachCluster->m_nodeList) {
            if (not eachNode->isSignificant) continue; // not significant
            if (pep2intMap.find(eachNode->peptide) == pep2intMap.end()) {
                // not found in map
                int counts = pep2intMap.size();
                pep2intMap[eachNode->peptide] = counts + 1;
            }
//            pep2intMap[eachNode->peptide] = counts + 1;
            cf.add_data(pep2intMap[eachNode->peptide]);
        }
        if (cf.totalInsertion() > 0) {
            bool verbosity = eachCluster->cluster_id == 603362;
            eachCluster->m_entropy = cf.getEntropy(true, verbosity);
            for (auto &each: pep2intMap) {
                string pepseq = each.first;
                int pep_identifier = each.second;
                int counts = cf.getFreq(pep_identifier);
                eachCluster->m_pepCounts[pepseq] = counts;
            }
        } else {
            eachCluster->m_entropy = -1;
        }
        // refresh peptide entroty
        // refresh peptide occurence counts.


    }
    // refersh peptide occurrence counts

}

Cluster *Clusters::mergeClusters(unordered_set<Cluster *> &neighborClusters) {
    if (neighborClusters.size() == 1) {
        return *neighborClusters.begin();
    } else if (neighborClusters.empty()) {
        std::cout << "Error: neighbor clusters should not be empty" << endl;
        exit(0);
        // empty!!!
        return nullptr;
    } else {
        int min_cluster_id = (*neighborClusters.begin())->cluster_id;
        for (auto &x: neighborClusters) {
            if (min_cluster_id > x->cluster_id) {
// if there is one neighbor, merge to it
// if there are more than one,  find the one with min_cluster_id, every others merge to it!
// This step updates the emptyClusterid, which can be reused later.
                min_cluster_id = x->cluster_id;
            }
        }
        // every one merge to this one;
        for (auto &x: neighborClusters) {
            if (x->cluster_id != min_cluster_id) {
                m_clusterPtrs.at(x->cluster_id)->mergeInto(m_clusterPtrs.at(min_cluster_id));
                m_emptyClusterIds.push_back(x->cluster_id);

            }
        }
        return m_clusterPtrs.at(min_cluster_id);
    }

}


void ClusterDBScan(vector<CNode> &DB, double min_dp, double min_neighbor) {
    // cluster of everything.
//
//    std::vector<CCluster> clusters;
//    for(auto P: DB){
//        // DBSCAN:
//        // for each node in database, it the point has been visited, skip
//        if(P.isvisited){continue;}
//        // for new node, search for its neighbors, within a window of min_dp
//        DB.rangeSearch(P, min_dp);
//
//        // as it is new, node, we create a new cluster for it.
//        int CurrentClusterID = clusters.size();
//        CCluster aCluster(CurrentClusterID);
//        clusters.push_back(aCluster);
//        // if the size of the anns is less than threshold, it is noise node;
//        if(anns.size() < min_neighbor){
//            P.add_label(CNode::MINORITY);
//            P.add_label(CNode::VISITED);
//            P.setClusterID(CurrentClusterID);
//        }
//
//        set<CNode> annset (anns.begin(), anns.end());
//        while(not annset.empty())
//        {
//            CNode Q=annset.front();
//            // from ANN set get a member. check if it has been visited!
//            if(Q.visited){
//                // Yes, Q has been visited:
//                Q.addAddClusterID(CurrentClusterID); // one node with multiple ID is possible.
//            }
//
//        }
//        if(ann.size()>min_neighbor){
//            aCluster.push_back(CNode(CNode::CMemberType::CORE, spec));
//        } else{
//
//        }
//
//    }
}