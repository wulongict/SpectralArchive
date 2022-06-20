//
// Created by wulong on 9/2/20.
//

#ifndef MYTOOL_CLUSTERING_H
#define MYTOOL_CLUSTERING_H

// API
// member: visitedNode
//bool isNodeVisited(nodeIdentifier node_id){
// return visitiedNode.isExist(node_id);
// }

//for i in range(len(archiveSpec)):
//    query = archiveSpec[i]
//
//    ListOfNodeIdentifiers = getTopANNs(query);
//
//    refinedListNodes = refineList(ListOfNodeIdentifiers,threshold0.8)
//
//    refineListOfNodes = visitedNodes.refineVisitedOnes(refinedList)

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <memory>
#include <map>


class ClusterNode {

public:
    bool isVisited;
    long id;
    int clusterId;
    int nodedepth;
    bool isSignificant;
    std::string peptide;

    ClusterNode();

    ClusterNode(long id);

    std::string str();
};

class Cluster {
// a collection of node pointers
public:
    std::unordered_set<Cluster *> m_neighborClusters;
    int cluster_id;
    std::vector<ClusterNode *> m_nodeList;
    double m_entropy;
    std::map<std::string, int> m_pepCounts;

    void addNode(ClusterNode *nodePtr);

    void addNeighborClusters(std::unordered_set<Cluster *> neighborClusters);

    void print(int topN = 3);

    std::string str(int topN = 3);

    void mergeInto(Cluster *target);

    friend std::ostream &operator<<(std::ostream &os, const Cluster &node);

    int size() const;

    int getDarkNodesNum() const;
};

class CAnnotationDB;

class CMzFileReader;

class Clusters {
    std::string m_summary;
    // Number of spectra to be clustered
    long m_nodeNum;
    std::vector<Cluster *> m_clusterPtrs;
    std::vector<ClusterNode> m_nodes;
    std::map<long, ClusterNode *> m_nodeId2NodePtr;
    std::shared_ptr<CMzFileReader> m_csa;
    std::vector<int> m_emptyClusterIds;

    void getSummary();

    void entropyOfEachCluster(); // Entropy of size; might be other sourced ...
    std::string getPeptideEntropy();

public:

    Clusters(long nodeNum, std::shared_ptr<CMzFileReader> csa);

    long getNodeNum() {
        return m_nodeNum;
    }

    void print();

    void makeClusters(CAnnotationDB &mdb, double maxDist, long start, long end, bool mergeNode, int maxSearchDepth,
                      bool usePepWithPTM, bool alwaysMerge);

    bool isNodeVisited(int idx) {
        return m_nodes.at(idx).isVisited;
    }

    void setNodeDepth(int idx, int depth) {
        if (idx >= 0 and idx < getNodeNum())
            m_nodes[idx].nodedepth = depth;
    }

    Cluster *mergeClusters(std::unordered_set<Cluster *> &neighborClusters);

    int BFS(std::unordered_set<ClusterNode *> &neighborNodes, Cluster &cluster, CAnnotationDB &mdb, double maxDist,
            int maxSearchDepth, bool usePepWithPTM, bool alwaysMerge);

    std::unordered_set<ClusterNode *>
    findNeighboringNodes(ClusterNode *node, CAnnotationDB &mdb, double maxDist, bool usePepWithPTM);

    std::unordered_set<Cluster *> findNeighboringClusters(std::unordered_set<ClusterNode *> &neighborNodes);

    Cluster *createCluster();

    int getNumCluster() {
        return m_clusterPtrs.size();
    }

    int getNonEmptyCluster() {
        int counts = 0;
        for (auto &x: m_clusterPtrs) {
            if (x->size() > 0) counts++;
        }
        return counts;
    }

    void saveAs(std::string filename);

    ~Clusters();

    void removeVisited(std::unordered_set<ClusterNode *> &neighborNodes);
};


// Another tempt of the DBScan algorithm
class CNode {

    enum CMemberType {
        CORE, BORDER, NOISE, MINORITY
    };
    CMemberType m_type;
    std::vector<long> clusterId;
    std::vector<long> neighbors;
    long queryid;
public:
    bool isvisited;

    CNode(CMemberType mt, long queryid_) : m_type(mt), queryid(queryid_) {}
};

void ClusterDBScan(std::vector<CNode> &DB, double min_dp, double min_neighbor);


#endif //MYTOOL_CLUSTERING_H
