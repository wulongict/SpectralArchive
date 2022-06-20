//
// Created by wulong on 4/16/19.
//

#include "CKMeans.h"
#include <algorithm>
#include <memory>
#include <numeric>

#include "../../../librarymsms/Util.h"
#include "../../../librarymsms/CThreadsPool.h"

using namespace Eigen;

CKMeans::CKMeans(vector<vector<double>> &data, vector<vector<double>> &centroids, int k, bool verbose) : m_k(k) {
    _data = new CData(data);
    m_num = _data->size();
    if (m_num == 0) {
        cout << "No data provided!" << endl;
        exit(0);
    }
    m_membership.assign(m_num, 0);
    _centroids = new CData(centroids);
    m_release = true;
    m_verbose = verbose;
}

CKMeans::CKMeans(ICData *data, ICData *centroids, bool verbose) {
    _data = data;
    m_num = _data->size();
    _centroids = centroids;
    m_k = _centroids->size();

    m_membership.assign(m_num, 0);
    m_release = false;
    m_verbose = verbose;
}

CKMeans::CKMeans(shared_ptr<ICData> data, shared_ptr<ICData> centroids, bool verbose) : CKMeans(data.get(),
                                                                                                centroids.get(),
                                                                                                verbose) {

}

void CKMeans::run() {
    double minRate = 0.001;
    int max_iteration = 100;
    int k = 0;
    double totaldist = numeric_limits<double>::max();
    double newtotal_dist;
    while (k < max_iteration) {
        k++;
        if (m_verbose)cout << "Iteration: " << k << endl;

        update_membership(newtotal_dist);
        double rate = 0;
        if (totaldist > 0) {
            rate = (totaldist - newtotal_dist) / totaldist;
        }
        if (m_verbose)
            cout << "total distance: " << totaldist << " ---> Decreased " << rate * 100 << "% to " << newtotal_dist
                 << endl;
        totaldist = newtotal_dist;
        if (rate < minRate or rate < 0) {
            break;
        }

        update_centroids();

    }
    cout << "K-Means finished with " << k << " iterations" << endl;

}

void CKMeans::display() {
    cout << "Sample" << endl;

    for (int i = 0; i < m_num; i++) {
        cout << "sample " << i << " in cluster:   " << m_membership[i] << ":";
        _data->getVec(i).display();
    }

    cout << "Centroid:" << endl;
    _centroids->display();

}


void CKMeans::update_centroid(int j) {
    vector<int> idxlist;
    vector<int> allidx(m_num, 0);
    std::iota(allidx.begin(), allidx.end(), 0);
    std::copy_if(allidx.begin(), allidx.end(), std::back_inserter(idxlist),
                 [&](const int a) { return m_membership[a] == j; });

    _data->average(idxlist, _centroids->getVec(j));

}

void CKMeans::update_centroids() {

    for (int i = 0; i < m_k; i++) {
        update_centroid(i);
    }
}


class CUpdateMembershipTask : public ICThreadTask {
    int m_id;
    CKMeans *m_kmeansobj;
    double &m_total_dist;
    mutex &m_totalDistMutex;
public:
    CUpdateMembershipTask(int i, CKMeans *km, double &total_dist, mutex &totalDistMutex) : m_total_dist(total_dist),
                                                                                           m_totalDistMutex(
                                                                                                   totalDistMutex) {
        m_id = i;
        m_kmeansobj = km;
    }

    void run() override {
        double val = 0;
        m_kmeansobj->update_membership(m_id, val);
        updateTotalDist(val);

    }

    void updateTotalDist(double val) {
        std::lock_guard<std::mutex> guard(m_totalDistMutex);
        m_total_dist += val;
    }
};


bool CKMeans::update_membership(double &totaldist) {
    mutex totalDistMutex;
    double total_dist = 0;
    vector<CUpdateMembershipTask> vtask;
    for (int i = 0; i < m_num; i++) {
        vtask.emplace_back(CUpdateMembershipTask(i, this, total_dist, totalDistMutex));
    }
    vector<ICThreadTask *> vtaskptrs;
    for (int i = 0; i < m_num; i++) {
        vtaskptrs.push_back(&vtask[i]);
    }

    CTaskPool ctp(vtaskptrs, false, m_verbose);
    int threadnum = getProperThreads() / 3 * 2 + 1;
    ctp.start(threadnum, "UpdateMembership");
    totaldist = total_dist;
    return true;


}

bool CKMeans::update_membership(int i, double &dist) {
    // update membership of data i
    bool updated = false;
    int newid = _centroids->getClusterID((*_data)[i]);
    dist = _centroids->getVec(newid).distance(_data->getVec(i));

    if (m_membership[i] != newid) {
        m_membership[i] = newid;
        updated = true;
    }
    return updated;
}

CKMeans::~CKMeans() {
    if (m_release) {
        delete _data;
        delete _centroids;

    }
}

void CKMeans::getResidual(vector<vector<double>> &res) {
    int dim = _data->getVec(0).size();
    res.assign(m_num, vector<double>(dim, 0));
    for (int i = 0; i < m_num; i++) {
        int k = m_membership[i];
        for (int j = 0; j < dim; j++) {
            res[i][j] = _data->getVec(i).get(j) - _centroids->getVec(k).get(j);
        }
    }

}

void CKMeans::getNearestId(vector<int> &ids) {
    if (_centroids == nullptr) {
        cout << "Error: the model is not trained yet!!" << endl;
    }
    cout << "start to get ids: " << _centroids->size() << endl;

    ids.assign(_centroids->size(), 0);
    iota(ids.begin(), ids.end(), 0);
    for (int i = 0; i < ids.size(); i++) {
        double dist = numeric_limits<double>::max();
        //int id = -1;
        for (int j = 0; j < _data->size(); j++) {
            if (m_membership[j] != i) continue;
            double tmp = _centroids->getVec(i).distance(_data->getVec(j));
            if (tmp < dist) {
                dist = tmp;
                ids[i] = j;
            }
        }

    }
}


void kmeans(ICData *data, ICData *centroids) {
    int maxIter = 10;
    int n = 0;
    vector<int> membership(data->size(), 0);
    bool updated = false;

    while (n < maxIter) {
        n++;
        // update membership
        for (int i = 0; i < data->size(); i++) {
            int id = centroids->getClusterID(data->getVec(i));
            if (id != membership[i]) {
                membership[i] = id;
                updated = true;
            }
        }

        if (not updated) {
            break;
        }

        // update centroid
        for (int i = 0; i < centroids->size(); i++) {
            vector<int> idxlist;
            for (int j = 0; j < data->size(); j++) {
                if (membership[j] == i) {
                    idxlist.push_back(j);
                }
            }
            data->average(idxlist, centroids->getVec(i));
        }
    }
    // centroids->display();
}

void ICKMeansVec::display() {
    cout << "[";
    for (int i = 0; i < size(); i++) {
        if (i > 0) cout << "\t";
        cout << get(i);

    }
    cout << "]" << endl;
}

double ICKMeansVec::distance(ICKMeansVec &other) {
    double sum = 0;
    for (int i = 0; i < other.size(); i++) {
        double diff = other[i] - get(i);
        sum += diff * diff;
    }
    return sqrt(sum);
}

ICKMeansVec &ICKMeansVec::operator+=(ICKMeansVec &other) {
    for (int i = 0; i < size(); i++) {
        set(i, get(i) + other[i]);
    }
    return *this;
}

void ICKMeansVec::initAsZero() {
    for (int i = 0; i < size(); i++) {
        set(i, 0);
    }
}

ICKMeansVec &ICKMeansVec::operator/=(double x) {

    if (x > 0.00000000001) {
        for (int i = 0; i < size(); i++) {
            double v = get(i);
            set(i, v / x);
        }

    }
    return *this;
}

double ICKMeansVec::operator[](int i) {
    return get(i);
}

void ICData::display() {
    cout << "[" << endl;
    for (int i = 0; i < size(); i++) {
        getVec(i).display();
    }
    cout << "]" << endl;
}

void ICData::average(vector<int> &idxlist, ICKMeansVec &out) {
    out.initAsZero();
    for (int i = 0; i < idxlist.size(); i++) {
        int j = idxlist[i];
        out += getVec(j);
    }
    out /= idxlist.size();

}

int ICData::getClusterID(ICKMeansVec &sample) {
    double maxdist = numeric_limits<double>::max();
    int id = -1;
    for (int i = 0; i < size(); i++) {
        double d = getVec(i).distance(sample);
        if (d < maxdist) {
            id = i;
            maxdist = d;
        }
    }

    return id;
}

ICKMeansVec &ICData::operator[](int i) { return getVec(i); }

double CVec::get(int i) {
    return m_vec[i];
}

int CVec::size() {
    return m_vec.size();
}

void CVec::set(int i, double val) {
    m_vec[i] = val;
}

CVec::CVec(vector<double> &v) : m_vec(v) {

}

ICKMeansVec &CData::getVec(int i) {
    return m_Vecs[i];
}

int CData::size() {
    return m_Vecs.size();
}

CData::CData(vector<vector<double>> &data) : m_data(data) {
    for (int i = 0; i < data.size(); i++) {
        m_Vecs.emplace_back(CVec(m_data[i]));
    }
}

CProductQuantization::CProductQuantization() : CProductQuantization(0, 0, 0) {}

CProductQuantization::CProductQuantization(int PQN, int IVF, int DIM) {
    m_subSpaceNum = PQN;
    m_clusterNum = IVF;
    m_vectorSize = DIM;

    m_indexCodeBookObj.setDim(PQN + 1);

    if (PQN == 0 or DIM % PQN != 0) {
        cout << "The zero means DIM " << DIM << " can not be divided by PQN " << PQN << "!" << endl;
    }

    m_initCenPtr = nullptr;
    m_data = nullptr;
    istrained = false;
    m_nprobe = 1;

    m_preCompDistPtr = nullptr;

}

vector<int> CProductQuantization::getIdx(int i) {
    int subspaceSize = m_vectorSize / m_subSpaceNum;
    vector<int> a(subspaceSize, 0);
    iota(a.begin(), a.end(), i * subspaceSize);
    return a;
}


void CProductQuantization::train(ICData *data, CPQParam &cpqParam) {
    chooseProperInitialVec(data, cpqParam);

    m_initCenPtr = make_shared<CData>(m_initCenMatrix);
    m_data = data;

    CKMeans kmeansobj(m_data, m_initCenPtr.get(), true);
    kmeansobj.run();

    vector<vector<double>> residual;

    kmeansobj.getResidual(residual);
    m_residualCenMatrix.assign(residual.begin(), residual.begin() + m_clusterNum);
    cout << "start kmeans on residuals of the first round" << endl;

    for (int i = 0; i < m_subSpaceNum; i++) {
        cout << "subspace : " << i << endl;
        vector<int> idx_selected = getIdx(i);
        shared_ptr<ICData> nv(new CDataSelected(residual, idx_selected));
        shared_ptr<ICData> nc(new CDataSelected(m_residualCenMatrix, idx_selected));
        CKMeans residualKmeans(nv, nc, false);
        residualKmeans.run();
    }
    istrained = true;
    cout << "kmeans is done on residuals" << endl;

    initPreComputeDist();

}

CProductQuantization::~CProductQuantization() {
    cout << "Destruction of Product Quantization object " << endl;
}


void CProductQuantization::search(float *queries, int querynum, CPQResult &res) {
    initPreComputeDist();

    int dim = getDim();
    if (m_nprobe > m_clusterNum or m_nprobe < 1) m_nprobe = 1;

    vector<int> idx_list;
    idx_list.reserve(size());
    vector<double> dist_list;
    dist_list.reserve(size());

    createIdsInEachBuckets();

    for (int i = 0; i < querynum; i++) {
        vector<double> dist(m_clusterNum, 0);
        shared_ptr<ICKMeansVec> a(new CVecFloatPtr(queries + i * dim, dim));
        for (int k = 0; k < dist.size(); k++) {
            dist[k] = m_initCenPtr->getVec(k).distance(*a);
        }
        vector<int> topBucketId(m_clusterNum, 0);
        iota(topBucketId.begin(), topBucketId.end(), 0);
        sort(topBucketId.begin(), topBucketId.end(), [&](const int x, const int y) { return dist[x] < dist[y]; });
        topBucketId.resize(m_nprobe);

        idx_list.resize(0);
        dist_list.resize(0);
        createPreComputeDist(a);

        SimpleTimer st("calculate dist");
        // use multiple thread here

        for (int j = 0; j < topBucketId.size(); j++) {
            calculateDistWithCodeBook(topBucketId[j], idx_list, dist_list);
        }
        res.add(idx_list, dist_list);

    }
}

vector<double> CProductQuantization::reconstruct(vector<int> &code) {
    vector<double> result(m_vectorSize, 0);
    vector<double> tmp(m_vectorSize, 0);
    shared_ptr<ICKMeansVec> x(new CVec(result));
    *x += m_initCenPtr->getVec(code[0]);
    for (int i = 0; i < m_subSpaceNum; i++) {
        vector<int> idx_selected = getIdx(i);
        shared_ptr<ICKMeansVec> y(new CVecSelected(result, idx_selected));
        shared_ptr<ICKMeansVec> z(new CVecSelected(m_residualCenMatrix[code[i + 1]], idx_selected));
        *y += *z;
    }
    return result;
}

void CProductQuantization::display(bool all) {
    cout << "--------------Product Quantization begin---------" << endl;
    cout << "dim = " << m_vectorSize << "\t\t// vector size" << endl;
    cout << "K = " << m_clusterNum << "\t\t// k for K-Means algorithm; Num of clusters" << endl;
    cout << "subspace num = " << m_subSpaceNum << "\t\t// Number of subspaces" << endl;
    // cout << "precomputeDist size = " << m_precomputeDist.size() << "\t\t the size of precomputed distance" << endl;

    if (all) {
        cout << "Centroid: ";
        shared_ptr<CData> p(new CData(m_initCenMatrix));  // start using make_shared!
        p->display();
        shared_ptr<CData> q(new CData(m_residualCenMatrix));
        cout << "Centroid of Residual Subspaces: ";
        q->display();
        cout << "Index Code Book" << endl;
        for (int j = 0; j < size(); j++) {
            for (int i = 0; i < m_subSpaceNum + 1; i++) {
                if (i != 0) cout << "\t";
                cout << getFromIndexCodeBook(j, i);
            }
            cout << endl;
        }
    }
    cout << "-------------------Product Quantization end--------" << endl;
}

void CProductQuantization::load(string filename) {
    cout << "start loading " << endl;
    bool first_centroid_batch = true;
    int idx_in_batch = 0;
    ifstream fin(filename, ios::in);
    string line;
    while (getline(fin, line)) {
        if (line.size() < 2) continue;
        else if (line.substr(0, 2) == "//") {
            continue;
        } else {
            int found = line.find_first_of("=");
            if (found != string::npos) {
                // find key value pairs
                string key = line.substr(0, found);
                string val = line.substr(found + 1);
                if (key == "version") {
                    string version = val;
                } else if (key == "dim") {
                    m_vectorSize = atoi(val.c_str());
                } else if (key == "subspace") {
                    m_subSpaceNum = atoi(val.c_str());
                } else if (key == "cluster") {
                    m_clusterNum = atoi(val.c_str());
                }

            } else {
                // no equal sign
                if (first_centroid_batch) {
                    if (m_initCenMatrix.size() == 0) {
                        m_initCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));
                    }

                    istringstream iss(line);
                    for (int i = 0; i < m_vectorSize; i++) {
                        iss >> m_initCenMatrix[idx_in_batch][i];
                    }
                    idx_in_batch++;
                    if (idx_in_batch >= m_clusterNum) {
                        first_centroid_batch = false;
                        idx_in_batch = 0;
                    }

                } else {
                    if (m_residualCenMatrix.size() == 0) {
                        m_residualCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));

                    }

                    istringstream iss(line);
                    for (int i = 0; i < m_vectorSize; i++) {
                        iss >> m_residualCenMatrix[idx_in_batch][i];
                    }
                    idx_in_batch++;
                    if (idx_in_batch >= m_clusterNum) first_centroid_batch = false;
                }
            }

        }
    }
    fin.close();
    if (m_residualCenMatrix.size() == m_clusterNum) {
        istrained = true;
    } else {
        cout << "Error: index is not trained!!!" << endl;
    }
    m_initCenPtr = make_shared<CData>(m_initCenMatrix);


    initPreComputeDist();
    //m_preCompDistPtr = make_shared<CPreComputeDist>(m_initCenMatrix, m_residualCenMatrix, m_subSpaceNum);
    //m_precomputeDist.assign(m_clusterNum*m_clusterNum*m_subSpaceNum,0);
    m_indexCodeBookObj.setDim(m_subSpaceNum + 1);
    //m_indexCodeBookObj.setDim(m_vectorSize);

    string codefile = filename + ".code";
    fin.open(codefile, ios::in);
    while (getline(fin, line)) {
        istringstream iss(line);
        resizeCodeIndex(1);

        for (int i = 0; i < getCodeBookDim(); i++) {
            int val;
            iss >> val;
            setIndexCodeBook(size() - 1, i, val);
        }
    }

    display(false);
}

void CProductQuantization::save(string filename) {
    // save the model we have trained;
    // m_initCenPtr
    ofstream fout(filename, ios::out);
    fout << "// Program Build on " << __DATE__ << endl;
    fout << "version=1.0" << endl;
    fout << "dim=" << m_vectorSize << endl;
    fout << "subspace=" << m_subSpaceNum << endl;
    fout << "cluster=" << m_clusterNum << endl;
    fout << "// centoroid " << endl;
    for (int i = 0; i < m_initCenPtr->size(); i++) {
        fout << m_initCenPtr->getVec(i) << endl;
    }
    fout << "// centroid PQN" << endl;
    // m_centroidsPQN
    shared_ptr<ICData> centroidPQN(new CData(m_residualCenMatrix));
    for (int i = 0; i < centroidPQN->size(); i++) {
        fout << centroidPQN->getVec(i) << endl;
    }
    cout << "saved to file " << filename << endl;
    fout.close();

    string codefile = filename + ".code";
    fout.open(codefile, ios::out);
    for (int i = 0; i < size(); i++) {
        for (int j = 0; j < getCodeBookDim(); j++) {
            if (j != 0) fout << "\t";
            fout << getFromIndexCodeBook(i, j);
        }
        fout << endl;
    }
    fout.close();
}


void CProductQuantization::add(int numspec, float *ptr, int d) {
    long len = size();
    resizeCodeIndex(numspec);
    Map<Matrix<float, Dynamic, Dynamic, RowMajor>> newSpec(ptr, numspec, d);

    Matrix<float, Dynamic, Dynamic> centroidsM(m_clusterNum, d);
    for (int i = 0; i < m_initCenMatrix.size(); i++) {
        for (int j = 0; j < m_initCenMatrix[0].size(); j++) {
            centroidsM(i, j) = m_initCenMatrix[i][j];
        }
    }

    auto centroidsNorm = centroidsM.rowwise().squaredNorm();
    auto dataNorm = newSpec.rowwise().squaredNorm();

    // Now do the multiplication
    MatrixXf res = newSpec * centroidsM.transpose();

    res *= -2.0f;
    auto finalres = (res.rowwise() + centroidsNorm.transpose()).colwise() + dataNorm;

    vector<int> z(numspec, 0);
    for (int i = 0; i < finalres.rows(); i++) {
        finalres.row(i).minCoeff(&z[i]);
        setIndexCodeBook(len + i, 0, z[i]);
    }

    // subspace matrix
    Matrix<float, Dynamic, Dynamic> Res_centroidsM(m_clusterNum, d);
    for (int i = 0; i < m_residualCenMatrix.size(); i++) {
        for (int j = 0; j < m_residualCenMatrix[0].size(); j++) {
            Res_centroidsM(i, j) = m_residualCenMatrix[i][j];
        }
    }

    MatrixXf resdualdata(newSpec.rows(), newSpec.cols());
    for (int i = 0; i < resdualdata.rows(); i++) {
        resdualdata.row(i) = newSpec.row(i) - centroidsM.row(z[i]);
    }

    for (int i = 0; i < m_subSpaceNum; i++) {
        // get the subspace id
        // step 1 get subspace matrix
        int subspaceCols = d / m_subSpaceNum;
        int startCol = i * subspaceCols;
        MatrixXf resMXf = resdualdata.block(0, startCol, numspec, subspaceCols);
        MatrixXf resCtrXf = Res_centroidsM.block(0, startCol, m_clusterNum, subspaceCols);

        auto centroidsNorm = resCtrXf.rowwise().squaredNorm();
        auto dataNorm = resMXf.rowwise().squaredNorm();

        MatrixXf resdp = resMXf * resCtrXf.transpose();

        resdp *= -2.0f;
        auto finalres = (resdp.rowwise() + centroidsNorm.transpose()).colwise() + dataNorm;

        vector<int> z(numspec, 0);
        for (int j = 0; j < finalres.rows(); j++) {
            finalres.row(j).minCoeff(&z[j]);
            setIndexCodeBook(len + j, i + 1, z[j]);
        }
    }
}


// single thread now
void CProductQuantization::calculateDistWithCodeBook(int idx, vector<int> &idx_list, vector<double> &dist_list) {
    if (m_idsInEachBuckets.size() <= idx) {
        cout << "size of ids in each bukcets " << m_idsInEachBuckets.size() << endl;
        throw "invalid idx";
    }

    for (int i = 0; i < m_idsInEachBuckets[idx].size(); i++) {
        int j = m_idsInEachBuckets[idx][i];
        double dist = 0;
        for (int k = 0; k < m_subSpaceNum; k++) {
            int subspace = getFromIndexCodeBook(j, k + 1);
            dist += getPreCompDist(idx, subspace, k);

        }
        idx_list.push_back(j);
        dist_list.push_back(dist);

    }


}

void CProductQuantization::createPreComputeDist(shared_ptr<ICKMeansVec> &a) {
    m_preCompDistPtr->create(a);
}

void CProductQuantization::initPreComputeDist() {
    if (m_precomputeDist.size() == 0) {
        // todo: to be deleted:
        m_precomputeDist.assign(m_clusterNum * m_clusterNum * m_subSpaceNum, 0);
    }
    if (istrained and m_preCompDistPtr == nullptr) {
        m_preCompDistPtr = make_shared<CPreComputeDist>(m_initCenMatrix, m_residualCenMatrix, m_subSpaceNum);
    }

}

float CProductQuantization::getPreCompDist(int iniCenId, int resCenId, int spaceId) {
    return m_preCompDistPtr->get(iniCenId, resCenId, spaceId);
//    int pos = iniCenId * m_clusterNum * m_subSpaceNum + resCenId * m_subSpaceNum + spaceId;
//    return m_precomputeDist.at(pos);
}

void CProductQuantization::createIdsInEachBuckets() {
    if (m_idsInEachBuckets.size() == 0) {
        m_idsInEachBuckets.assign(m_clusterNum, vector<long>());
        for (int i = 0; i < m_indexCodeBookObj.size(); i++) {
            int idx = m_indexCodeBookObj.get(i, 0);
            m_idsInEachBuckets[idx].push_back(i);
        }

    }
}

void CProductQuantization::chooseProperInitialVec(ICData *data, CPQParam &cpqParam) {

    if (cpqParam._option == CPQParam::INITIAL_ONES) {
        useFirstFewAsCentroid(data);
    } else if (cpqParam._option == CPQParam::RANDOM_SAMPLES) {
        useRandomVecAsCentroid(data);
    } else if (cpqParam._option == CPQParam::RANDOM) {
        useGaussianRandomVecAsCentroid(data);
    } else if (cpqParam._option == CPQParam::INITIAL_ONES) {
        useRandomProjection(data);
    } else {
        cout << "Invalid option: " << cpqParam._option << endl;
    }
}

void CProductQuantization::useRandomVecAsCentroid(ICData *data) {
    cout << "initializing centroids with randomly chosen vectors " << endl;
    if (m_initCenMatrix.size() == 0) {
        m_initCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));
    }
    int len = data->size();
    vector<int> idx(len, 0);
    iota(idx.begin(), idx.end(), 0);

    gcc5shuffle(idx.begin(), idx.end(), mt19937(42));
    for (int i = 0; i < m_clusterNum; i++) {
        for (int j = 0; j < m_vectorSize; j++) {
            m_initCenMatrix[i][j] = data->getVec(idx[i]).get(j);
        }
    }

}

void CProductQuantization::useGaussianRandomVecAsCentroid(ICData *data) {
    cout << "initializing centroids with random vector of uniform r.v. " << endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    if (m_initCenMatrix.size() == 0) {
        m_initCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));
    }

    for (int i = 0; i < m_clusterNum; i++) {
        for (int j = 0; j < m_vectorSize; j++) {
            m_initCenMatrix[i][j] = distribution(generator);
        }
        // normalization
        double norm = 0;
        for (int j = 0; j < m_vectorSize; j++) {
            norm += m_initCenMatrix[i][j] * m_initCenMatrix[i][j];
        }
        norm = sqrt(norm);

        // divide by norm
        for (int j = 0; j < m_vectorSize; j++) {
            m_initCenMatrix[i][j] /= norm;
        }
    }
}

void CProductQuantization::useRandomProjection(ICData *data) {
    cout << "initializing centroids with random projection" << endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int d = 2;
    if (m_clusterNum < d) d = m_clusterNum;
    vector<vector<double>> projection(d, vector<double>(m_vectorSize, 0));

    for (int i = 0; i < d; i++) {
        for (int j = 0; j < m_vectorSize; j++) {
            projection[i][j] = distribution(generator);
        }
    }
    cout << "we have projection vectors " << endl;
    // now do projection
    vector<vector<double>> newdata(data->size(), vector<double>(d, 0));
    for (int i = 0; i < data->size(); i++) {
        for (int j = 0; j < d; j++) {
            double val = 0;
            for (int k = 0; k < m_vectorSize; k++) {
                val += (*data)[i][k] * projection[j][k];
            }
            newdata[i][j] = val;
        }
    }
    cout << "we just created new data with projection" << endl;
    // after projection, we have new data
    vector<vector<double>> centroids(m_clusterNum, vector<double>(d, 0));
    for (int i = 0; i < m_clusterNum; i++) {
        centroids[i].assign(newdata[i].begin(), newdata[i].end());
    }
    cout << "we have created the centroid with new data" << endl;
    CData ndata(newdata);
    CData cen(centroids);
    CKMeans kmeansjob(&ndata, &cen, false);
    kmeansjob.run();
    cout << "We are going to do k means now " << endl;
    vector<int> ids;
    kmeansjob.getNearestId(ids);
    cout << "those ids: " << endl;
    for (auto each: ids) cout << " " << each;
    cout << endl;

    if (m_initCenMatrix.size() == 0) {
        m_initCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));
    }

    cout << "now start copy data " << endl;
    for (int i = 0; i < m_clusterNum; i++) {
        for (int j = 0; j < m_vectorSize; j++) {
            m_initCenMatrix[i][j] = (*data)[ids[i]][j];
        }
    }
    cout << "use the ids we have to init the matrix of centroid" << endl;
}

void CProductQuantization::useFirstFewAsCentroid(ICData *data) {
    cout << "initializing centroids with the first few vector " << endl;
    if (m_initCenMatrix.size() == 0) {
        m_initCenMatrix.assign(m_clusterNum, vector<double>(m_vectorSize, 0));
    }
    for (int i = 0; i < m_clusterNum; i++) {
        for (int j = 0; j < m_vectorSize; j++) {
            m_initCenMatrix[i][j] = data->getVec(i).get(j);
        }
    }
}


CDataSelected::CDataSelected(vector<vector<double>> &data,
                             vector<int> idx_selected) : m_data(data) {
    m_idx_selected = idx_selected;
    for (int i = 0; i < data.size(); i++) {
        m_Vecs.emplace_back(CVecSelected(data[i], m_idx_selected));
    }
}

ICKMeansVec &CDataSelected::getVec(int i) {
    return m_Vecs[i];
}

int CDataSelected::size() {
    return m_Vecs.size();
}

// Question: does vector of vector has a fixed ptr?
// todo: ------ change the initialization step
void CPQResult::add(vector<int> &idx, vector<double> &dist) {
    SimpleTimer st("Add index to result");
    m_res.emplace_back(vector<CPQResultEntry>());
    for (int i = 0; i < idx.size(); i++) {
        m_res.back().emplace_back(CPQResultEntry(idx[i], dist[i]));
    }
    vector<CPQResultEntry> &r = m_res.back();
    partial_sort(r.begin(), r.begin() + m_topN, r.end(),
                 [](const CPQResultEntry &a, const CPQResultEntry &b) { return a.m_dist < b.m_dist; });
    r.resize(m_topN);
}

void CPQResult::display() {
    for (auto each: m_res) {
        for (auto eachEntry: each) {
            cout << eachEntry << endl;
        }
    }

}

void CPQResult::add(vector<CPQResultEntry> &x) {
    m_res.push_back(x);
}

const vector<CPQResultEntry> &CPQResult::get(int i) const {
    return m_res[i];
}

int CPQResult::size() const {
    return m_res.size();
}

bool CPQResult::operator==(const CPQResult &other) const {
    if (size() != other.size()) return false;

    for (int i = 0; i < other.size(); i++) {
        const vector<CPQResultEntry> &a = other.get(i), &b = get(i);
        if (a.size() != b.size()) return false;

        for (int j = 0; j < a.size(); j++) {
            if (a[j].m_idx != b[j].m_idx or a[j].m_dist != b[j].m_dist) return false;
        }
    }
    return true;
}

bool CPQResult::equal(const CPQResult &other) {
    return *this == other;
}

void CPQResult::flattern(vector<float> &dist, vector<long> &ind, int d) {
    for (int i = 0; i < m_res.size(); i++) {
        int k = i * d;
        const vector<CPQResultEntry> &cpqr_entry = m_res[i];
        for (int j = 0; j < cpqr_entry.size(); j++) {
            ind[k] = cpqr_entry[j].m_idx;
            dist[k] = cpqr_entry[j].m_dist;
            k++;
        }
    }
}

CVecSelected::CVecSelected(vector<double> &v, vector<int> &idx_selected) : m_vec(v), m_idx_selected(idx_selected) {

}

double CVecSelected::get(int i) {
    return m_vec[m_idx_selected[i]];
}

void CVecSelected::set(int i, double val) {
    m_vec[m_idx_selected[i]] = val;
}

int CVecSelected::size() {
    return m_idx_selected.size();
}

CPQResultEntry::CPQResultEntry(int idx, double dist) {
    m_idx = idx;
    m_dist = dist;
}

CPQResultEntry::CPQResultEntry() {
    m_idx = -1;
    m_dist = numeric_limits<double>::max();
}

ostream &operator<<(ostream &os, const CPQResultEntry &cpqrEntry) {
    os << cpqrEntry.m_idx << "\t" << cpqrEntry.m_dist;
    return os;
}

void CPQResultEntry::set(int idx, double dist) {
    m_idx = idx;
    m_dist = dist;
}

void CPreComputeDist::create(shared_ptr<ICKMeansVec> &a) {
    for (int idx = 0; idx < m_clusterNum; idx++) {
        for (int k = 0; k < m_dim; k++) {
            m_resM(idx, k) = a->get(k) - m_initCenMatrix[idx][k];//m_initCenPtr->getVec(idx).get(k);
            //resCenM(idx, k) = m_residualCenMatrix[idx][k];
        }
    }

    // Now consider the distance between each subspace;
    int subsize = m_dim / m_subspaceNum;
    for (int i = 0; i < m_subspaceNum; i++) {
        MatrixXf subM = m_resM.block(0, i * subsize, m_clusterNum, subsize); // the sub matrix of res
        MatrixXf subCenM = m_resCenM.block(0, i * subsize, m_clusterNum, subsize);

        // Norm
        auto subMnorm = subM.rowwise().squaredNorm();
        auto subCenMnorm = subCenM.rowwise().squaredNorm();

        // Now distance
        MatrixXf dp = subM * subCenM.transpose();
        dp *= -2.0f;

        m_preDist[i] = (dp.colwise() + subMnorm).rowwise() + subCenMnorm.transpose();
//        for(int idx = 0; idx< m_clusterNum; idx ++)
//        {
//            for(int j = 0; j < m_clusterNum; j ++)
//            {
//                setPreCompDist(idx, j, i, finalres(idx,j));
//                // m_precomputeDist[idx*m_clusterNum*m_subSpaceNum+j*m_subSpaceNum+i]=finalres(idx,j);
//
//            }
//        }

    }
}

CPreComputeDist::CPreComputeDist(vector<vector<double>> &initCen, vector<vector<double>> &resCen, int subspaceNum)
        : m_initCenMatrix(initCen), m_resCenMatrix(resCen) {
    m_clusterNum = m_initCenMatrix.size();
    if (m_clusterNum == 0) {
        cout << "Error: invalid cluster num: " << m_clusterNum << endl;
    }
    m_dim = m_initCenMatrix[0].size();
    m_subspaceNum = subspaceNum;
    m_preDist.resize(m_subspaceNum);


    SimpleTimer st("precompute distance");

    m_resM.resize(m_clusterNum, m_dim), m_resCenM.resize(m_clusterNum, m_dim);
    for (int idx = 0; idx < m_clusterNum; idx++) {
        for (int k = 0; k < m_dim; k++) {
            //resM(idx, k) = a->get(k) - m_initCenPtr->getVec(idx).get(k);
            m_resCenM(idx, k) = m_resCenMatrix[idx][k];
        }
    }

}

void CPQParam::display() {
    cout << "init_option = " << _option << endl;
}
