//
// Created by wulong on 10/3/19.
//

#include <random>
#include "pvalueCalc.h"
#include "../../../librarymsms/CThreadsPool.h"
#include "CTailEstimation.h"
#include "Visual.h"
#include "../../../librarymsms/CFragmentIndex.h"
#include "../../../librarymsms/ConcreteQueries.h"
#include "../../../librarymsms/CArchiveSearchReply.h"
#include "../../../librarymsms/CMzFileReader.h"
#include "../../../librarymsms/Util.h"
#include "../../../librarymsms/ProteomicsDataTypes.h"


// not used now
class DPTask : public ICThreadTask {
    uint16_t *m_ptr;
    long m_query_idx;
    CAnnSpectra *m_qr;
    vector<long> &m_indexlist;
    CMzFileReader *m_p;
    long m_ArchiveSpecNum;
public:
    DPTask(uint16_t *ptr, long query_idx, CAnnSpectra *qr, vector<long> &indexlist, CMzFileReader *p, long totalSpec)
            : m_indexlist(indexlist), m_p(p) {
        m_ArchiveSpecNum = totalSpec;
        m_ptr = ptr;
        m_qr = qr;
        m_query_idx = query_idx;
    }

    void run() {
        // todo: we can use fragment index to accelerate it
        // todo: fragment index to be used here also!!!!!

        shared_ptr<ICHistdata> histscore = m_p->distributionPartialGeneral(1, true, m_indexlist, m_ptr); // to be fixed
        CTailEstimation cte(histscore, false);
        bool good = cte.buildLinearRegressionModel(true,
                                                   string("Q_") + to_string(m_query_idx) + "_partial_",
                                                   false, false);
        if (not good) {
            cout << "This is always good!!!! I feel " << good << endl << " " << m_query_idx << endl;
            return;
        }
        // todo: avoid some of the calculation
        float fraction = 1.0 * m_indexlist.size() / m_ArchiveSpecNum;
        for (int m = 0; m < m_qr->m_anns.size(); m++) {
            double x = m_qr->m_anns[m].dist;
            double dp = 1 - x * x / 2;
            m_qr->m_anns[m].setpvaluePartial(cte.getTopHitpvalue(dp, fraction, false));
        }
    }
};

class DPTaskFragIdx : public ICThreadTask {
public:

    struct Option {

        bool save_background_score;
        bool plot;
        bool verbose;
        int m_tol;

        Option() {
            save_background_score = false;
            verbose = false;
            plot = false;
            m_tol = 1;

        }
    };

private:
    uint16_t *m_ptr;
    long m_query_idx;  // used for generate the file name
    CAnnSpectra *m_qr;
    long m_ArchiveSpecNum;
    CFragIndex *m_cfi;
    bool m_plot;
    bool m_verbose;
    Option m_opt;
    vector<int> m_dpscores;
    SLinearRegressionModel m_model;
public:
    vector<int> getDPs(){
        return m_dpscores;
    }
    SLinearRegressionModel getLRModel(){
        return m_model;
    }
    DPTaskFragIdx(uint16_t *ptr, long query_idx, CAnnSpectra *qr, long totalSpec, CFragIndex *cfi, Option opt) {
        m_ArchiveSpecNum = totalSpec;
        m_ptr = ptr;
        m_qr = qr;
        m_query_idx = query_idx;
        m_cfi = cfi;
        m_plot = opt.plot;
//        if (m_query_idx==57) m_plot = true;
        m_verbose = opt.verbose;
        m_opt = opt;

    }

    void run() {
//        vector<int> dpscores;
        int tol = m_opt.m_tol;
        if(m_verbose){
            cout << "start dp score " << endl;
            CMzSpec(m_ptr,50).display();
        }
        m_cfi->dpscores(m_ptr, m_dpscores, true, tol);
        if(m_verbose)cout << "finished dp score" << endl;
        if (m_opt.save_background_score) {
            string filename = to_string("dpscore_10000_random_query_", "", m_query_idx,
                                        "_tol_", tol, ".txt");
            File::saveas(m_dpscores, filename, false);
        }
        if (m_plot) {
            if(m_verbose)cout << "m_plot" << endl;
            // so, it is because of thread safe....
            shared_ptr<ICId2RankInten> id2rank = m_cfi->getId2Inten();
            vector<double> scores(m_dpscores.begin(), m_dpscores.end());

            float maxScore = id2rank->getMaxScore();
            for (int i = 0; i < scores.size(); i++) {
                scores[i] /= maxScore;
            }
            vector<double> scores_sqrt_dp(scores.begin(), scores.end());
            for (int i = 0; i < scores_sqrt_dp.size(); i++) {
                scores_sqrt_dp[i] = sqrt(scores_sqrt_dp[i]);
            }

            shared_ptr<ICHistdata> histscoresqrt = ICMzFile::score2histogram(scores_sqrt_dp, 10000);
            CTailEstimation cte2(histscoresqrt, m_verbose);
            bool good2 = cte2.buildLinearRegressionModel(true,
                                                         string("Q_") + to_string(m_query_idx) +
                                                         "_sqrt_partial_",
                                                         false, m_plot);
            // scores...
            CVisual::SImageSize imsize;
            CVisual::gnuplot_histogram(scores, -1, 0.01, 1, 100, 40);
            string basename = to_string("Q_", "_", m_query_idx, "partial_tol", tol, "bgdphist.png");//
            // string("Q_") + to_string(m_query_idx) + "_partial_";
            CVisual::gnuplotWrapper info, info2;
            info.set_filename(basename + "_bgdphist.png")
                    .set_minmax(0.01, 1)
                    .set_width(imsize.m_wpx).set_height(imsize.m_hpx)
                    .set_xlabel("dot product")
                    .set_ylabel("frequency");
            CVisual::gnuplot_histogram_topng(scores, 100, 0.01, 1, info);

            // scores...
            CVisual::gnuplot_histogram(scores_sqrt_dp, -1, 0.01, 1, 100, 40);
            string basename_sqrt_dp = to_string("Q_", "_", m_query_idx,
                                                "sqrtdp", "partial_tol", tol, "bgdphist.png");
            //string("Q_") + to_string(m_query_idx) + "_sqrtdp_partial_";

            info2.set_filename(basename_sqrt_dp)
                    .set_minmax(0.01, 1)
                    .set_xlabel("sqrt_dp")
                    .set_ylabel("frequency")
                    .set_width(imsize.m_wpx)
                    .set_height(imsize.m_hpx);
            CVisual::gnuplot_histogram_topng(scores_sqrt_dp, 100, 0.01, 1, info2);
        }

        if(m_verbose){cout << "histogram" << endl;}
        shared_ptr<ICHistdata> histscore = ICMzFile::score2histogram(m_dpscores, 10000);

        if(m_verbose){cout << "tail estimation" << endl;}
        if(histscore==nullptr){
            return;
        }
        CTailEstimation cte(histscore, m_verbose);
        bool good = cte.buildLinearRegressionModel(true,
                                                   string("Q_") + to_string(m_query_idx) + "_partial_",
                                                   false, m_plot, &m_model);
        if (not good) {
            return;
        }
        float fraction = 1.0 * m_cfi->getSpecNum() / m_ArchiveSpecNum;

        for (int m = 0; m < m_qr->m_anns.size(); m++) {
            double x = m_qr->m_anns[m].dist;
            double dp = 1 - x * x / 2;
            m_qr->m_anns[m].setpvaluePartial(cte.getTopHitpvalue(dp, fraction, m_verbose));
        }
    }
};

class DPTaskGroup : public ICThreadTask {
    vector<ICThreadTask *> m_taskgroup;
public:
    DPTaskGroup(vector<DPTask> &vTask, int start, int batchsize) {
        for (int i = start; i < start + batchsize and i < vTask.size(); i++) {
            add(&vTask[i]);
        }
    }

    DPTaskGroup(vector<DPTaskFragIdx> &vTask, int start, int batchsize) {
        for (int i = start; i < start + batchsize and i < vTask.size(); i++) {
            add(&vTask[i]);
        }
    }

    void add(ICThreadTask *p) {
        m_taskgroup.push_back(p);
    }

    void run() {
        for_each(m_taskgroup.begin(), m_taskgroup.end(), [](ICThreadTask *x) { x->run(); });
    }
};

// Update the p-value of a list of query, and corresponding ANNs
// Input:
// indexListForQuery, (idx)
// query (Peaks)
//
// Output:
// p-values (updated in vqr)
//
// Note: multiple threads supported
void CPValueCalculator::run(vector<long> &idxListForQuery, ICQuery &query, CArxivSearchResult &vqr, bool plot, bool verbose,
                            vector<int> *ptrDP, SLinearRegressionModel *ptrLR , vector<vector<int> >*ptrDP_q, vector<SLinearRegressionModel> *ptrLR_q){
//    SimpleTimer st("p-value task");
    if(verbose) cout << "Using CPU calculate p-value: plot and verbose: " << std::boolalpha << "plot: " << plot << "  verbose: " << verbose << endl;
    vector<DPTaskFragIdx> vTaskfragidx;
    vector<ICThreadTask *> vthreadtask;
    DPTaskFragIdx::Option opt;
    opt.save_background_score = m_saveBackgroundScore;
    opt.verbose = verbose;
    opt.plot = plot;
    opt.m_tol = m_tol;

    // create list of tasks
    for (int j = 0; j < idxListForQuery.size(); j++) {
        uint16_t *ptr = query.getPtrUint16(j);
        // creating single jobs:
        // DPTaskFragIdx(uint16_t *ptr, long query_idx, CAnnSpectra *qr, long totalSpec, CFragIndex *cfi, Option opt)
        vTaskfragidx.emplace_back(ptr, idxListForQuery[j], vqr.get(j), m_total_spec_num, m_FragIdxPtr.get(), opt);
    }

    // put tasks as list of pointers
    for (int i = 0; i < vTaskfragidx.size(); i++) {
        vthreadtask.push_back(&vTaskfragidx[i]);
    }
//    for_each(vTaskfragidx.begin(), vTaskfragidx.end(),[](auto it))


    if(verbose) cout << "Final task: " << vthreadtask.size() << " Thread num is : " << m_threadnum << endl;
    CTaskPool::issueTasks(vthreadtask,verbose, false,m_threadnum, "p-value");


    // if only one query, keep the DP and the LR.
    if (idxListForQuery.size() == 1) {
        if (ptrDP != nullptr) {
            *ptrDP = vTaskfragidx[0].getDPs();
        }
        if (ptrLR != nullptr) {
            *ptrLR = vTaskfragidx[0].getLRModel();
        }
    }
    for(int i = 0; i < idxListForQuery.size(); i ++){
        if(ptrDP_q != nullptr){
            (*ptrDP_q)[i] = vTaskfragidx[i].getDPs();
        }
        if(ptrLR_q != nullptr){
            (*ptrLR_q)[i] = vTaskfragidx[i].getLRModel();
        }
    }
}

void CPValueCalculator::buildFragIndex(shared_ptr<ICMzFile> &csa, bool verbose) {
    const int PEAK_NUM_PER_SPEC = csa->getPeakNumPerSpec();
    m_FragIdxPtr = make_shared<CFragIndex>(PEAK_NUM_PER_SPEC);
    m_FragIdxPtr->buildIndex(csa, m_indexlist);
    if(verbose) cout << "frag index created" << endl;
}


CPValueCalculator::CPValueCalculator(long totalSpecNum, int threadnum, int seed, int tol, bool saveBackgroundScore,
                                     bool verbose) {
    m_saveBackgroundScore = saveBackgroundScore;
    m_tol = tol;
    if(verbose)cout << "[background spectra] Using seed " << seed << " for permutation" << endl;
    m_threadnum = threadnum;
    m_total_spec_num = totalSpecNum;
    if (threadnum == -1) {
        m_threadnum = getProperThreads() / 2 + 1;
        if(verbose)cout << "m_threadnum = " << m_threadnum << endl;
    }

    IndexListBuilder indexlist_builder;
    indexlist_builder.build(m_indexlist, totalSpecNum, seed);

    m_FragIdxPtr = nullptr;
}

CPValueCalculator::~CPValueCalculator() {

}

void CPValueCalculator::calcDotProducts(uint16_t *spec, vector<int> &score, bool normalize) {

    m_FragIdxPtr->dpscores(spec, score, normalize, m_tol);

}

CPvalueMultiCalculator::CPvalueMultiCalculator(int numCalc, long totalnum, int threadnum,
                                               shared_ptr<ICMzFile> &pScorer, int startseed, int tol) {
    m_tol = tol;
    for (int i = 0; i < numCalc; i++) {
        m_pvalueCalc.push_back(make_shared<CPValueCalculator>(totalnum, threadnum, i + startseed, m_tol,false,false));
        m_pvalueCalc[i]->buildFragIndex(pScorer, false);
    }
}
//
//vector<vector<int>> ptrDPs(numOfpValues, vector<int>());
//vector<vector<vector<int>>> ptrDP_all_pv_query(numOfpValues, vector<vector<int>>(idxlist.size(), vector<int>()));
//vector<SLinearRegressionModel> slRMs(numOfpValues, SLinearRegressionModel());
//vector<vector<SLinearRegressionModel>> slRMs_pv_query(numOfpValues, vector<SLinearRegressionModel>(idxlist.size(), SLinearRegressionModel()));
//

void CPvalueMultiCalculator::run(vector<long> &idxListForQuery, ICQuery &query, CArxivSearchResult &vqr, bool plot,
                                 bool verbose,vector<vector<int>> *ptrDPs,
                                 vector<SLinearRegressionModel> *ptrLRs,
                                 vector<vector<vector<int>>> *ptrDP_all_pv_query,
                                 vector<vector<SLinearRegressionModel>> *slRMs_pv_query) {
    if(verbose) cout << std::boolalpha << "plot: " << plot << "  verbose: " << verbose << endl;
    for(int i = 0 ; i < m_pvalueCalc.size(); i ++)
    {
        vector<int> *ptrDP = nullptr;
        SLinearRegressionModel * ptrLR = nullptr;
        if(ptrDPs != nullptr){
            ptrDP = &(*ptrDPs)[i];
        }
        if(ptrLRs != nullptr){
            ptrLR = &(*ptrLRs)[i];
        }
        vector<vector<int>> *ptrDP_q = nullptr;
        vector<SLinearRegressionModel> *ptrLR_q = nullptr;

        if(ptrDP_all_pv_query != nullptr){
            ptrDP_q = &(*ptrDP_all_pv_query)[i];
        }
        if(slRMs_pv_query != nullptr){
            ptrLR_q = &(*slRMs_pv_query)[i];
        }

        m_pvalueCalc[i]->run(idxListForQuery, query, vqr,plot, verbose,ptrDP, ptrLR, ptrDP_q, ptrLR_q);
    }
//    for (auto x: m_pvalueCalc) x->run(idxListForQuery, query, vqr, plot, verbose, ptrDP, ptrLR);
}

void IndexListBuilder::build(vector<long> &m_indexlist, long totalSpecNum, int seed) {
    // todo: the following index should be used with param: csa in run()
    m_indexlist.assign(totalSpecNum, 0);
    iota(m_indexlist.begin(), m_indexlist.end(), 0);
    if(m_fraction < 1 and m_fraction >0){
        gcc5shuffle(m_indexlist.begin(), m_indexlist.end(), mt19937_64(seed));
        long minNum = 5000;
        long newLen = m_fraction * totalSpecNum;
        if (newLen > m_MaxNum) {
            // newLen out of the window
            //--10000]---newLen-
            m_fraction = m_MaxNum*1.0/totalSpecNum;
        }
        if( newLen < minNum)
        {
            // ---newLen---[5000---
            m_fraction = min(m_MaxNum*1.0/totalSpecNum, 1.0);
        }
        m_indexlist.resize(long(m_fraction * totalSpecNum));
        cout << "[background spectra] Num: " << m_indexlist.size() << "\t fraction: " << m_fraction << endl;
    }else{
        cout << "[Error] invalid m_fraction, will use all index " << endl;
    }


}
