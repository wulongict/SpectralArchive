//
// Created by wulong on 2/28/19.
//



#include "ICMzFile.h"
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <cassert>
#include "ProteomicsDataTypes.h"
#include "Util.h"
#include "PeakList.h"
//#include "../ArchiveSearch/CTailEstimation.h"
#include "ICHistdata.h"

using namespace std;

long
ICMzFile::calculate_dot_product_sorted(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore) {
    if (mzTopN > getPeakNumPerSpec() or mzTopN < 1) {
        // Invalid value for TopN found, use Default:50 instead
        mzTopN = getPeakNumPerSpec();
    }
    long score = 0;
    uint16_t *x = getSpecBy(queryX), *y = getSpecBy(queryY);

    vector<int> yi(mzTopN), xi(mzTopN);
    iota(yi.begin(), yi.end(), 0);
    sort(yi.begin(), yi.end(), [y](const int &a, const int &b) { return y[a] < y[b]; });
    iota(xi.begin(), xi.end(), 0);
    sort(xi.begin(), xi.end(), [x](const int &a, const int &b) { return x[a] < x[b]; });
    int i = 0, j = 0;

    uint64_t recX = 0, recY = 0;

    while (i < mzTopN and j < mzTopN) {
        uint16_t alpha = x[xi[i]];
        uint16_t beta = y[yi[j]];
        if (alpha == 0) {
            i++;
            continue;
        }
        if (beta == 0) {
            j++;
            continue;
        }
        // we have minimal of x and minimal of y
        if (alpha < beta) //x is smaller
        {
            if (beta - alpha <= tol) {
                // beta in alpha tolerance
                // Match!!
                // ------alpha-tol----------alpha-----------alpha+tol---------
                //----------------|-----------*-----^--------|-----------------
                //----------------------------------|<-beta-------------------
                // Keep looking for bigger one! ?? NO;  we are not sure whetherr alpha is the largest one!
                //cout << score << " -<-- " << m_PeakNumPerSpectrum -xi[i] << " x " << m_PeakNumPerSpectrum-yi[j] << " --> -- ";
                score += prescore[xi[i]][yi[j]];
                //cout << score << endl;
                j++;
                i++;

//                    cout << score << " a=" << alpha << " b=" << beta << " xi="
//                        << xi[i] << " yi=" << yi[j] << " i="  << i <<" j=" << j << endl;
            } else {
                i++; // beta is out of windows of [alpha +/- tol]
            }
        } else// beta<=alpha
        {
            if (alpha - beta <= tol) {
                // Beta is in windows, but less than alpha
                // ------alpha-tol----------alpha-----------alpha+tol---------
                //----------------|---^--------*--------------|-----------------
                //--------------------|<-beta-------------------
                //cout << score << " -<-- " << m_PeakNumPerSpectrum -xi[i] << " x " << m_PeakNumPerSpectrum-yi[j] << " --> -- ";
                score += prescore[xi[i]][yi[j]];
                //cout << score << endl;
                j++;
                i++;
//                cout << score << " -<-- " << m_PeakNumPerSpectrum -xi[i] << " x " << m_PeakNumPerSpectrum-yi[j] << endl;
//                    cout << score << " a=" << alpha << " b=" << beta << " xi="
//                         << xi[i] << " yi=" << yi[j] << " i="  << i <<" j=" << j << endl;
            } else {
                j++; // alpha is too big, beta need to catch up;
            }
        }
    }

    return score;
}

long ICMzFile::calculate_dot_product_norec_short_if_else(long queryX, long queryY, int tol,
        int mzTopN, vector<vector<int>> &prescore) {
    if (mzTopN > getPeakNumPerSpec() or mzTopN < 1) {
        mzTopN = getPeakNumPerSpec();
    }
    long score = 0;
    uint16_t *x = getSpecBy(queryX), *y = getSpecBy(queryY);

    // todo: the algorithm can be improved!
    // faster: by put the score there
    for (int i = 0; i < mzTopN; i++) {
        if (x[i] == 0) {
            break;
        }

        for (int j = 0; j < mzTopN; j++) {
            if (y[j] == 0) {
                break;
            }
            if (x[i] > y[j] && x[i] - y[j] <= tol) // x[i] bigger make this shortter might be faster!todo
            {// add = Aug 1st 2018
                score += prescore[i][j];
                break;

            } else if (x[i] <= y[j] && y[j] - x[i] <= tol) { //add = Aug 1st 2018
                score += prescore[i][j];
                break;
            }
        }
    }
    return score;
}

long
ICMzFile::calculate_dot_product_norec(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore) {
    if (mzTopN > getPeakNumPerSpec() or mzTopN < 1) {
        // Invalid value for TopN found, use Default:50 instead
        mzTopN = getPeakNumPerSpec();
    }
    long score = 0;
    uint16_t *x = getSpecBy(queryX), *y = getSpecBy(queryY);

    // todo: the algorithm can be improved!
    // faster: by put the score there
    for (int i = 0; i < mzTopN; i++) {
        if (x[i] == 0) {
            break;
        }
        for (int j = 0; j < mzTopN; j++) {
            if (y[j] == 0) {
                break;
            }
            if (x[i] > y[j]) // x[i] bigger make this shortter might be faster!todo
            {
                if (x[i] - y[j] <= tol) { // add = Aug 1st 2018
                    score += prescore[i][j];
                    break;
                }
            } else {
                if (y[j] - x[i] <= tol) { // add = Aug 1st 2018
                    score += prescore[i][j];
                    break;
                }
            }
        }
    }
    return score;
}

long ICMzFile::calculate_dot_product(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore) {
    if (mzTopN > getPeakNumPerSpec() or mzTopN < 1) {
        // Invalid value for TopN found, use Default:50 instead
        mzTopN = getPeakNumPerSpec();
    }
    long score = 0;
    uint16_t *x = getSpecBy(queryX), *y = getSpecBy(queryY);

    uint64_t rec = 0;

    // todo: the algorithm can be improved!
    for (int i = 0; i < mzTopN; i++) {
        if (x[i] == 0) break;
        for (int j = 0; j < mzTopN; j++) {
            if (y[j] == 0) {
                break;
            }

            if (rec & (1L << j)) {
                continue;
            }
            if (x[i] > y[j]) // x[i] bigger
            {
                if (x[i] - y[j] < tol) {// add = on Aug 1st 2018; to be done on GPU todo:
                    score += prescore[i][j];
                    rec = (rec | (1L << j));
                    break;
                }
            } else {
                if (y[j] - x[i] < tol) { // add = on Aug 1st 2018; to be done on GPU todo:
                    score += prescore[i][j];
                    rec = (rec | (1L << j));
                    break;
                }
            }
        }
    }
    return score;
}

long ICMzFile::calculate_dot_product_with_hist(long queryX, long queryY, int tol, int mzTopN, vector<vector<int>> &prescore, vector<double> &diff) {

    uint16_t *x = getSpecBy(queryX), *y = getSpecBy(queryY);
    return calculate_dot_product_with_hist(tol, mzTopN, prescore, diff, x, y, false);
}
struct weightDiff
{
    double diff;
    double weight;
    double mz1;
    double mz2;
    double rank1;
    double rank2;
    weightDiff(double d, double w, double m1, double m2, double r1, double r2):diff(d),weight(w), mz1(m1), mz2(m2), rank1(r1), rank2(r2){
//        print();
    }
    void print() const{
        print_on_fixed_width(10,mz1, rank1, mz2, rank2, diff, weight,"\n");
    }
};

long ICMzFile::calculate_dot_product_with_hist(int tol, int mzTopN, vector<vector<int>> &prescore, vector<double> &diff,
                                               uint16_t *x, uint16_t *y, bool print_matched_pks=false) const {
    const double ratio = 2000.0/65535;

    if (mzTopN > getPeakNumPerSpec() or mzTopN < 1) {
        mzTopN = getPeakNumPerSpec();
    }
    long score = 0;

    uint64_t rec = 0, recx=0;

    string matched_peaks ="matched peaks\n\t";

    // todo: the algorithm can be improved!
    int cnt_matched = 0;
    for (int i = 0; i < mzTopN; i++) {

        if (x[i] == 0) break;
        for (int j = 0; j < mzTopN; j++) {
            if (y[j] == 0) {
                break;
            }

            if (rec & (1L << j)) {
                continue;
            }

            if (x[i] > y[j]) // x[i] bigger
            {
                if (x[i] - y[j] <= tol) {// add = on Aug 1st 2018; to be done on GPU todo:
                   matched_peaks =  to_string(matched_peaks,"\t",50-i,x[i]*ratio,50-j,y[j]*ratio,"--",i,j,"\n");
                    score += prescore[i][j];
                    rec = (rec | (1L << j));
                    recx = (recx| (1L<<i));
                    cnt_matched ++;
                    break;
                }
            } else {
                if (y[j] - x[i] <= tol) { // add = on Aug 1st 2018; to be done on GPU todo:
                    matched_peaks = to_string(matched_peaks,"\t",50-i,x[i]*ratio,50-j,y[j]*ratio,"--",i,j,"\n");
                    score += prescore[i][j];
                    rec = (rec | (1L << j));
                    recx = (recx| (1L<<i));
                    cnt_matched ++;
                    break;
                }
            }
        }
    }
    if(print_matched_pks) cout << matched_peaks << endl;

    vector<int> x_r, y_r;
    vector<int> x_w, y_w;

    const int NUM_PEAK = 50;
    for(int i = 0; i < mzTopN; i ++)  {
        if (rec & (1L << i)) {
           // continue;
        } else if(y[i]>0){
            y_r.push_back(y[i]);
            y_w.push_back(NUM_PEAK-i);
        }

        if (recx & (1L << i)) {
         //   continue;
        } else if(x[i]>0){
            x_r.push_back(x[i]);
            x_w.push_back(NUM_PEAK-i);
        }
    }

    cout << "peaks remaining: " << x_r.size() << "\t" << y_r.size() << endl;
    cout << "peaks matched: " << cnt_matched << endl;
    vector<weightDiff> wdiff;
    if(!x_r.empty() and y_r.size()!=0)  {

        for(int i = 0; i < x_r.size(); i ++)
        {
            for(int j = 0; j < y_r.size(); j ++)  {
                double d = (x_r[i]-y_r[j])*ratio;
                diff.push_back(d);
                wdiff.emplace_back(weightDiff(d,x_w[i]*y_w[j],x_r[i]*ratio,y_r[j]*ratio,x_w[i], y_w[j]));
            }
        }
        sort(diff.begin(), diff.end());
        sort(wdiff.begin(), wdiff.end(),[](const weightDiff& a, const weightDiff &b){return a.diff < b.diff;});


        for(int i = 0; i < wdiff.size()-1; i ++)  {
            if(i>0 and fabs(wdiff[i].diff-wdiff[i-1].diff)<EPSILON) continue;
            int counts = 1;
            double sum = wdiff[i].weight;
            weightDiff &w = wdiff[i];
            string info = to_string("","\t","\tmz1\trank1\t -- \tmz2\trank2\tdiff\tweight\n",w.mz1, w.rank1," -- ",w.mz2, w.rank2, w.diff , w.weight, "\n");

            int j = i +1;
            while(j < wdiff.size())  {
                if(wdiff[j].diff-wdiff[i].diff>0.05) break;
                else  {
                    counts ++;
                    sum += wdiff[j].weight;
                    weightDiff &w1 = wdiff[j];
                    info = to_string(info,"\t",w1.mz1, w1.rank1," -- ",w1.mz2, w1.rank2, w1.diff , w1.weight, "\n");
                }
                j++;
            }

            if(counts>2 and sum > 3000 )  {
                string tmp = to_string("range["," ",diff[i],'-',diff[j-1], "]",
                        "FREQ:", counts, "OLD:", score, "NEW", sum, "Total:",
                        score+sum, "UP", sum/score);

                CANSIConsole a;
                if(counts>=5) cout << a.getColorStr(tmp, CANSIConsole::BLUE) << endl;
                else{
                    cout << tmp  << endl;
                }


                if(print_matched_pks)cout << info << endl;
            }
        }
    }
    return score;
}

// use this one
void ICMzFile::getQueryVecByMZArrayIndex(int queryindex, vector<float> &v, bool useFlankingBins) {
    const int PEAKNUM_PER_SEPC = 50; // todo: should not be fixed!
    const int dim = 4096; // todo: should not be fixed!

    uint16_t *p = getSpecBy(queryindex);
	assert(p != nullptr);

    vector<uint16_t> queryspec(p, p + PEAKNUM_PER_SEPC);

    vector<double> mz(PEAKNUM_PER_SEPC, 0), intensity(PEAKNUM_PER_SEPC, 0);
    const double MAX_PEAK_NUM = UINT16_MAX;
    const double MAX_PEAK_MZ = 2000; // todo: should not be fixed!
    for (int i = 0; i < PEAKNUM_PER_SEPC; i++) {
        if (queryspec[i] == 0) {
            mz.resize(i);
            intensity.resize(i);
            break;
        }
        mz[i] = queryspec[i];
        mz[i] /= MAX_PEAK_NUM;
        mz[i] *= MAX_PEAK_MZ;
        intensity[i] = PEAKNUM_PER_SEPC - i;
    }

    BinningPeakList bpl(mz, intensity, useFlankingBins);
    getFloatVecPaddingZeros(&bpl, dim, false, v);
}

void ICMzFile::get_vector_form(long idx, int tol, vector<int> &vecform) {
    uint16_t *x = getSpecBy(idx);
	get_vector_form(x, tol,  vecform);
}


void ICMzFile::get_vector_form(uint16_t * x, int tol, vector<int> &vecform) const {
	const int PEAKNUM_PER_SEPC = getPeakNumPerSpec();
	int veclen = UINT16_MAX + 1;
	vecform.assign(veclen, 0);

	for (int i = 0; i < PEAKNUM_PER_SEPC; i++) {
		if (x[i] == 0) break;
		int j = x[i] >= tol ? x[i] - tol : 0;
		int max_j = x[i] >= veclen - tol ? veclen - tol : x[i] + tol;
		while (j <= max_j)		{
			if (j != 0 and vecform[j]==0) vecform[j] += (PEAKNUM_PER_SEPC - i);  // fixed another diff
			j++;
		}
	}
}

// this calculation is not smart
long ICMzFile::calculate_dot_product_with_vecfrom(long queryX, vector<int> &vecformY, int mzTopN, bool debug,
                                                  long debug_index) {
        long s = 0;
    try {
        const int PEAKNUM_PER_SPEC = getPeakNumPerSpec();
        uint16_t *x = getSpecBy(queryX);
        // we will use index 1-50 to get value, so we need 51 values. as C++ array index is 0-based
        uint16_t used[51]={0};
        for (int i = 0; i < mzTopN; i++) {
            if(x[i] == 0) { // last peak, no more peaks.
                break;
            }
            if (vecformY[x[i]] == 0 or used[vecformY[x[i]]] == 1) {
                // peaks dost not matched or already matched with another peak.
                continue;
            }
            // new unmatched peak get matched, set the flag on intensity value.
            used[vecformY[x[i]]] = 1;

            // get the multiplication of two intensities.
            s += (PEAKNUM_PER_SPEC - i) * vecformY[x[i]];
            if(debug and queryX == debug_index)  {
                cout <<"peak i = " << i << ": " <<  (PEAKNUM_PER_SPEC - i) << " x " <<  vecformY[x[i]] << "-->" <<s << endl;
            }
        }

    }
    catch(const exception &ex){
        cout << "Error: what  " << ex.what() << endl;
    }
    catch(...)
    {
        cout << "Error: " << queryX << " -- " << debug << endl;
    }

    return s;
}

float *ICMzFile::buildNormalizedQueries(int dim, const vector<long> &candidates, bool useFlankingBins, int numQuery) {
    float *vquery = new float[numQuery * dim];
	assert(vquery != nullptr);
    for (int j = 0; j < candidates.size(); j++) {
        vector<float> tmp_query;
        getQueryVecByMZArrayIndex(candidates[j], tmp_query, useFlankingBins);
        
        L2Normalization(tmp_query, dim); // todo: which function to use?
        copy(tmp_query.begin(), tmp_query.end(), vquery + dim * j);  // copy from tmp to vquery
    }
    return vquery;
}

void ICMzFile::get_prescoreMatrix(vector<vector<int>> &prescorematrix) const {
    const int PeakNum = getPeakNumPerSpec();
    prescorematrix.assign(PeakNum, vector<int>(PeakNum, 0));
    for (int i = 0; i < PeakNum; i++) {
        for (int j = 0; j < PeakNum; j++) {
            prescorematrix[i][j] = (PeakNum - i) * (PeakNum - j);
        }
    }
}

void ICMzFile::displayspec(vector<uint16_t> &spec) {
    const double MAX_SPEC_MZ = 2000;
    const int PeakNum = spec.size();
    vector<double> mz(PeakNum, 0), intensity(PeakNum, 0);
    for (int i = 0; i < PeakNum; i++) {
        mz[i] = spec[i] * 1.0 / UINT16_MAX * MAX_SPEC_MZ;
        intensity[i] = PeakNum - i;
    }
    vector<int> idx(PeakNum, 0);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&mz](const int &x, const int &y) { return mz[x] < mz[y]; });
    cout << "id\tmz\tintensity\n-------------------------" << endl;
    for (int i = 0; i < PeakNum; i++) {
        cout << i << "\t" << mz[idx[i]] << "\t" << intensity[idx[i]] << endl;
    }
    cout << "--------------------------" << endl;
}

float ICMzFile::getSquaredNorm(long queryindex) {
    uint16_t *p = getSpecBy(queryindex);
    return getSquaredNorm(p);
}

float ICMzFile::getSquaredNorm(uint16_t *p) {
    const int PEAKNUM_PER_SEPC = getPeakNumPerSpec();
    float norm = ICMzFile::getSquaredNorm(p, PEAKNUM_PER_SEPC);
    return norm;
}

float ICMzFile::getSquaredNorm(const uint16_t *p, const int PEAKNUM_PER_SEPC)  {
    float norm = 0;
    for (int i = 0; i < PEAKNUM_PER_SEPC; i++) {
        if (p[i] == 0) {
            break;
        } else {
            norm += (PEAKNUM_PER_SEPC - i) * (PEAKNUM_PER_SEPC - i);
        }
    }
    return norm;
}

vector<int> ICMzFile::distributionPartial(int tol, bool normalize, vector<long> &indexlist,
                                           uint16_t *queryspec)
{
    vector<int> scores;
    scorePartiallyWithVecForm(getPeakNumPerSpec(), tol, 32, normalize, indexlist, queryspec, scores);
    string filename = "all_scores.txt";
    File::saveas(scores, filename,false);

    return  ICMzFile::score_to_histogram(scores);
}

shared_ptr<ICHistdata> ICMzFile::distributionPartialGeneral(int tol, bool normalize,
                                                            vector<long> &indexlist, uint16_t *queryspec) {
    vector<int> scores;
    scorePartiallyWithVecForm(50, tol, 32, normalize, indexlist, queryspec, scores);
    string filename = "all_scores.txt";
    File::saveas(scores, filename,false);
    return ICMzFile::score2histogram(scores,10000);

}

vector<int> ICMzFile::distributionAll(int tol, long queryindex, bool normalize) {
    vector<long> all_ids(getSpecNum(),0);
    iota(all_ids.begin(), all_ids.end(),0);
    uint16_t *queryspec=getSpecBy(queryindex);
    return distributionPartial(tol, normalize,all_ids, queryspec);

}

vector<int> ICMzFile::score_to_histogram(vector<int> &scores) {
    const int MAX_TOP50_COS = 42925;
    int left = 0, right = MAX_TOP50_COS;
    auto it = find_if(scores.begin(), scores.end(), [&](const int x){return x<left or x > right;});
    if(it != scores.end() )    {
        throw runtime_error("out of range: 0 ~ 42925");
    }
    vector<int> histogram(MAX_TOP50_COS + 1,0);

    for_each(scores.begin(), scores.end(),[&](const int &x){ histogram[x]++;  });
    return histogram;
}



void ICMzFile::getCompactForm(const double *mz, const double *intensity, vector<uint16_t> &newspec) const {
    const int nPeakNum = getPeakNumPerSpec();
    const double MAX_PEAK_MZ = 2000.0;
    vector<int> idx(nPeakNum);
    iota(idx.begin(), idx.end(), 0);
    getsortindex(idx, intensity);
    for (int i = 0; i < nPeakNum; i++) {
        double t = mz[idx[i]] * 1.0 / MAX_PEAK_MZ * UINT16_MAX;
        newspec[i] = t > UINT16_MAX ? UINT16_MAX : t;
    }
}

void ICMzFile::getsortindex(vector<int> &idx, const double *intensity) const {
    sort(idx.begin(), idx.end(), [intensity]
                 (const int &a, const int &b) -> bool {
             return intensity[a] > intensity[b];
         }
    );

}

// choose normalize the dot product is better.
void ICMzFile::scorePartiallyWithVecForm(int mzTopN, int tol, int blockSize, bool normalize, vector<long> &indexlist,
                                         uint16_t *queryspec, vector<int> &scores) {
    if(indexlist.empty())
    {
        cout << "Empty index list in " << __FUNCTION__ << endl;
    }
    calcDotProduct(mzTopN, tol, queryspec, blockSize, indexlist, scores);

    // here is normalized the score into the range of 42925.
    if(normalize)  {
        double querynorm = getSquaredNorm(queryspec);
        const int MAX_SCORE = 42925;

        for (long i = 0; i < scores.size(); i++) {
            double s = scores[i];
            if (s > EPSILON) {
                s /= sqrt(getSquaredNorm(indexlist[i]));
                s /= sqrt(querynorm);
            }

            s *= MAX_SCORE;
            if (s < 0 or int(s) > MAX_SCORE ) {
                cout << "Normalize Partial Scores: Invalid score s: " << s << " on index: " << indexlist[i] << endl;
                cout << "origin: " << scores[i] << " norm: " << querynorm << " & " << getSquaredNorm(indexlist[i]) << endl;

            }
            scores[i] = int(s) > MAX_SCORE ? MAX_SCORE : int(s);
        }
    }
}

void ICMzFile::dist(long query_index, vector<long> &ind, int tol, vector<float> &distances,vector<int> &dpscore) {
    if(ind.empty())    {
        cout << "Empty index list in " << __FUNCTION__ << endl;
    }
    uint16_t * queryspec = getSpecBy(query_index);
    distOnSpec(queryspec, ind, tol, distances,dpscore);
}

// the dot product is normalized, then the score is converted to L2 distance.
void ICMzFile::distOnSpec(uint16_t *queryspec, vector<long> &ind, int tol, vector<float> &dist,vector<int> &dpscore) {
    if(ind.empty())
    {
//        cout << "Empty index list in " << __FUNCTION__ << endl;
        return;
    }
    //        const double NORMALIZER_TOP50_COS = 42925;
//        int tol = 15;    // mass tolerance: 0.45 Th
    // : unit bins; 0 ~ 2000 is divided into 65535 bins, so each bin is : 0.03 Th
    // it is not good enough for HCD, but already good enough; 15 bins is about 0.45 Th which is low mass accuracy
    int useTopNPeak = 50;
    int blockSize = 32;

    vector<long> indexlist(ind.begin(), ind.end()); // todo: why I make a copy of the list ???

    calcDotProduct(useTopNPeak, tol, queryspec, blockSize, indexlist, dpscore);

    float querySquaredNorm = getSquaredNorm(queryspec);
    dist.assign(indexlist.size(), 0);

    for (int i = 0; i < indexlist.size(); i++) {
        float candSquaredNorm = getSquaredNorm(indexlist[i]);
        dist[i] = dpscore[i];

        if (candSquaredNorm > EPSILON and querySquaredNorm > EPSILON) {
            // dp normalized as cosine score.
            dist[i] /= sqrt(candSquaredNorm * querySquaredNorm);
        }
        dist[i] = sqrt(2.0) * sqrt(1.0 - dist[i] < EPSILON ? 0 : 1.0 - dist[i]);
    }

}

int ICMzFile::getPeakNum(long queryindex) {
    return getMzSpec(queryindex).getPeakNum();
}

CMzSpec ICMzFile::getMzSpec(long queryindex) {
    return CMzSpec(getSpecBy(queryindex), getPeakNumPerSpec());
}

string ICMzFile::getClassName() {return "ICMzFile";}

shared_ptr<ICMzFile> createScorer(bool use_gpu, string mzXMLListFileName) {
    if(use_gpu)
    {
#ifdef __CUDACC__
        cout << "defined __CUDACC__  " << __CUDACC__ << endl;
        return  make_shared<CUDAScore>(mzXMLListFileName + ".mz");
#endif
    }
    return nullptr;
}

shared_ptr<ICMzFile> ICMzFactory::getMzFileObject() {
    assert(m_csa!=nullptr);
//    cout << "pointer of csa " << m_csa << endl;
//    int num = m_csa->getPeakNum(0L);
//    float norm = m_csa->getSquaredNorm(0L);
//    cout << "initialized m_csa: num =" << num << "\t norm="<< norm << endl;
    return m_csa;
}

string ICMzFactory::getmzfilename() {return mzfilename;}
