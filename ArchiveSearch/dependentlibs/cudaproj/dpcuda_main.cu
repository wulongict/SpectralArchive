#include "dpcuda.h"
#include "../../../librarymsms/Util.h"
#include <list>
#include "../../../librarymsms/ICQuery.h"
#include "../../../librarymsms/ProteomicsDataTypes.h"



void sort_for_top_1000(uint16_t *score, int specnum) {
    SimpleTimer st("sorting");
    vector<int> idx(specnum, 0);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [score](const int &a, const int &b) { return score[a] > score[b]; });
    st.stop();
    cout << "rank\tindex\tScore" << endl;
    for (int i = 0; i < 10; i++) {
        cout << i << "\t" << idx[i] << "\t" << score[idx[i]] << endl;
    }
}

void buketsort_for_top_1000(uint16_t *score, int specnum, vector<int> &topscore, int outputtopN) {
    SimpleTimer st("bucket sorting");
    vector<list<int>> idx(1 + UINT16_MAX, list<int>());
    for (int i = 0; i < specnum; i++) {
        idx[score[i]].push_back(i);
    }
    st.stop();
    int count = 0;

    cout << "rank\tindex\tScore" << endl;
    for (int i = 42925; count < outputtopN and i > 0; i--) {
        if (idx[i].size() != 0) {
            for (list<int>::iterator it = idx[i].begin(); count < outputtopN and it != idx[i].end(); it++) {
                cout <<count << "\t"<<  *it << "\t" << score[*it] << endl;
                count++;
                topscore.push_back(*it);
            }
        }
    }
    cout << "BUCKET for top 1000 DONE--------------------" << endl;

}

int main(void) {
    int topn = 50;
    int tol = 15;
    string mzfile = "/data/wulong/data/honeybee/paulmzXMLs.txt.mz";
    CUDAScore cudadp(mzfile);   //OK
    for (int j = 32; j <= 32; j += 32) {
        SimpleTimer st("topn=50 Query blocksize=" + to_string(j));
        for(int i = 1312081; i < 1312089; i ++)
        {
            cudadp.scoreAllVecForm(topn, tol, i, j,true);

        }
        st.stop();
        cout << "score 0 " << cudadp.getscore()[0] << endl;
    }
    vector<int> topscoreindex;
    int topn_rescore = 1;
    buketsort_for_top_1000(cudadp.getscore(), cudadp.getSpecNum(), topscoreindex, topn_rescore);
    vector<long> indexlist = {1312081,1342694,1306944,1340167,1478103,1304360,1428028,1345310,1347586,1498139};
    cout << "--- Start for the query of index list! --------" << endl;
    cudadp.queryFastOnIndex(50, 8, 1312081, 32, indexlist);
    cout << "------End of query of index list ------" << endl;
    for(int i = 0; i < indexlist.size(); i ++)    {
        cout << indexlist[i]       << " " << cudadp.getscore()[i] << endl;
    }

    return 0;
}