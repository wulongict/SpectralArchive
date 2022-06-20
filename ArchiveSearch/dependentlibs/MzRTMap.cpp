//
// Created by wulong on 9/1/20.
//


#include "../../librarymsms/ProteomicsDataTypes.h"
#include "../../librarymsms/Util.h"
#include "../../External/gnuplot-iostream/gnuplot-iostream.h"
#include "MzRTMap.h"
#include <iostream>
using namespace std;

CMzRtMap::CMzRtMap(DataFile &df) {
    initializeMap(df);
}

void CMzRtMap::initializeMap(DataFile &df) {
    double maxRT = df.getSpectrum(df.getSpectrumNum() - 1)->m_rt_in_sec;
    double maxMz = 2000;
//        vector<CMzRtMap::peak> x;
    cout << "maxRT is " << maxRT << endl;
    int total = 0;
    for (int i = 0; i < df.getSpectrumNum(); i++) {
        CSpectrum *spec = df.getSpectrum(i);
        if (spec->getMSLevel() != 1) {
            continue;
        }
// get MS1 spectrum
        double rt = spec->m_rt_in_sec;
        int peaknum = spec->getPeakNum();
        total += peaknum;
        for (int j = 0; j < peaknum; j++) {
            double mz, intensity;
            spec->getOnePeak(mz, intensity, j);
            x.emplace_back(mz, rt, intensity);
        }
    }
    maxRT = floor(maxRT / 100) * 100 + 100;
    cout << "total peak count " << x.size() << "and maxRT is " << maxRT << endl;
}

CMzRtMap::CMzRtMap(string filename) {
    DataFile df(filename,0,-1);
    initializeMap(df);
//    double maxRT = df.getSpectrum(df.getSpectrumNum()-1)->m_rt_in_sec;
//    double maxMz = 2000;
////        vector<CMzRtMap::peak> x;
//    cout << "maxRT is " << maxRT << endl;
//    int total = 0;
//    for(int i = 0; i < df.getSpectrumNum(); i ++)
//    {
//        CSpectrum * spec = df.getSpectrum(i);
//        if(spec->getMSLevel() != 1){
//            continue;
//        }
//        // get MS1 spectrum
//        double rt = spec->m_rt_in_sec;
//        int peaknum = spec->getPeakNum();
//        total += peaknum;
//        for(int j =0; j < peaknum; j ++)
//        {
//            double mz, intensity;
//            spec->getOnePeak(mz, intensity, j);
//            x.emplace_back(mz,rt,intensity);
//        }
//    }
//    maxRT = floor(maxRT/100)*100+100;
//    cout << "total peak count " << x.size() << "and maxRT is " << maxRT<< endl;
}

void CMzRtMap::projection(CMzRtMap::projectionParam &p) {
    SimpleTimer st("projection");
    mzrtMap.resize(0);
    vector<vector<peak*>> pic(p.xPixel*p.yPixel,vector<peak*>());
    cout << "x size: " << x.size() << endl;
    for(int i = 0; i < x.size(); i ++)
    {
        if(x[i].rt<p.minRT) continue;
        if(x[i].rt>p.maxRT) break;
        if(x[i].mz<p.minMz or x[i].mz >p.maxMz) continue;

        // the mz and rt are in the window
        int xPos = floor((x[i].rt-p.minRT)/(p.maxRT-p.minRT)*p.xPixel);
        int yPos = floor((x[i].mz-p.minMz)/(p.maxMz-p.minMz)*p.yPixel);
//        cout << "found dot " << i  << endl;
        pic[xPos * p.yPixel + yPos].push_back(&x[i]);// += x[i].intensity;

    }
    vector<peak> outpic;

    for(int i = 0; i < p.xPixel; i ++)
    {
        for(int j = 0; j < p.yPixel; j ++)
        {
            int k = i * p.yPixel + j;
            if(pic[k].empty()) continue;
            peak tmp(0,0,0);
            for(int l = 0; l < pic[k].size(); l ++)
            {
                tmp.mz += pic[k][l]->mz;
                tmp.intensity += pic[k][l]->intensity;
                tmp.rt += pic[k][l]->rt;
            }
            tmp.mz = tmp.mz/pic[k].size();
            tmp.rt = tmp.rt/pic[k].size();

            outpic.push_back(tmp);
            mzrtMap.push_back(vector<double>({tmp.rt,tmp.mz,tmp.intensity,log(tmp.intensity+1)}));
        }
    }
    cout << "pic size: " << outpic.size() << endl;
    cout << mzrtMap.size() << endl;
}

void CMzRtMap::gnuplotView(CMzRtMap::projectionParam &p, bool svg) {
    SimpleTimer st("gnuplot");
    Gnuplot gnuplot;
    if(svg){
        gnuplot << "set terminal svg" << endl;
    }
    gnuplot << "set style fill transparent solid 0.2 noborder" << endl;
    gnuplot << "set palette defined (0 \"white\", 0.2 \"light-gray\", 0.4 \"blue\", 0.6 \"green\", 0.8 \"yellow\", 1 \"red\" )" << endl;
    gnuplot << "set logscale cb" << endl; // set log scale color bar
    gnuplot << "set format cb \"10^{%L}\"" << endl; // set the color bar label in 10^n format.
    gnuplot << "set xrange [" << p.minRT <<":" << p.maxRT <<"]" << endl;
    gnuplot << "set yrange [" << p.minMz <<":" << p.maxMz <<"]" << endl;
    gnuplot << "set xlabel \'RT(sec)\'"  << endl;
    gnuplot << "set ylabel \'Mz(Th)\'" << endl;
    // square: pt 5; circle: pt 7; triangle: pt 9
    gnuplot << "plot '-' using 1:2:3 with points pt 5 ps 0.4 palette notitle" << endl;
    gnuplot.send1d(mzrtMap);
}

void CMzRtMap::projectView(CMzRtMap::projectionParam &p, bool svg) {


    projection(p);
    // visualization of outpic
    gnuplotView(p,svg);
}

std::string CMzRtMap::toJsonStr() {
    ostringstream  oss;
    if(mzrtMap.empty()){
        cout << "no map available" << endl;
    }
    else{
        oss << "[" << endl;
        for(int i = 0; i < mzrtMap.size(); i ++)
        {
            if(i>0) oss << "," ;
            vector<double> &t = mzrtMap[i];
            oss << "{\"rt\":" << t[0] << ", \"mz\": " << t[1] << ", \"inten\": " << t[2] << ", \"loginten\": " << t[3] << "}" << endl;
        }
        oss << "]" << endl;
    }
    return oss.str();
}
