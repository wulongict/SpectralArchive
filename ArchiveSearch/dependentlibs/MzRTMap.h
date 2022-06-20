//
// Created by wulong on 9/1/20.
//

#ifndef MYTOOL_MZRTMAP_H
#define MYTOOL_MZRTMAP_H
#include <vector>
#include <string>

class DataFile;
class CMzRtMap{
public:
    struct projectionParam{
        double minMz;
        double maxMz;
        double minRT;
        double maxRT;
        int xPixel;
        int yPixel;
        projectionParam(double _minMz, double _maxMz, double _minRT, double _maxRT, int _xPixel, int _yPixel){
            minMz = _minMz;
            maxMz = _maxMz;
            minRT = _minRT;
            maxRT = _maxRT;
            xPixel = _xPixel;
            yPixel = _yPixel;
        }
    };
    struct peak{
        double mz;
        double rt;
        double intensity;
        peak(double _mz, double _rt, double _intensity):mz(_mz),intensity(_intensity),rt(_rt){}
    };
private:
//    projectionParam m_paramProj;
    std::vector<peak> x;
    std::vector<std::vector<double> > mzrtMap;
public:
    CMzRtMap(std::string filename);
    CMzRtMap(DataFile &df);
    void projection(projectionParam &p);
    void gnuplotView(projectionParam &p, bool svg);
    void projectView(projectionParam &p,bool svg);
    std::string toJsonStr();

    void initializeMap(DataFile &df) ;
};
#endif //MYTOOL_MZRTMAP_H
