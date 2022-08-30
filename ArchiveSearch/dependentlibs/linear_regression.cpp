//
// Created by wulong on 8/25/22.
//

#include <iostream>
#include <vector>
#include <cmath>

//
//       S_Y S_{XX}-S_X S_{XY}
// b   = ---------------------------------
//       n  S_{XX} -  S_X S_X
//
//        n S_{XY} - S_Y S_X
// k   =   ---------------------------
//        n  S_{XX} -  S_X S_X

// 1D linear regression problem.
// y = kx +b
bool simple_1D_linear_regression(std::vector<double> x, std::vector<double> y, double &k, double &b, double &r2){
    bool ret = false;
    int n = x.size();
    double sx = 0, sy=0, sxx=0, sxy=0;
    for(int i = 0; i < x.size(); i ++){
        sx += x[i];
        sy += y[i];
        sxy += x[i]*y[i];
        sxx += x[i] * x[i];
    }

    double denominator = n * sxx - sx * sx;
    double EPSILON = 1e-8;
    if(fabs(denominator)<EPSILON){
        ret = false;
    }else{
        b = (sy * sxx - sx * sxy) / denominator;
        k = (n * sxy - sy * sx)/ denominator;
        // calculate r2
        double y_bar = sy / n;
        double SSR = 0, SST = 0;
        for(int i = 0; i < x.size(); i ++){
            double t = y[i]-(k * x[i] + b), tt = y[i]-y_bar;
            SSR += t*t;
            SST += tt * tt;
        }
        std::cout << SSR << "\t" << SST << std::endl;
        r2 = 1- SSR/SST; // is it possible that SST = 0?
        ret = true;
    }

    return ret;

}


#include <algorithm>
int main(int argc, char *argv[]){
    std::vector<double> x={1,2,3};
    std::vector<double> y={4.00,7.0,10.001};//{1.0000001,1.0,1.000001};//
    // rescale y
    double sy = 0;
    for(int i = 0; i < y.size(); i ++){
        sy += y[i];
    }
    double ybar = sy/ y.size();
    for(int i = 0; i < y.size(); i ++){
        y[i] -= ybar;
        std::cout << y[i] << std::endl;
    }

    auto elm = std::minmax_element(y.begin(), y.end());
    std::cout << *elm.first << " from \t to " << *elm.second << std::endl;
    const auto [min, max] = std::minmax_element(begin(y), end(y));
    std::cout << *min << "\t" << *max << std::endl;

    double miny = *min;
    double maxy = *max;


    for(int i = 0; i < y.size(); i ++){
        std::cout <<"yi, min, max " << y[i] << "\t" <<  *(elm.first) << " ? changed? " << *(elm.second) << std::endl;
        std::cout << (y[i]-*(elm.first))<< "\t" << (*(elm.second) - *(elm.first)) << std::endl;
        std::cout << (y[i]-*(elm.first))<< "\t" << (*(elm.second) - *(elm.first)) << std::endl;
        std::cout << (y[i]-*(elm.first))<< "\t" << (*(elm.second) - *(elm.first)) << std::endl;
        y[i] = (y[i]-miny)/(maxy - miny);
        std::cout <<i << "\t" <<  y[i] << std::endl;
        std::cout <<"yi, min, max " << y[i] << "\t" <<  *(elm.first) << " ? changed? " << *(elm.second) << std::endl;
    }

//    y={4.00,7.0,10.001};
    double k, b, r2;
    simple_1D_linear_regression(x, y, k, b, r2);
    std::cout << k << "\t"<< b  << "\t"<< r2 << std::endl;



    return 0;
}
