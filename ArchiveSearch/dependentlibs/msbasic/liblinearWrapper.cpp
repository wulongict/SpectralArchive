//
// Created by wulong on 1/18/18.
//

#ifdef LIBLINEAR_MT
#include "../../../External/liblinear-multicore-2.42/linear.h"
#else
#include "../../../External/liblinear-2.11/linear.h"
#endif
#include "liblinearWrapper.h"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>
#define INF HUGE_VAL
using namespace std;

void print_prob(problem *prob) {
    cout << "--Prob -----" << endl
         << "bias: " << prob->bias << " \n"
         << "l: " << prob->l << " \n"
         << "n: " << prob->n << " \n"
         << "x: " << prob->x << " \n"
         << "y: " << prob->y << endl;
    cout << "-- Prob ------" << endl;
}

void free_problem(problem *prob) {
    for(int i = 0; i < prob->l; i ++)
    {
        if (prob->x[i] != nullptr)        {
        free(prob->x[i]);
//        prob->x[i] == NULL;
        }  else   {
            cout << "empty pointer already prob->x[" << i << "]=" << prob->x[i] << endl;
        }
    }
    if (prob->x != nullptr)  {
    free(prob->x);
//        prob->x = NULL;
    } else  {
        cout << "empty pointer already " << prob->x << endl;
    }

    if (prob->y != nullptr)   {
    free(prob->y);
//        prob->y = NULL;
    } else {
        cout << "empty pointer already " << prob->y << endl;
     }

}

void initialize_problem(problem *prob, int size) {
    prob->bias = -1; //
    prob->l = size; // how many samples
    prob->y = (double *)malloc(prob->l * sizeof(double));

    if (prob->y == nullptr)  {
        cout << "invalid pointer prob->y " << __FILE__ << ":" << __LINE__ << endl;
        exit(0);
    }

    prob->x = (feature_node**)malloc ((prob->l) * sizeof(feature_node*));
    if (prob->x == nullptr) {
        cout << "invalid pointer prob->x " << __FILE__ << ": " << __LINE__ << endl;
        exit(0);
    }

    for (int i= 0; i < prob->l ; i ++)  {
        prob->x[i] = nullptr;
    }

}

void make_default_param_liblinear(parameter &param)  {// default values
//    param.solver_type = L2R_L2LOSS_SVC_DUAL; // -s 1
    param.solver_type = L2R_LR; // -s 0
    param.C = 1;
    param.eps = INF; // see setting below
    param.p = 0.1;
    param.nr_weight = 0;
    param.weight_label = nullptr;
    param.weight = nullptr;
    param.init_sol = nullptr;

#ifdef LIBLINEAR_MT
        // in this version, param has more members.
        // [YES] int solver_type;
        //
        //	/* these are for training only */
        //	[YES] double eps;             /* stopping tolerance */
        //	[YES] double C;
        //	int nr_thread;
        param.nr_thread = 10;
        //	[YES] int nr_weight;
        //	[YES] int *weight_label;
        //	[YES] double* weight;
        //	[YES] double p;
        //	double nu;
        param.nu = 0.5;
        //	[YES]double *init_sol;
        //	int regularize_bias;
        param.regularize_bias=1;


#endif

    if (param.eps == INF) {
        switch (param.solver_type) {
            case L2R_LR:
            case L2R_L2LOSS_SVC:
                param.eps = 0.01;
                break;
            case L2R_L2LOSS_SVR:
                param.eps = 0.001;
                break;
            case L2R_L2LOSS_SVC_DUAL:
            case L2R_L1LOSS_SVC_DUAL:
            case MCSVM_CS:
            case L2R_LR_DUAL:
                param.eps = 0.1;
                break;
            case L1R_L2LOSS_SVC:
            case L1R_LR:
                param.eps = 0.01;
                break;
            case L2R_L1LOSS_SVR_DUAL:
            case L2R_L2LOSS_SVR_DUAL:
                param.eps = 0.1;
                break;
        }
    }
}

