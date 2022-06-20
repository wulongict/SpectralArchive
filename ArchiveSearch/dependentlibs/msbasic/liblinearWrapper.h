//
// Created by wulong on 1/18/18.
//

#ifndef MYTOOL_LIBLINEARWRAPPER_H
#define MYTOOL_LIBLINEARWRAPPER_H

class problem;
class parameter;

void print_prob(problem *prob);
void free_problem(problem *prob) ;

void initialize_problem(problem *prob, int size);

void make_default_param_liblinear(parameter &param);

//// todo: fix the bug later.
//void do_predict(string model_file, vector<double> &probs, vector<int> &labels)
//{
//    struct feature_node *x;
//    int max_nr_attr = 64;
//
//    vector<CPeakPair> samples;
//    model * model_=load_model(model_file.c_str());
//
//    int nr_class=get_nr_class(model_);
//    double *prob_estimates=NULL;
//    prob_estimates = (double *) malloc(nr_class*sizeof(double));
//
//    if(!check_probability_model(model_))
//    {
//        fprintf(stderr, "probability output is only supported for logistic regression\n");
//        exit(1);
//    }
//
//
//    for(auto eachrow : samples)
//    {
//        int i = 0;
//        double predict_label = 0;
//        int inst_max_index = 0; // strtol gives 0 if wrong format
//
//        for(auto eachitem : eachrow)
//        {
//            if(i>=max_nr_attr-2)	// need one more for index = -1
//            {
//                max_nr_attr *= 2;
//                x = (feature_node *) realloc(x,max_nr_attr*sizeof(struct feature_node));
//            }
//            x[i].index = eachitem.index;
//            x[i].value = eachitem.value;
//        }
//        x[i].index = -1;
//        predict_label = predict_probability(model_,x,prob_estimates);
//
//    }
//    free(prob_estimates);
//    free(x);
//}

#endif //MYTOOL_LIBLINEARWRAPPER_H
