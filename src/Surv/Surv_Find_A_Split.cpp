//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Survival
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// my header file
# include "..//survForest.h"
# include "..//Utility//utility.h"

void Surv_Find_A_Split(int* splitVar,
                       double* splitVal,
                       const std::vector<  colvec > X,
                       const ivec Y,
                       const ivec Censor,
                       const ivec Ncat,
                       const vec Interval,
                       const PARAMETERS* myPara,
                       const vec subjectweight,
                       ivec useObs,
                       const int node_n,
                       vec variableweight,
                       ivec variableindex,
                       const int P,
                       int &counter)
{

  int N = myPara->N;
  int nmin = myPara->nmin;
  int nmin_control = myPara->nmin_control;
  int nmin_failure = myPara->nmin_failure;
  int mtry = myPara->mtry;
  int use_sub_weight = myPara->use_sub_weight;
  int use_var_weight = myPara->use_var_weight;
  int nsplit = myPara->nsplit;
  int split_gen = myPara->split_gen;
  int split_rule = myPara->split_rule;
  int alpha = myPara->alpha;

  int j;

  int mincount = max(nmin, (int) alpha*node_n);
  //double minweight = (double) max(nmin, (int) alpha*node_n) / N;

  // collapse Y into contiguous integers

  int timepoints = 0;
  ivec Y_collapse(node_n);
  ivec Censor_collapse(node_n);

  collapse(Y, Censor, Y_collapse, Censor_collapse, useObs, node_n, timepoints);
  
  int temp_var;
  double temp_val;
  double temp_score;
  double best_score = -1;

  // calculate node information

  int mtry_remaining = max(1, min(mtry, P));

  ivec var_used(P);
  var_used.fill(0);

  for (j = 0; j < mtry_remaining; j++)
  {
    temp_val = 0;
    temp_score = -1;

    temp_var = sample_rotate(variableindex, variableweight, j, P);

    counter++;

    double temp_vw = 1.0;
    if(use_var_weight) temp_vw=variableweight[j]; 

    if (Ncat[temp_var] > 1)
    {
      if (use_sub_weight)
      {
        Surv_One_Split_Cat_W(&temp_val, &temp_score, (const ivec) useObs, node_n, X[temp_var], temp_vw, Y_collapse, Censor_collapse, subjectweight, Ncat[temp_var],
                             timepoints, split_gen, split_rule, nsplit, nmin, alpha);
      }else{

        Surv_One_Split_Cat(&temp_val, &temp_score, (const ivec) useObs, node_n, X[temp_var], temp_vw, Y_collapse, Censor_collapse, Ncat[temp_var],
                           timepoints, split_gen, split_rule, nsplit, mincount);

      }

    }else{

      if (use_sub_weight)
      {
        Surv_One_Split_Cont_W(&temp_val, &temp_score, (const ivec) useObs, node_n, X[temp_var], temp_vw, Y_collapse, Censor_collapse, subjectweight,
                              timepoints, split_gen, split_rule, nsplit, nmin, alpha, nmin_control, nmin_failure);

      }else{
        Surv_One_Split_Cont(&temp_val, &temp_score, (const ivec) useObs, node_n, X[temp_var], temp_vw, Y_collapse, Censor_collapse,
                            timepoints, split_gen, split_rule, nsplit, mincount, nmin_control, nmin_failure);
      }
    }

    if (use_var_weight and split_rule<3)
      temp_score = temp_score*temp_vw;

    // update the score
    if (temp_score > 0 && temp_score > best_score)
    {
      best_score = temp_score;
      *splitVar = temp_var;
      *splitVal = temp_val;
    }
  }

  return;
}


// collapse Y into contiguous integers, Y will always be in an increasing order

void collapse(const ivec Y, const ivec Censor, ivec &Y_collapse, ivec &Censor_collapse, const ivec useObs, int node_n, int &timepoints)
{
  int i=0;

  while(Y[useObs[i]] == 0){
    Y_collapse[i] = 0;
    Censor_collapse[i] = Censor[useObs[i]];
    i++;
  }

  int hold_y = Y[useObs[i]];
  int counter = 1;
  int j;
  int this_y;
  int have_fail;

  for (; i < node_n; i++)
  {

    if (Y[useObs[i]] > hold_y) // a new Y value should we create failure time?
    {
      have_fail = 0;
      this_y = Y[useObs[i]];
      for (j = i; j<node_n; j++)
      {
        if (Y[useObs[j]] > this_y)
          break;

        if (Censor[useObs[j]] == 1)
        {
          have_fail = 1;
          break;
        }
      }

      if (have_fail == 1)
      {
        counter++;
        hold_y = Y[useObs[i]];
      }
    }

    Y_collapse[i] = counter;
    Censor_collapse[i] = Censor[useObs[i]];
  }

  timepoints = counter;
}
