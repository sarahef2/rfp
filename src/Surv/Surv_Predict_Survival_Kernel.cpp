//  **********************************************************************
//
//    Survival Forests (survForest)
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 3
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public
//    License along with this program; if not, write to the Free
//    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301, USA.
//
//  **********************************************************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
//# include <Rdefines.h>
//# include <R.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void PredictSurvivalKernel(const std::vector< colvec > &X,
                           const ivec &Y,
                           const ivec &Censor,
                           const ivec &Ncat,
                           const vec &subjectweight,
                           const std::vector< mat > &tree_matrix,
                           const imat &ObsTrack,
                           imat &ObsTerminal,
                           const std::vector< std::vector< ivec > > &NodeRegi,
                           mat &surv_matrix,
                           const PARAMETERS* myPara,
                           int testN, //The number of observations to predict
                           const ivec &use_obs, //The index of the obserations to predict
                           const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                           const ivec &perm_j, //The permuted variable (for variable importance)
                           int oob_only, //Should prediction be run on only oob obs?
                           bool InTrainSet,
                           int use_cores){
  

  int N = myPara->N;
  int Nfail = myPara->Nfail;
  int ntrees = myPara->ntrees;
  if (tree_matrix.size() == 1) ntrees=1;
  int verbose = myPara->verbose;
  int use_sub_weight = myPara->use_sub_weight;

  // parallel computing... set cores
  
  use_cores = imax(1, use_cores);

  int haveCores = omp_get_max_threads();
  
  if(use_cores > haveCores)
  {
    if (verbose) Rprintf("Do not have %i cores, use maximum %i cores. \n", use_cores, haveCores);
    use_cores = haveCores;
  }
  
  mat remove_matrix(testN,Nfail+1);
  remove_matrix.fill(0);

 #pragma omp parallel for schedule(guided) num_threads(use_cores)
  for (int i = 0; i < testN; i++)
  {
    vec weights(N);
    weights.fill(0);
    int j, nt;

    for(nt = 0; nt < ntrees; nt++){
      bool Checkj;
      int perm_ind_tmp = perm_ind;
      if(oob_only){
        CheckVar(tree_matrix[nt], perm_ind, Checkj);
        if(not Checkj){
          perm_ind_tmp = -1;
        }
      }
      if(ObsTrack(use_obs[i],nt) < 1 or !oob_only){ //If the observation is not used to build the tree OR if we want to use all obs, not just oob obs (so oob_only=FALSE)
        if(use_sub_weight){
          Get_Kernel_Weights_w(use_obs[i], X, Ncat, tree_matrix[nt], ObsTerminal,
                               NodeRegi[nt], subjectweight, weights, N, InTrainSet,  perm_ind_tmp, perm_j, i, nt);
        }else{
          Get_Kernel_Weights(use_obs[i], X, Ncat, tree_matrix[nt], ObsTerminal,
                             NodeRegi[nt], weights, N, InTrainSet, perm_ind_tmp, perm_j, i, nt);
        }
      }
    }
    
    double weights_sum = 0;

    for (j = 0; j < N; j++)
    {
      remove_matrix(i,Y[j]) += weights[j];//[i][Y[j]]

      if (Censor[j] == 1)
        surv_matrix(i,Y[j]) += weights[j];

      weights_sum += weights[j];
    }

    surv_matrix(i,0) = 1;
    weights_sum -= remove_matrix(i,0);
    
    // KM survival function
    for (j = 1; j <= Nfail; j++)
    {
      if(weights_sum > 0){
        surv_matrix(i,j) = surv_matrix(i,j-1) * (1 - surv_matrix(i,j)/weights_sum);
      }else{
        //Or stop
        surv_matrix(i,j) = 0;
      }
      weights_sum -= remove_matrix(i,j);
    }
  }

return;

}
  