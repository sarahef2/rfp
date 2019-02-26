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
# include <math.h>       /* isnan, sqrt */
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
//# include <Rdefines.h>
//# include <R.h>
//# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void Variable_Importance(const std::vector<  colvec > &X,
                         const ivec &Y,
                         const ivec &Censor,
                         const ivec &Ncat,
                         const PARAMETERS* myPara,
                         const vec &subjectweight,
                         std::vector< mat > &tree_matrix,
                         const ivec &subj_id,
                         const int &N,
                         const int &P,
                         imat &ObsTrack,
                         imat &ObsTerminal,
                         std::vector< std::vector< ivec > > &NodeRegi,
                         vec &VarImp,
                         int &use_cores,
                         mat &oob_surv_matrix,
                         vec &oob_residuals){
    int Nfail = myPara->Nfail; //Used in variable importance calculations
    ivec tmp(N);
    tmp.fill(0);
    PredictSurvivalKernel((const std::vector< colvec >) X,
                          Y,
                          Censor,
                          Ncat,
                          subjectweight,
                          (const std::vector< mat >) tree_matrix,
                          (const imat) ObsTrack,
                          ObsTerminal,
                          (const std::vector< std::vector< ivec > >) NodeRegi,
                          oob_surv_matrix,
                          (const PARAMETERS*) myPara,
                          N,
                          (const ivec) subj_id,
                          (const int) -1,
                          (const ivec) tmp,
                          1,
                          true,
                          use_cores);
    dev_resid(Censor, Y, subj_id, N, oob_surv_matrix, oob_residuals);
    double Dev_MSE = 0;
    for(int d=0; d < N; d++){
      Dev_MSE += oob_residuals[d]*oob_residuals[d];
    }
    
    //May want to make this variable?
    int nsim = 1;
    
    for(int j=0; j < P; j++){
      vec imp_sim(nsim);
      imp_sim.fill(0);
      for(int k=0; k < nsim; k++){
        ivec perm_j = subj_id;//oobagObs;
        permute_i(perm_j, N);
        
        mat surv_matrix_perm(N,Nfail+1);
        surv_matrix_perm.fill(0);
        
        PredictSurvivalKernel((const std::vector< colvec >) X,
                              Y,
                              Censor,
                              Ncat,
                              subjectweight,
                              (const std::vector< mat >) tree_matrix, //tree_matrix, //Figure out
                              (const imat) ObsTrack,
                              ObsTerminal,
                              (const std::vector< std::vector< ivec > >) NodeRegi, //NodeRegi, //Figure out
                              surv_matrix_perm,
                              (const PARAMETERS*) myPara,
                              N,
                              (const ivec) subj_id,
                              (const int) j,
                              (const ivec) perm_j,
                              1,
                              true,
                              use_cores);
        vec resids_perm(N);
        resids_perm.fill(0);
        dev_resid(Censor, Y, subj_id, N, surv_matrix_perm, resids_perm);
        double Dev_MSE_perm = 0;
        for(int d=0; d < N; d++){
          Dev_MSE_perm += resids_perm[d]*resids_perm[d];
        }
        //if(Dev_MSE_perm!=Dev_MSE_perm) Rcout << "Dev_MSE_perm is missing var "<< j<<std::endl;;
        imp_sim[k] = Dev_MSE_perm;
        //if(imp_sim[k]!=imp_sim[k]) Rcout << "imp_sim[k] is missing var "<< j<<std::endl;;
      }
      for(int m = 0; m < nsim; m++){
        VarImp[j] += imp_sim[m] / nsim;
        //if(VarImp[j]!=VarImp[j]) Rcout << "A) VarImp[j] is missing var "<< j<< " imp_sim[m] "<< imp_sim[m]<<" nsim "<< nsim<<std::endl;;
      }
      //if(VarImp[j]!=VarImp[j]) Rcout << "B) VarImp[j] is missing var "<< j<<std::endl;;
      
      VarImp[j] = VarImp[j] / Dev_MSE - 1;
    }
}