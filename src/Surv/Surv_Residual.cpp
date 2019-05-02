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

void martin_resid(const ivec &Censor, const ivec &Y, const ivec &obs, const int &Nb, mat &surv_matrix, vec &MResid){


  for(int i=0; i < Nb; i++){ //For each observation included
    MResid[i] = Censor[obs[i]] - (-log(surv_matrix(i,Y[obs[i]]))); //Find the Martingale residual
  }
}

void dev_resid(const ivec &Censor, const ivec &Y, const ivec &obs, const int &Nb, mat &surv_matrix, vec &DResid){
  vec MResid(Nb);
  MResid.fill(0);
  martin_resid(Censor, Y, obs, Nb, surv_matrix, MResid);

  double minMR = 0;
  for(int i=0; i < Nb; i++){ //For each observation included
    if(std::isinf(-MResid[i])){
      if(minMR == 0){
        for(int j=0; j < Nb; j++){
          if(MResid[j] < minMR &&  std::isinf(-MResid[i])){
            minMR = MResid[j];
          }
        }
      }
      //Rcout << "-infty res replace with " << minMR << std::endl;;
      MResid[i] = minMR;
    }
    DResid[i] = ((MResid[i] > 0) - (MResid[i] < 0))*sqrt(-2*(MResid[i]+Censor[obs[i]]*log(Censor[obs[i]]-MResid[i]))); //Find the Deviance residual
    if(Censor[obs[i]]==0 and MResid[i]==0) DResid[i] = 0;
    //if(DResid[i]!=DResid[i]) Rcout << "Missing Residual: "<< MResid[i]<< " "<< Censor[obs[i]]<< " "<< Y[obs[i]]<<" "<<surv_matrix(i,Y[obs[i]])<<" " << MResid[i]+Censor[obs[i]]*log((Censor[obs[i]]-MResid[i])+.00000001) << std::endl;;
  }

  double maxDR = 0;
  double minDR = 0;
  for(int i=0; i < Nb; i++){
    if(std::isinf(-DResid[i]) || DResid[i] != DResid[i]) {
      if(minDR == 0){
        for(int j=0; j < Nb; j++){
          if(DResid[j] < minDR &&  std::isinf(-DResid[i])){
            minDR = DResid[j];
          }
        }
      }
      DResid[i] = minDR;
    }
    if(std::isinf(DResid[i])) {
      if(maxDR == 0){
        for(int j=0; j < Nb; j++){
          if(DResid[j] > maxDR &&  std::isinf(DResid[i])){
            maxDR = DResid[j];
          }
        }
      }
      //Rcout << "+infty res replace with " << maxDR << std::endl;;
      DResid[i] = maxDR;
    }
  }

}
