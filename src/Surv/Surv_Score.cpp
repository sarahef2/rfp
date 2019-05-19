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


// log rank and sup log rank for equal weight version

double logrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints)
{
  double numerator = 0;
  double denominator = 0;
  double tempscore = -1;

  // calculate the logrank for this split
  LeftN -= Left_Count_Censor[0];
  AllN -= Left_Count_Censor[0] - Right_Count_Censor[0];

  for (int j = 1; j <= timepoints && AllN > 1; j++)
  {
    numerator += LeftN*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllN - Left_Count_Fail[j];
    denominator += LeftN*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllN*(1 - LeftN / AllN)*(AllN - Left_Count_Fail[j] - Right_Count_Fail[j]) / (AllN - 1);

    if (denominator > 0)
      tempscore = numerator*numerator / denominator;
    //else{
    //  Rcout<<"LeftN "<<LeftN<<" Left_Count_Fail[j] "<<Left_Count_Fail[j]<<" Right_Count_Fail[j] "<<Right_Count_Fail[j]<<" AllN "<<AllN<<std::endl;;
    //  break;
    //}

    LeftN -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllN -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];

  }

  return tempscore;
}



double suplogrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints)
{
  double numerator = 0;
  double denominator = 0;
  double tempscore = -1;

  // calculate the logrank for this split
  LeftN -= Left_Count_Censor[0];
  AllN -= Left_Count_Censor[0] - Right_Count_Censor[0];

  for (int j = 1; j <= timepoints && AllN > 1; j++)
  {
    numerator += LeftN*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllN - Left_Count_Fail[j];
    denominator += LeftN*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllN*(1 - LeftN / AllN)*(AllN - Left_Count_Fail[j] - Right_Count_Fail[j]) / (AllN - 1);

    if (denominator > 0)
      tempscore = max(numerator*numerator / denominator, tempscore);
    //else
    //  break;

    LeftN -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllN -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}


// log rank and sup log rank for subject weighted version

double logrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftWeights, double AllWeights, int timepoints)
{
  double numerator = 0;
  double denominator = 0;
  double tempscore = -1;

  // calculate the logrank for this split
  LeftWeights -= Left_Count_Censor[0];
  AllWeights -= Left_Count_Censor[0] - Right_Count_Censor[0];

  for (int j = 1; j <= timepoints; j++)
  {
    numerator += LeftWeights*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllWeights - Left_Count_Fail[j];
    denominator += LeftWeights*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllWeights*(1 - LeftWeights / AllWeights)*(AllWeights - Left_Count_Fail[j] - Right_Count_Fail[j]) / AllWeights;

    if (denominator > WeightTH)
      tempscore = numerator*numerator / denominator;
    //else // due to precision loss, this might be 0 already, so just stop here
    //  break;

    LeftWeights -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllWeights -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}



double suplogrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftWeights, double AllWeights, int timepoints)
{
  double numerator = 0;
  double denominator = 0;
  double tempscore = -1;

  // calculate the logrank for this split
  LeftWeights -= Left_Count_Censor[0];
  AllWeights -= Left_Count_Censor[0] - Right_Count_Censor[0];

  for (int j = 1; j <= timepoints; j++)
  {
    numerator += LeftWeights*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllWeights - Left_Count_Fail[j];
    denominator += LeftWeights*(Left_Count_Fail[j] + Right_Count_Fail[j]) / AllWeights*(1 - LeftWeights / AllWeights)*(AllWeights - Left_Count_Fail[j] - Right_Count_Fail[j]) / AllWeights;

    if (denominator > WeightTH)
      tempscore = max(numerator*numerator / denominator, tempscore);
    //else // due to precision loss, this might be 0 already, so just stop here
    //  break;

    LeftWeights -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllWeights -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}

// Likelihood for equal weight version

vec haz(ivec Count_Fail, ivec Count_Censor, double N, int timepoints){
  int at_risk=N;
  vec haz(timepoints);
  haz.fill(0);
  
  for(int k = 0; (k < timepoints)&& at_risk>0; k++){
    haz[k] = ((double) Count_Fail[k+1])/at_risk;
    at_risk -= Count_Fail[k+1];
    at_risk -= Count_Censor[k+1];
  }
  
  return haz;
}

double loglik(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints, double w, vec &lambda0){
  if(timepoints==1) return -1;
  //Rewrite so that if w=1, skip one of the initializations.  
  vec lambdaLtmp = haz(Left_Count_Fail, Left_Count_Censor, LeftN, timepoints);
  vec lambdaRtmp = haz(Right_Count_Fail, Right_Count_Censor, AllN-LeftN, timepoints);
  
  vec lambdaL(timepoints+1);
  lambdaL.fill(0);
  vec lambdaR(timepoints+1);
  lambdaR.fill(0);
  
  for(int i = 0; i < timepoints; i++){
    lambdaL[i] = (lambdaLtmp[i]-lambda0[i])*w + lambda0[i];//Change w to 0.0001 for everyone, calc change in likelihood, and then divide by the w.  Each variable has its own gradient, and we penalize the gradient
    lambdaR[i] = (lambdaRtmp[i]-lambda0[i])*w + lambda0[i];
  }
  
  double loglik = 0;
  double tempL = 0;
  double tempR = 0;
  
  for(int j = 1; j <= timepoints; j++){
    tempL += lambdaL[j-1];
    tempR += lambdaR[j-1];
    loglik -= (Left_Count_Fail[j]+Left_Count_Censor[j])*tempL;
    loglik -= (Right_Count_Fail[j]+Right_Count_Censor[j])*tempR;
    if(lambdaL[j-1]>0) loglik += Left_Count_Fail[j]*log(lambdaL[j-1]);
    if(lambdaR[j-1]>0) loglik += Right_Count_Fail[j]*log(lambdaR[j-1]);
  }
  
  return -1/loglik;
}

// Likelihood for subject weighted version

vec haz_w(vec Count_Fail, vec Count_Censor, double Nw, int timepoints){
  double at_risk=Nw;
  vec haz(timepoints);
  haz.fill(0);
  
  for(int k = 0; (k < timepoints)&& at_risk>0; k++){
    haz[k] = ((double) Count_Fail[k+1])/at_risk;
    at_risk -= Count_Fail[k+1];
    at_risk -= Count_Censor[k+1];
  }
  
  return haz;
}

double loglik_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftN, double AllN, int timepoints, double w, vec &lambda0){
  vec lambdaLtmp = haz_w(Left_Count_Fail, Left_Count_Censor, LeftN, timepoints);
  vec lambdaRtmp = haz_w(Right_Count_Fail, Right_Count_Censor, AllN-LeftN, timepoints);
  
  vec lambdaL(timepoints+1);
  lambdaL.fill(0);
  vec lambdaR(timepoints+1);
  lambdaR.fill(0);
  
  for(int i = 0; i < timepoints; i++){
    lambdaL[i] = (lambdaLtmp[i]-lambda0[i])*w + lambda0[i];
    lambdaR[i] = (lambdaRtmp[i]-lambda0[i])*w + lambda0[i];
  }
  
  double loglik = 0;
  double tempL = 0;
  double tempR = 0;
  
  for(int j = 1; j <= timepoints; j++){
    tempL += lambdaL[j-1];
    tempR += lambdaR[j-1];
    if(tempL>0) loglik -= (Left_Count_Fail[j]+Left_Count_Censor[j])*tempL;
    if(tempR>0) loglik -= (Right_Count_Fail[j]+Right_Count_Censor[j])*tempR;
    if(lambdaL[j-1]>0) loglik += Left_Count_Fail[j]*log(lambdaL[j-1]);
    if(lambdaR[j-1]>0) loglik += Right_Count_Fail[j]*log(lambdaR[j-1]);
  }
  
  return -1/loglik;
}