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

//# include <Rdefines.h>
//# include <R.h>
//# include <Rmath.h>
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

// Likelihood for equal weight version

vec haz(ivec Count_Fail, ivec Count_Censor, int timepoints){
  ivec risk_set(timepoints+1);
  risk_set.fill(0);
  vec haz(timepoints);
  haz.fill(0);
  
  //Rcout <<"Count_Fail: "<< Count_Fail<<std::endl;;
  for(int k = 0; k <= timepoints; k++){
    for(int j=0; j <= k; j++){
      risk_set[j] += Count_Fail[k];
      risk_set[j] += Count_Censor[k];
    }
  }

  for(int k = 1; k <= timepoints; k++){
    if(risk_set[k]==0) risk_set[k] = 1;
    haz[k-1] = ((double) Count_Fail[k])/risk_set[k];
  }
  
  //Rcout <<"hazard: "<< haz<<std::endl;;
  
  return haz;
}

double loglik(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, int timepoints, double w){
  ivec Count_Fail(timepoints+1);
  Count_Fail.fill(0);
  ivec Count_Censor(timepoints+1);
  Count_Censor.fill(0);
  
  for(int i = 0; i <= timepoints; i++){
    Count_Fail[i] = Left_Count_Fail[i] + Right_Count_Fail[i];
    Count_Censor[i] = Left_Count_Censor[i] + Right_Count_Censor[i];
  }
  
  vec lambda0 = haz(Count_Fail, Count_Censor, timepoints);
  vec lambdaLtmp = haz(Left_Count_Fail, Left_Count_Censor, timepoints);
  vec lambdaRtmp = haz(Right_Count_Fail, Right_Count_Censor, timepoints);
  
  vec lambdaL(timepoints+1);
  lambdaL.fill(0);
  vec lambdaR(timepoints+1);
  lambdaR.fill(0);
  
  for(int i = 0; i < timepoints; i++){
    lambdaL[i] = (lambdaLtmp[i]-lambda0[i])*w + lambda0[i];
    lambdaR[i] = (lambdaRtmp[i]-lambda0[i])*w + lambda0[i];
  }
  
  double loglik = 0;
  double temp = 0;
  
  for(int j = 1; j <= timepoints; j++){
    temp = 0;
    for(int k = j; k<=timepoints; k++){
      temp += lambdaL[k-1]*(Left_Count_Fail[k]+Left_Count_Censor[k]) +
        lambdaR[k-1]*(Right_Count_Fail[k]+Right_Count_Censor[k]);
    }
    if(temp>0) loglik -= (Count_Censor[j]+Count_Fail[j])*log(temp);
    if(lambdaL[j-1]>0) loglik += Left_Count_Fail[j]*log(lambdaL[j-1]);
    if(lambdaR[j-1]>0) loglik += Right_Count_Fail[j]*log(lambdaR[j-1]);
  }
  //Rcout <<"Left_Count_Fail[1]: "<< Left_Count_Fail[1]<< " Lhaz[1]: "<<lambdaL[1]<<std::endl;;
  //Rcout <<"Right_Count_Fail[1]: "<< Right_Count_Fail[1]<< " Rhaz[1]: "<<lambdaR[1]<<std::endl;;
  //Rcout <<"loglik: "<< -loglik<<std::endl;;

  return -loglik;
}

vec haz_w(vec Count_Fail, vec Count_Censor, int timepoints){
  vec risk_set(timepoints+1);
  risk_set.fill(0);
  vec haz(timepoints);
  haz.fill(0);
  
  for(int k = 0; k <= timepoints; k++){
    for(int j=0; j <= k; j++){
      risk_set[j] += Count_Fail[k];
      risk_set[j] += Count_Censor[k];
    }
  }
  
  for(int k = 1; k <= timepoints; k++){
    if(risk_set[k]==0) risk_set[k] = 1;
    haz[k-1] = 1.0/risk_set[k];
  }
  return haz;
}

double loglik_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, int timepoints, double w){
  vec Count_Fail(timepoints+1);
  Count_Fail.fill(0);
  vec Count_Censor(timepoints+1);
  Count_Censor.fill(0);
  
  for(int i = 0; i < timepoints; i++){
    Count_Fail[i] = Left_Count_Fail[i] + Right_Count_Fail[i];
    Count_Censor[i] = Left_Count_Censor[i] + Right_Count_Censor[i];
  }
  
  vec lambda0 = haz_w(Count_Fail, Count_Censor, timepoints);
  vec lambdaLtmp = haz_w(Left_Count_Fail, Left_Count_Censor, timepoints);
  vec lambdaRtmp = haz_w(Right_Count_Fail, Right_Count_Censor, timepoints);
  
  vec lambdaL(timepoints+1);
  lambdaL.fill(0);
  vec lambdaR(timepoints+1);
  lambdaR.fill(0);
  
  for(int i = 0; i < timepoints; i++){
    lambdaL[i] = (lambdaL[i]-lambda0[i])*w + lambda0[i];
    lambdaR[i] = (lambdaR[i]-lambda0[i])*w + lambda0[i];
  }
  
  double loglik = 0;
  double temp = 0;
  
  for(int j = 0; j < timepoints; j++){
    temp = 0;
    for(int k = j; k<timepoints; k++){
      temp = temp + lambdaL[k]*(Left_Count_Fail[k]+Left_Count_Censor[k]) + 
        lambdaR[k]*(Right_Count_Fail[k]+Right_Count_Censor[k]);
    }
    if(temp>0) loglik = loglik - (Count_Censor[j]+Count_Fail[j])*log(temp);
    if(Left_Count_Fail[j]*lambdaL[j]>0) loglik= loglik + log(Left_Count_Fail[j]*lambdaL[j]);
    if(Right_Count_Fail[j]*lambdaR[j]>0) loglik = loglik + log(Right_Count_Fail[j]*lambdaR[j]);
  }
  return -loglik;
}

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
      tempscore = dmax(numerator*numerator / denominator, tempscore);
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
      tempscore = dmax(numerator*numerator / denominator, tempscore);
    //else // due to precision loss, this might be 0 already, so just stop here
    //  break;

    LeftWeights -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllWeights -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}
