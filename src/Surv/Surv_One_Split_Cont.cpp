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

void Surv_One_Split_Cont(double* cut,
                         double* score,
                         const ivec &useObs,
                         int node_n,
                         const colvec &x,
                         double &varw, //The weight for this variable
                         const ivec &Y, // y should be called by index i
                         const ivec &Censor, // censor should be called by index i
                         int &timepoints,
                         int &split_gen,
                         int &split_rule,
                         int &nsplit,
                         int &nmin,
                         int &nmin_control,
                         int &nmin_failure)
{
  //Rcout << "Starting new split"<<std::endl;;
  ivec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  ivec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  ivec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);
  ivec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);

  int i, k;

  double temp_score;
  double temp_cut;
  vec lambda0(timepoints);
  lambda0.fill(0);
  double loglik0;
  
  if(split_rule==3){//Finding the base hazard once
    if(varw < 1){
      ivec Count_Fail(timepoints+1);
      Count_Fail.fill(0);
      ivec Count_Censor(timepoints+1);
      Count_Censor.fill(0);
      
      for(i=0; i<node_n; i++){
        if (Censor[i] == 1)
          Count_Fail[Y[i]]++;
        else
          Count_Censor[Y[i]]++;
      }
      lambda0 = haz(Count_Fail, Count_Censor, node_n, timepoints);
      loglik0 = loglik(Count_Fail, Count_Censor, Right_Count_Fail, Right_Count_Censor, node_n, node_n, timepoints, 0, lambda0);
    }
  }
  
  if (split_gen == 1) // random split
  {
    //R_DBP("run random splitting rule");

    for (k = 0; k < nsplit; k++)
    {
      temp_cut = x[useObs[random_in_range(0, node_n)]];
      temp_score = -1;

      double LeftN = 0;

      for (i = 0; i<node_n; i++)
      {
        if (x[useObs[i]] <= temp_cut) // go left
        {
          if (Censor[i] == 1)
            Left_Count_Fail[Y[i]]++;
          else
            Left_Count_Censor[Y[i]]++;

          LeftN++;
        }else{  // go right
          if (Censor[i] == 1)
            Right_Count_Fail[Y[i]]++;
          else
            Right_Count_Censor[Y[i]]++;
        }
      }

      if (LeftN > 0 && LeftN < node_n)
      {
        if (split_rule == 1)
          temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else if (split_rule==2)
          temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else
          temp_score = loglik(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints, varw, lambda0);
      }

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = temp_cut;
      }
    }

    return;
  }

  // rank and best split need to copy x
  vec xtemp(node_n);
  ivec index(node_n);
  ivec censortemp(node_n);

  //R_DBP("Calculate sizes\n");
  for (i = 0; i < node_n; i++)
  {
    xtemp[i] = x[useObs[i]];
    index[i] = i;
  }

  //R_DBP("Sort\n");
  std::sort(index.begin(), index.end(), [&xtemp](size_t i, size_t j) {return (xtemp[i] < xtemp[j]);});

  for(i=0; i<node_n; i++){
    xtemp[i] = x[useObs[index[i]]];
    censortemp[i] = Censor[index[i]];
    //Rcout << "Y[index[" <<i<<"]]="<<Y[index[i]] << std::endl;
  }

  //Rcout << "xtemp = " << xtemp << std::endl;

  int lowindex = 0;
  int highindex = node_n - 1;

  //R_DBP("Control size\n");
  if(nmin_control){
    if(nmin_failure){
      int failcount = 0;
      int fi = 0; //index of observations to count up to the required number of failures
      if(sum(censortemp) < nmin){
        return;//Stop running function- there are too few failures in this node already
      }else{
        failcount = failcount + censortemp[fi];
        while(failcount < nmin){
          fi++;
          failcount = failcount + censortemp[fi];
        }
        lowindex = fi;
        fi = node_n - 1;
        failcount = censortemp[fi];
        while(failcount < nmin){
          fi--;
          failcount += censortemp[fi];
        }
        highindex = fi;
      }
    }else{
      lowindex = nmin - 1;
      highindex = node_n - 1 - lowindex;
    }
  }

  //R_DBP("Ties\n");
  // check for ties
  while((xtemp[lowindex] == xtemp[lowindex+1]) & (lowindex < highindex)) lowindex ++;
  while((xtemp[highindex] == xtemp[highindex+1]) & (lowindex < highindex)) highindex --;
  if ((lowindex == highindex) & (xtemp[lowindex] == xtemp[lowindex+1])) return;

  if (split_gen == 2) // rank split
  {

    int temp_rank;

    for (k = 0; k < nsplit; k++)
    {
      temp_rank = random_in_range(lowindex, highindex+1);
      temp_score = -1;

      // while in the middle of a sequence of ties, either move up or move down
      if (xtemp[temp_rank] == xtemp[temp_rank+1])
      {
        if(R::runif(0, 1) > 0.5)
        {
          do{temp_rank++;} while (xtemp[temp_rank] == xtemp[temp_rank+1]);
        }else{
          do{temp_rank--;} while (xtemp[temp_rank] == xtemp[temp_rank+1]);
        }
      }

      for (i = 0; i<=temp_rank; i++)
      {
        if (Censor[index[i]] == 1){
          Left_Count_Fail[Y[index[i]]]++;
        }else{
          Left_Count_Censor[Y[index[i]]]++;
        }
      }

      for (i = temp_rank+1; i<node_n; i++)
      {
        if (Censor[index[i]] == 1)
          Right_Count_Fail[Y[index[i]]]++;
        else
          Right_Count_Censor[Y[index[i]]]++;
      }

      if (split_rule == 1)
        temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, temp_rank+1, node_n, timepoints);
      else if(split_rule == 2)
        temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, temp_rank+1, node_n, timepoints);
      else
        temp_score = loglik(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, temp_rank+1, node_n, timepoints, varw, lambda0);
      
      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[temp_rank] + xtemp[temp_rank+1])/2;
      }

    }
  }

  //Rcout << "low index is " << lowindex << "; high index is " << highindex << " split rule is " << split_rule << std::endl;


  //R_DBP("Best split\n");
  if (split_gen == 3) // best split
  {

    //R_DBP("run best splitting rule");

    //R_DBP("Place initial\n");
    // place initial
    for (i = 0; i<=lowindex; i++)
    {
      if (Censor[index[i]] == 1){
        Left_Count_Fail[Y[index[i]]]++;
      }else{
        Left_Count_Censor[Y[index[i]]]++;
      }
    }

    for (i = lowindex+1; i<node_n; i++)
    {
      if (Censor[index[i]] == 1)
        Right_Count_Fail[Y[index[i]]]++;
      else
        Right_Count_Censor[Y[index[i]]]++;
    }

    //R_DBP("Move up and calculate score\n");
    // move up and calculate score
    for (i = lowindex; i <= highindex; i++)
    {
      // if ties
      while (xtemp[i] == xtemp[i+1]){
        i++;

        if (Censor[index[i]] == 1)
        {
          Left_Count_Fail[Y[index[i]]]++;
          Right_Count_Fail[Y[index[i]]]--;
        }else{
          Left_Count_Censor[Y[index[i]]]++;
          Right_Count_Censor[Y[index[i]]]--;
        }
      }

      // get score
      if (split_rule == 1)
        temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, i+1, node_n, timepoints);
      else if(split_rule == 2)
        temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, i+1, node_n, timepoints);
      else
        temp_score = (loglik(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, i+1, node_n, timepoints, varw, lambda0));
      //Rcout<<"Cut: "<<xtemp[i]<<" Score: " << temp_score <<" Y[index[i]]: "<<Y[index[i]]<<" index[i]: "<<index[i]<<std::endl;;

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[i] + xtemp[i+1])/2;
      }

      //Rcout << "cut is " << (xtemp[i] + xtemp[i+1])/2 << "; temp score is " << temp_score << "; temp cut is " << *cut << std::endl;

      // get next ready

      if(i+1 <= highindex){
        if (Censor[index[i+1]] == 1)
        {
          Left_Count_Fail[Y[index[i+1]]]++;
          Right_Count_Fail[Y[index[i+1]]]--;
        }else{
          Left_Count_Censor[Y[index[i+1]]]++;
          Right_Count_Censor[Y[index[i+1]]]--;
        }
      }
    }
  }

  //Rcout << " \n" << std::endl;

  return;
}




void Surv_One_Split_Cont_W(double* cut,
                           double* score,
                           const ivec &useObs,
                           int node_n,
                           const vec &x,
                           double &varw, //The weight for this variable
                           const ivec &Y, // y should be called by index i
                           const ivec &Censor, // censor should be called by index i
                           const vec &subjectweight, // subjectweight should be called by index useObs[i]
                           int &timepoints,
                           int &split_gen,
                           int &split_rule,
                           int &nsplit,
                           int &nmin,
                           int &alpha,
                           int &nmin_control,
                           int &nmin_failure)
{

  vec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  vec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  vec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);
  vec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);

  int i, k;

  double temp_score;
  double temp_cut;

  int obs;
  double LeftWeights = 0;
  double RightWeights = 0;
  double tempweight = 0;
  vec lambda0(timepoints);
  lambda0.fill(0);
  
  if(split_rule==3){//Finding the base hazard once
    if(varw < 1){
      vec Count_Fail(timepoints+1);
      Count_Fail.fill(0);
      vec Count_Censor(timepoints+1);
      Count_Censor.fill(0);
      double allweight=0;
      
      for(i=0; i<node_n; i++){
        obs = useObs[i];
        tempweight = subjectweight[obs];
        allweight += tempweight;
        
        if (Censor[i] == 1)
          Count_Fail[Y[i]] += tempweight;
        else
          Count_Censor[Y[i]] += tempweight;
      }
      lambda0 = haz_w(Count_Fail, Count_Censor, allweight, timepoints);
    }
  }
  
    if (split_gen == 1) // random split
  {

    for (k = 0; k < nsplit; k++)
    {
      temp_cut = x[useObs[random_in_range(0, node_n)]];
      temp_score = -1;

      LeftWeights = 0;
      RightWeights = 0;

      for (i = 0; i<node_n; i++)
      {
        obs = useObs[i];
        tempweight = subjectweight[obs];

        if (x[obs] <= temp_cut) // go left
        {
          if (Censor[i] == 1)
            Left_Count_Fail[Y[i]] += tempweight;
          else
            Left_Count_Censor[Y[i]] += tempweight;

          LeftWeights += subjectweight[obs];
        }else{  // go right
          if (Censor[i] == 1)
            Right_Count_Fail[Y[i]] += tempweight;
          else
            Right_Count_Censor[Y[i]] += tempweight;

          RightWeights += tempweight;
        }
      }

      if (LeftWeights > 0 && RightWeights > 0)
      {
        if (split_rule == 1)
          temp_score = logrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else if(split_rule == 2)
          temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else
          temp_score = loglik_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints, varw, lambda0);
      }

      if (temp_score > *score)//replaced score* with score
      {
        *score = temp_score;//replaced *score with score
        *cut = temp_cut;//replaced *cut with cut
      }
    }

    return;
  }

  // rank and best split need to copy x

  vec xtemp(node_n);
  ivec index(node_n);

  for (i = 0; i < node_n; i++)
  {
    xtemp[i] = x[useObs[i]];
    index[i] = i;
  }

  std::sort(index.begin(), index.end(), [&xtemp](size_t i, size_t j) {return (xtemp[i] < xtemp[j]);});

  ivec censortemp(node_n);
  
  for (i = 0; i < node_n; i++)
  {
    xtemp[i] = x[useObs[index[i]]];
    censortemp[i] = Censor[index[i]];
  }

  int lowindex = 0;
  int highindex = node_n - 1;
  
  if(nmin_control){
    if(nmin_failure){
      //Needs to be updated to use weights!!!
      int failcount = 0;
      int fi = 0; //index of observations to count up to the required number of failures
      if(sum(censortemp) < nmin){
        return;//Stop running function- there are too few failures in this node already
      }else{
        failcount = failcount + censortemp[fi];
        while(failcount < nmin){
          fi++;
          failcount = failcount + censortemp[fi];
        }
        lowindex = fi;
        fi = node_n - 1;
        failcount = censortemp[fi];
        while(failcount < nmin){
          fi--;
          failcount += censortemp[fi];
        }
        highindex = fi;
      }
    }else{
      int lowindex = max(nmin, (int) alpha*node_n) - 1;
      int highindex = node_n - 1 - lowindex;
    }
  }

  // check for ties
  while((xtemp[lowindex] == xtemp[lowindex+1]) & (lowindex < highindex)) lowindex ++;
  while((xtemp[highindex] == xtemp[highindex+1]) & (lowindex < highindex)) highindex --;
  if ((lowindex == highindex) & (xtemp[lowindex] == xtemp[lowindex+1])) return;

  if (split_gen == 2) // rank split
  {


    int temp_rank;

    for (k = 0; k < nsplit; k++)
    {
      temp_rank = random_in_range(lowindex, highindex+1);
      temp_score = -1;

      // while in the middle of a sequence of ties, either move up or move down
      if (xtemp[temp_rank] == xtemp[temp_rank+1])
      {
        if(R::runif(0, 1) > 0.5)
        {
          do{temp_rank++;} while (xtemp[temp_rank] == xtemp[temp_rank+1]);
        }else{
          do{temp_rank--;} while (xtemp[temp_rank] == xtemp[temp_rank+1]);
        }
      }

      LeftWeights = 0;
      RightWeights = 0;

      for (i = 0; i<=temp_rank; i++)
      {
        obs = index[i];
        if (Censor[obs] == 1)
          Left_Count_Fail[Y[obs]] += subjectweight[useObs[obs]];
        else
          Left_Count_Censor[Y[obs]] += subjectweight[useObs[obs]];

        LeftWeights += subjectweight[useObs[obs]];
      }

      for (i = temp_rank+1; i<node_n; i++)
      {
        obs = index[i];
        if (Censor[obs] == 1)
          Right_Count_Fail[Y[obs]] += subjectweight[useObs[obs]];
        else
          Right_Count_Censor[Y[obs]] += subjectweight[useObs[obs]];

        RightWeights += subjectweight[useObs[obs]];
      }

      if (split_rule == 1)
        temp_score = logrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);
      else if(split_rule == 2)
        temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);
      else
        temp_score = loglik_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints, varw, lambda0);
      
      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[temp_rank] + xtemp[temp_rank+1])/2;
      }
    }
  }

  if (split_gen == 3) // best split
  {


    LeftWeights = 0;
    RightWeights = 0;

    // place initial
    for (i = 0; i<=lowindex; i++)
    {

      obs = index[i];
      tempweight = subjectweight[useObs[obs]];

      if (Censor[obs] == 1){
        Left_Count_Fail[Y[obs]] += tempweight;
      }else
        Left_Count_Censor[Y[obs]] += tempweight;

      LeftWeights += tempweight;
    }

    for (i = lowindex+1; i<node_n; i++)
    {
      obs = index[i];
      tempweight = subjectweight[useObs[obs]];

      if (Censor[obs] == 1)
        Right_Count_Fail[Y[obs]] += tempweight;
      else
        Right_Count_Censor[Y[obs]] += tempweight;

      RightWeights += tempweight;
    }

    // move up and calculate score
    for (i = lowindex; i <= highindex; i++)
    {
      // if ties
      while (xtemp[i] == xtemp[i+1]){
        i++;

        obs = index[i];
        tempweight = subjectweight[useObs[obs]];

        if (Censor[obs] == 1)
        {
          Left_Count_Fail[Y[obs]] += tempweight;
          Right_Count_Fail[Y[obs]] -= tempweight;
        }else{
          Left_Count_Censor[Y[obs]] += tempweight;
          Right_Count_Censor[Y[obs]] -= tempweight;
        }

        LeftWeights += tempweight;
        RightWeights -= tempweight;
      }

      // get score
      if (split_rule == 1)
        temp_score = logrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);
      else if(split_rule == 2)
        temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);
      else
        temp_score = loglik_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints, varw, lambda0);
      
      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[i] + xtemp[i+1])/2;;
      }

      // get next ready
      obs = index[i+1];
      tempweight = subjectweight[useObs[obs]];

      if (Censor[obs] == 1)
      {
        Left_Count_Fail[Y[obs]] += tempweight;
        Right_Count_Fail[Y[obs]] -= tempweight;
      }else{
        Left_Count_Censor[Y[obs]] += tempweight;
        Right_Count_Censor[Y[obs]] -= tempweight;
      }

      LeftWeights += tempweight;
      RightWeights -= tempweight;

    }
  }

  return;

}
