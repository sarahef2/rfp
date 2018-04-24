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
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
using namespace Rcpp;
# include <Rmath.h>

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void Surv_One_Split_Cont(double* cut,
                         double* score,
                         const ivec useObs,
                         int node_n,
                         //const double* x, // x should be called by index useObs[i]
                         const colvec x,
                         const ivec Y, // y should be called by index i
                         const ivec Censor, // censor should be called by index i
                         int timepoints,
                         int split_gen,
                         int split_rule,
                         int nsplit,
                         int mincount)
{
  //int *Left_Count_Fail = (int *) malloc((timepoints+1)* sizeof(int));
  //int *Left_Count_Fail = new int[timepoints+1];
  ivec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  //int *Left_Count_Censor = (int *) malloc((timepoints+1)* sizeof(int));
  ivec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  //int *Right_Count_Fail = (int *) malloc((timepoints+1)* sizeof(int));
  ivec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);
  //int *Right_Count_Censor = (int *) malloc((timepoints+1)* sizeof(int));
  ivec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);

  int i, k;

  double temp_score;
  double temp_cut;

  if (split_gen == 1) // random split
  {

    //R_DBP("run random splitting rule");

    for (k = 0; k < nsplit; k++)
    {
      temp_cut = x[useObs[random_in_range(0, node_n)]];
      temp_score = -1;

      double LeftN = 0;

      //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(int));
      //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(int));
      //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(int));
      //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(int));

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
        else
          temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
      }

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = temp_cut;
      }
    }

    //free(Left_Count_Fail);
    //delete[] Left_Count_Fail;
    //free(Left_Count_Censor);
    //delete[] Left_Count_Censor;
    //free(Right_Count_Fail);
    //delete[] Right_Count_Fail;
    //free(Right_Count_Censor);
    //delete[] Right_Count_Censor;
    return;
  }

  // rank and best split need to copy x

  //double* xtemp = (double *) malloc(node_n * sizeof(double));
  vec xtemp(node_n);
  //int* index = (int *) malloc(node_n * sizeof(int));
  uvec index(node_n);

  for (i = 0; i < node_n; i++)
  {
    xtemp[i] = x[useObs[i]];
    index[i] = i;
  }

  //qSort_dindex(xtemp, 0, node_n-1, index);
  index = sort_index(xtemp);
  xtemp = sort(xtemp);

  int lowindex = mincount - 1;
  int highindex = node_n - 1 - lowindex;

  // check for ties
  while((xtemp[lowindex] == xtemp[lowindex+1]) & (lowindex < highindex)) lowindex ++;
  while((xtemp[highindex] == xtemp[highindex+1]) & (lowindex < highindex)) highindex --;
  if ((lowindex == highindex) & (xtemp[lowindex] == xtemp[lowindex+1])) return;

  if (split_gen == 2) // rank split
  {

    //R_DBP("run rank splitting rule");
    int temp_rank;

    for (k = 0; k < nsplit; k++)
    {
      temp_rank = random_in_range(lowindex, highindex+1);
      temp_score = -1;

      //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(int));
      //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(int));
      //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(int));
      //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(int));

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

      //R_DBP(" Rank is %i out of %i \n", temp_rank, node_n);
      //for (int j = 0; j < timepoints + 1; j ++)
      //  R_DBP("%i, %i, %i, %i \n", Left_Count_Fail[j], Left_Count_Censor[j], Right_Count_Fail[j], Right_Count_Censor[j]);

      if (split_rule == 1)
        temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, temp_rank+1, node_n, timepoints);
      else
        temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, temp_rank+1, node_n, timepoints);

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[temp_rank] + xtemp[temp_rank+1])/2;
      }
    }
  }



  if (split_gen == 3) // best split
  {

    //R_DBP("run best splitting rule");

    //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(int));
    //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(int));
    //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(int));
    //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(int));

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
      else
        temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, i+1, node_n, timepoints);

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[i] + xtemp[i+1])/2;;
      }

      // get next ready

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

  //free(Left_Count_Fail);
  //delete[] Left_Count_Fail;
  //free(Left_Count_Censor);
  //delete[] Left_Count_Censor;
  //free(Right_Count_Fail);
  //delete[] Right_Count_Fail;
  //free(Right_Count_Censor);
  //delete[] Right_Count_Censor;

  //free(xtemp);
  //delete[] xtemp;
  //free(index);
  //delete[] index;
  return;
}




void Surv_One_Split_Cont_W(double* cut,
                           double* score,
                           const ivec useObs,
                           int node_n,
                           //const double* x, // x should be called by index useObs[i]
                           const vec x,
                           const ivec Y, // y should be called by index i
                           const ivec Censor, // censor should be called by index i
                           const vec subjectweight, // subjectweight should be called by index useObs[i]
                           int timepoints,
                           int split_gen,
                           int split_rule,
                           int nsplit,
                           int nmin,
                           int alpha)
{

  //double *Left_Count_Fail = (double *) malloc((timepoints+1)* sizeof(double));
  vec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  //double *Left_Count_Censor = (double *) malloc((timepoints+1)* sizeof(double));
  vec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  //double *Right_Count_Fail = (double *) malloc((timepoints+1)* sizeof(double));
  vec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);
  //double *Right_Count_Censor = (double *) malloc((timepoints+1)* sizeof(double));
  vec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);
  
  int i, k;

  double temp_score;
  double temp_cut;

  int obs;
  double LeftWeights = 0;
  double RightWeights = 0;
  double tempweight = 0;

    if (split_gen == 1) // random split
  {

    for (k = 0; k < nsplit; k++)
    {
      temp_cut = x[useObs[random_in_range(0, node_n)]];
      temp_score = -1;

      LeftWeights = 0;
      RightWeights = 0;

      //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(double));
      //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(double));
      //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(double));
      //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(double));

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
        else
          temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
      }

      if (temp_score > *score)//replaced score* with score
      {
        *score = temp_score;//replaced *score with score
        *cut = temp_cut;//replaced *cut with cut
      }
    }

    //free(Left_Count_Fail);
    //delete[] Left_Count_Fail;
    //free(Left_Count_Censor);
    //delete[] Left_Count_Censor;
    //free(Right_Count_Fail);
    //delete[] Right_Count_Fail;
    //free(Right_Count_Censor);
    //delete[] Right_Count_Censor;
    return;
  }

  // rank and best split need to copy x

  //double* xtemp = (double *) malloc(node_n * sizeof(double));
  vec xtemp(node_n);
  //int* index = (int *) malloc(node_n * sizeof(int));
  uvec index(node_n);

  for (i = 0; i < node_n; i++)
  {
    xtemp[i] = x[useObs[i]];
    index[i] = i;
  }

  //qSort_dindex(xtemp, 0, node_n-1, index);
  index = sort_index(xtemp);
  xtemp = sort(xtemp);

  int lowindex = imax(nmin, (int) alpha*node_n) - 1;
  int highindex = node_n - 1 - lowindex;

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

      //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(double));
      //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(double));
      //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(double));
      //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(double));

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
      else
        temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);

      if (temp_score > *score)
      {
        *score = temp_score;
        *cut = (xtemp[temp_rank] + xtemp[temp_rank+1])/2;
      }
    }
  }

  if (split_gen == 3) // best split
  {

    
    //memset(Left_Count_Fail, 0, (timepoints+1)*sizeof(double));
    //memset(Left_Count_Censor, 0, (timepoints+1)*sizeof(double));
    //memset(Right_Count_Fail, 0, (timepoints+1)*sizeof(double));
    //memset(Right_Count_Censor, 0, (timepoints+1)*sizeof(double));

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
      else
        temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights + RightWeights, timepoints);

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

  //free(Left_Count_Fail);
  //delete[] Left_Count_Fail;
  //free(Left_Count_Censor);
  //delete[] Left_Count_Censor;
  //free(Right_Count_Fail);
  //delete[] Right_Count_Fail;
  //free(Right_Count_Censor);
  //delete[] Right_Count_Censor;
  
  //free(xtemp);
  //delete[] xtemp;
  //free(index);
  //delete[] index;


  return;

}