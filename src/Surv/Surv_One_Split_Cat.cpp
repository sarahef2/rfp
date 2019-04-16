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
# include <cstring>
# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void Surv_One_Split_Cat(double* cut,
                        double* score,
                        const ivec useObs,
                        int node_n,
                        const vec x,
                        const ivec Y, // y should be called by index i
                        const ivec Censor, // censor should be called by index i
                        int ncat,
                        int timepoints,
                        int split_gen,
                        int split_rule,
                        int nsplit,
                        int mincount)
{

  int i, k, j, temp_cat;
  double LeftN;

  double temp_score = -1;

  // summerize all categories in this node
  SURVCAT* Cat_Count = new SURVCAT[ncat];
    
  //Moved since the function couldn't find them
  ivec goright(ncat);
  goright.fill(0);
  ivec tempRight(ncat);
  tempRight.fill(0);
  ivec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  ivec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  ivec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);
  ivec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);
    

  for (i=0; i< ncat; i++)
  {
    Cat_Count[i].cat = i;
    Cat_Count[i].f = 0;
    Cat_Count[i].c = 0;
    Cat_Count[i].flist = ivec(timepoints+1);
    Cat_Count[i].clist = ivec(timepoints+1);
  }

  // record each category
  for (i=0; i<node_n; i++)
  {
    temp_cat = (int) x[useObs[i]] -1;

    if (Censor[i] == 0)
    {
      Cat_Count[temp_cat].clist[Y[i]]++;
      Cat_Count[temp_cat].c++;
    }else{
      Cat_Count[temp_cat].flist[Y[i]]++;
      Cat_Count[temp_cat].f++;
    }
  }

  // put nonzero categories to the front

  int true_ncat = ncat;
  for (i=0; i < true_ncat; i++)
  {
    if (Cat_Count[i].f + Cat_Count[i].c <= 0)
    {
      swap_SURVCAT(&Cat_Count[i], &Cat_Count[true_ncat-1]);
      true_ncat--;
      i--;
    }
  }

  // put categories with nonzero failures to the front

  int true_ncat_f = true_ncat;

  for (i =0; i < true_ncat_f; i++)
  {
    if (Cat_Count[i].f <= 0)
    {
      swap_SURVCAT(&Cat_Count[i], &Cat_Count[true_ncat_f-1]);
      true_ncat_f--;
      i--;
    }
  }

  if (true_ncat_f <= 1)
    goto NothingToFind;

  // if there are too many categories, save some computation by switching to random
  if (split_gen == 3 && true_ncat > 6)
  {
    split_gen = 1;
    nsplit = 64;
  }
  
  // if there are too many nsplit, switch to best
  if (split_gen < 3 && nsplit > pow(2, true_ncat-1) -1 )
  {
    split_gen = 3;
  }
  
  // random split
  if (split_gen == 1 || split_gen == 2)
  {
    
    for (k =0; k<nsplit; k++)
    {
      temp_score = -1;
      LeftN = 0;
      
      // generate a random split that puts some cat into right
      for (i=0; i<true_ncat_f; i++)
        tempRight[i] = 1;
      
      //memset(tempRight, 0, random_in_range(1, true_ncat_f)*sizeof(int));
      permute_i(tempRight, true_ncat_f);
      
      // We should be checking if there will be enough observations for each side,
      // but I will skip it. Might need to fix it later
      
      // for the nonfailure categories, randomly assign
      for (i=true_ncat_f; i<ncat; i++)
        tempRight[i] = (int) unif_rand()>0.5;
      
      for (i = 0; i < true_ncat_f; i++)
      {
        if (tempRight[i] == 0) // go left
        {
          for (j = 1; j < (timepoints+1); j++)
          {
            Left_Count_Censor[j] += Cat_Count[i].clist[j];
            Left_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          LeftN += Cat_Count[i].f + Cat_Count[i].c;
        }else{ // go right
          for (j = 1; j < (timepoints+1); j++)
          {
            Right_Count_Censor[j] += Cat_Count[i].clist[j];
            Right_Count_Fail[j] += Cat_Count[i].flist[j];
          }
        }
      }
      
      if (LeftN > 0 && LeftN < node_n)
      {
        if (split_rule == 1)
          temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else if (split_rule == 2)
          temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else 
          temp_score = loglik(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, timepoints, varw);
      }
      
      if (temp_score > *score and split_rule<3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
          
          *score = temp_score;
          *cut = pack(ncat, goright);
      }
      if (temp_score < *score and split_rule==3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
        
        *score = temp_score;
        *cut = pack(ncat, goright);
      }
    }
  }
  
  // best split
  
  if (split_gen == 3)
  {

    nsplit = pow(2, true_ncat-1) -1;
    
    for (k = 0; k < nsplit; k++)
    {
      
      LeftN = 0;
      temp_score = -1;
      
      // fit the next possible split rule using binary vectors
      tempRight[0]++;

      for (i = 0; i < true_ncat - 1; i++)
      {
        if (tempRight[i] == 1)
        {
          tempRight[i] = 0;
          tempRight[i+1] ++;
        }
      }
      
      for (i = 0; i < true_ncat; i++)
      {
        if (tempRight[i] == 0)
        {
          for (j = 1; j < (timepoints+1); j++)
          {
            Left_Count_Censor[j] += Cat_Count[i].clist[j];
            Left_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          LeftN += Cat_Count[i].f + Cat_Count[i].c;
        }else{
          for (j = 1; j < (timepoints+1); j++)
          {
            Right_Count_Censor[j] += Cat_Count[i].clist[j];
            Right_Count_Fail[j] += Cat_Count[i].flist[j];
          }
        }
      }

      if (LeftN > 0 && LeftN < node_n)
      {
        if (split_rule == 1)
          temp_score = logrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else if(split_rule == 2)
          temp_score = suplogrank(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftN, node_n, timepoints);
        else 
          temp_score = loglik(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, timepoints, varw);
          
      }
      
      if (temp_score > *score and split_rule<3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
          
          *score = temp_score;
          *cut = pack(ncat, goright);
      }
      if (temp_score < *score and split_rule==3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
        
        *score = temp_score;
        *cut = pack(ncat, goright);
      }
      
    }
  }

  NothingToFind: ;


  delete[] Cat_Count;

}

void Surv_One_Split_Cat_W(double* cut,
                          double* score,
                          const ivec useObs,
                          int node_n,
                          const colvec x,
                          const ivec Y, // y should be called by index i
                          const ivec Censor, // censor should be called by index i
                          const vec subjectweight,
                          int ncat,
                          int timepoints,
                          int split_gen,
                          int split_rule,
                          int nsplit,
                          int nmin,
                          int alpha)
{

  int i, k, j, temp_cat;
  double LeftWeights = 0;
  double RightWeights = 0;

  double temp_score = -1;
  
  ivec goright(ncat);
  goright.fill(0);
  ivec tempRight(ncat);
  tempRight.fill(0);
  vec Left_Count_Censor(timepoints+1);
  Left_Count_Censor.fill(0);
  vec Left_Count_Fail(timepoints+1);
  Left_Count_Fail.fill(0);
  vec Right_Count_Censor(timepoints+1);
  Right_Count_Censor.fill(0);
  vec Right_Count_Fail(timepoints+1);
  Right_Count_Fail.fill(0);


  // summerize all categories in this node
  SURVCAT_w* Cat_Count = new SURVCAT_w[ncat];

  for (i=0; i< ncat; i++)
  {
    Cat_Count[i].cat = i;
    Cat_Count[i].f = 0;
    Cat_Count[i].c = 0;
    Cat_Count[i].flist = vec(timepoints+1);
    Cat_Count[i].clist = vec(timepoints+1);
  }

  // record each category
  for (i=0; i<node_n; i++)
  {
    temp_cat = (int) x[useObs[i]] -1;

    if (Censor[i] == 0)
    {
      Cat_Count[temp_cat].clist[Y[i]] += subjectweight[useObs[i]];
      Cat_Count[temp_cat].c += subjectweight[useObs[i]];
    }else{
      Cat_Count[temp_cat].flist[Y[i]] += subjectweight[useObs[i]];
      Cat_Count[temp_cat].f += subjectweight[useObs[i]];
    }
  }

  // put nonzero categories to the front

  int true_ncat = ncat;
  for (i=0; i < true_ncat; i++)
  {
    if (Cat_Count[i].f + Cat_Count[i].c <= 0)
    {
      swap_SURVCAT_w(&Cat_Count[i], &Cat_Count[true_ncat-1]);
      true_ncat--;
      i--;
    }
  }

  // put categories with nonzero failures to the front

  int true_ncat_f = true_ncat;

  for (i =0; i < true_ncat_f; i++)
  {
    if (Cat_Count[i].f <= 0)
    {
      swap_SURVCAT_w(&Cat_Count[i], &Cat_Count[true_ncat_f-1]);
      true_ncat_f--;
      i--;
    }
  }

  if (true_ncat_f <= 1)
    goto NothingToFind;

  // if there are too many categories, save some computation by switching to random
  if (split_gen == 3 && true_ncat > 6)
  {
    split_gen = 1;
    nsplit = 64;
  }

  // if there are too many nsplit, switch to best
  if (split_gen < 3 && nsplit > pow(2, true_ncat-1) -1 )
  {
    split_gen = 3;
  }

  // random split
  if (split_gen == 1 || split_gen == 2)
  {

    for (k =0; k<nsplit; k++)
    {
      temp_score = -1;
      LeftWeights = 0;
      RightWeights = 0;

      // generate a random split that puts some cat into right
      for (i=0; i<true_ncat_f; i++)
        tempRight[i] = 1;

      permute_i(tempRight, true_ncat_f);

      // We should be checking if there will be enough observations for each side,
      // but I will skip it. Might need to fix it later

      // for the nonfailure categories, randomly assign
      for (i=true_ncat_f; i<ncat; i++)
        tempRight[i] = (int) unif_rand()>0.5;

      for (i = 0; i < true_ncat_f; i++)
      {
        if (tempRight[i] == 0) // go left
        {
          for (j = 1; j < (timepoints+1); j++)
          {
            Left_Count_Censor[j] += Cat_Count[i].clist[j];
            Left_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          LeftWeights += Cat_Count[i].f + Cat_Count[i].c;
        }else{ // go right
          for (j = 1; j < (timepoints+1); j++)
          {
            Right_Count_Censor[j] += Cat_Count[i].clist[j];
            Right_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          RightWeights += Cat_Count[i].f + Cat_Count[i].c;
        }
      }

      if (LeftWeights > 0 && RightWeights >0)
      {
        if (split_rule == 1)
          temp_score = logrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else if(split_rule == 2)
          temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else 
          temp_score = loglik_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, timepoints, varw);
      }

      if (temp_score > *score and split_rule<3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;

          *score = temp_score;
          *cut = pack(ncat, goright);
      }
      if (temp_score < *score and split_rule==3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
          
          *score = temp_score;
          *cut = pack(ncat, goright);
      }
    }
  }

  // best split

  if (split_gen == 3)
  {

    nsplit = pow(2, true_ncat-1) -1;

    for (k = 0; k < nsplit; k++)
    {

      LeftWeights = 0;
      RightWeights = 0;
      temp_score = -1;

      // fit the next possible split rule using binary vectors
      tempRight[0]++;

      for (i = 0; i < true_ncat - 1; i++)
      {
        if (tempRight[i] == 1)
        {
          tempRight[i] = 0;
          tempRight[i+1] ++;
        }
      }

      for (i = 0; i < true_ncat; i++)
      {
        if (tempRight[i] == 0)
        {
          for (j = 1; j < (timepoints+1); j++)
          {
            Left_Count_Censor[j] += Cat_Count[i].clist[j];
            Left_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          LeftWeights += Cat_Count[i].f + Cat_Count[i].c;
        }else{
          for (j = 1; j < (timepoints+1); j++)
          {
            Right_Count_Censor[j] += Cat_Count[i].clist[j];
            Right_Count_Fail[j] += Cat_Count[i].flist[j];
          }
          RightWeights += Cat_Count[i].f + Cat_Count[i].c;
        }
      }

      if (LeftWeights > 0 && RightWeights >0)
      {
        if (split_rule == 1)
          temp_score = logrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else if(split_rule == 2)
          temp_score = suplogrank_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, LeftWeights, LeftWeights+RightWeights, timepoints);
        else 
          temp_score = loglik_w(Left_Count_Fail, Left_Count_Censor, Right_Count_Fail, Right_Count_Censor, timepoints, varw);
      }

      if (temp_score > *score  and split_rule<3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;

          *score = temp_score;
          *cut = pack(ncat, goright);
      }
      if (temp_score < *score and split_rule==3)
      {
        for (i = 0; i< ncat; i ++)
          if (tempRight[i] == 1)
            goright[Cat_Count[i].cat] = 1;
          
          *score = temp_score;
          *cut = pack(ncat, goright);
      }
    }
  }

  NothingToFind:

  delete[] Cat_Count;

}
