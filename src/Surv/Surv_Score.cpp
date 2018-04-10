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
# include <Rmath.h>
#include <Rcpp.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

// log rank and sup log rank for equal weight version

double logrank(int* Left_Count_Fail, int* Left_Count_Censor, int* Right_Count_Fail, int* Right_Count_Censor, double LeftN, double AllN, int timepoints)
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
    else
      break;

    LeftN -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllN -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];

  }

  return tempscore;
}



double suplogrank(int* Left_Count_Fail, int* Left_Count_Censor, int* Right_Count_Fail, int* Right_Count_Censor, double LeftN, double AllN, int timepoints)
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
    else
      break;

    LeftN -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllN -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}


// log rank and sup log rank for subject weighted version

double logrank_w(double* Left_Count_Fail, double* Left_Count_Censor, double* Right_Count_Fail, double* Right_Count_Censor, double LeftWeights, double AllWeights, int timepoints)
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
    else // due to precision loss, this might be 0 already, so just stop here
      break;

    LeftWeights -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllWeights -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}



double suplogrank_w(double* Left_Count_Fail, double* Left_Count_Censor, double* Right_Count_Fail, double* Right_Count_Censor, double LeftWeights, double AllWeights, int timepoints)
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
    else // due to precision loss, this might be 0 already, so just stop here
      break;

    LeftWeights -= Left_Count_Fail[j] + Left_Count_Censor[j];
    AllWeights -= Left_Count_Fail[j] + Left_Count_Censor[j] + Right_Count_Fail[j] + Right_Count_Censor[j];
  }

  return tempscore;
}
