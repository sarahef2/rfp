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
//# include <Rmath.h>
using namespace Rcpp;
using namespace arma;

// my header file
# include "utilities.h"


// debug

void printLog(const char* mode, const char* x, const int n1, const double n2)
{
  FILE* pFile = fopen("survForest_log.txt", mode);

  if(pFile != NULL)
    fprintf(pFile, x, n1, n2);

  fclose(pFile);
  return;
}

// get the list element named str, or return NULL

SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < Rf_length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
    return elmt;
}

// read-in all parameters

void copyParameters(PARAMETERS* myPara, SEXP list)
{
  myPara->N = INTEGER(getListElement(list, "n"))[0];
  myPara->P = INTEGER(getListElement(list, "p"))[0];
  myPara->Nfail = INTEGER(getListElement(list, "nfail"))[0];
  myPara->ntrees = INTEGER(getListElement(list, "ntrees"))[0];
  myPara->mtry = INTEGER(getListElement(list, "mtry"))[0];
  myPara->split_gen = INTEGER(getListElement(list, "split.gen"))[0];
  myPara->split_rule = INTEGER(getListElement(list, "split.rule"))[0];
  myPara->nsplit= INTEGER(getListElement(list, "nsplit"))[0];
  myPara->nmin = INTEGER(getListElement(list, "nmin"))[0];
  myPara->alpha = REAL(getListElement(list, "alpha"))[0];
  myPara->replacement = INTEGER(getListElement(list, "replacement"))[0];
  myPara->resample_prob = REAL(getListElement(list, "resample.prob"))[0];
  myPara->use_sub_weight = INTEGER(getListElement(list, "use.sub.weight"))[0];
  myPara->use_var_weight = INTEGER(getListElement(list, "use.var.weight"))[0];
  myPara->importance = INTEGER(getListElement(list, "importance"))[0];
  myPara->nimpute = INTEGER(getListElement(list, "nimpute"))[0];
  myPara->verbose = INTEGER(getListElement(list, "verbose"))[0];

  return;
}

// print all parameters


void printParameters(PARAMETERS* myPara)
{
  Rprintf("\nSurvival Forest (survForest) model fitting summary: ------------------\n");
  Rprintf("Data number of observations:                               n = %i \n", myPara->N);
  Rprintf("Data number of features:                                   p = %i \n", myPara->P);
  Rprintf("Data number of failures:                               nfail = %i \n", myPara->Nfail);
  Rprintf("Number of trees:                                      ntrees = %i \n", myPara->ntrees);
  Rprintf("Number of variables try at each split:                  mtry = %i \n", myPara->mtry);
  Rprintf("Splitting point generating method:                 split.gen = %s \n", myPara->split_gen == 1 ? "Random" : myPara->split_gen == 2 ? "Rank" : "Best");
  Rprintf("Splitting rule method:                            split.rule = %s \n", myPara->split_rule == 1 ? "logrank" : "suplogrank");
  if (myPara->split_gen < 3)
    Rprintf("Number of random splits:                              nsplit = %i \n", myPara->nsplit);
  Rprintf("Minimum terminal node size:                             nmin = %i \n", myPara->nmin);
  Rprintf("Minimum proportion of each child node:                 alpha = %1.2f%%\n", myPara->alpha*100);
  Rprintf("Sample with replacement:                         replacement = %s \n", myPara->replacement ? "Yes" : "No");
  Rprintf("Re-sampling proportion:                        resample.prob = %2.1f%% \n", myPara->resample_prob*100);
  Rprintf("Subject weights used:                         use.sub.weight = %s \n", myPara->use_sub_weight ? "Yes" : "No");
  Rprintf("Variable weights used:                        use.var.weight = %s \n", myPara->use_var_weight ? "Yes" : "No");
  Rprintf("Variable importance calculated:                   importance = %s \n", myPara->importance ? "Yes" : "No");
  Rprintf("Number of imputations for importance:                nimpute = %i \n", myPara->nimpute);
  Rprintf("------------------------------------------------------------------------\n");
}


// tree node functions

// tree functions
int TreeSize(TREENODE *root)
{
  if (root->Var == -1)
    return 1;
  else{
    int count = 1;
    count += TreeSize(root->Left);
    count += TreeSize(root->Right);
    return count;
  }
}

void CheckVar(mat tree_matrix_nt, int j, bool& Check){
  int tree_depth = tree_matrix_nt.n_rows;
  Check = false;
  for(int i=0; i<tree_depth; i++){
    if((int) tree_matrix_nt(i,0)-1 == j){
      Check = true;
      break;
    }
  }
}






// other utility functions

void standardize(
    vec &x, int n)
{
  double sumx=0;
  int i;

  for (i=0; i<n; i++)
    sumx += x[i];
  
  for (i=0; i<n; i++)
    x[i] = x[i]/sumx;
}

// random sample a value
int random_in_range(int min, int max)
{
  if (min == max)
    return min;

  double u;
  do {u = R::runif((double) min, (double) max);} while (u <= min || u >= max);
  return (int) u; // generates integers from min to max-1
}

// sample with weight from the tail of a vector and rotate the sampled one to the front

int sample_rotate(ivec &index, vec &weights, int start, int end)
{
  // sample without replacement, will put sampled one to the front
  // will order weights vector in the same way
  int a = random_in_range(start, end);
  int temp = index[a];
  double tempw = weights[a];

  index[a] = index[start];
  index[start] = temp;

  weights[a] = weights[start];
  weights[start] = tempw;

  return temp;
}

// random sample with weight

int weighted_sample(const double* x, int n)
{
  double a = R::runif(0, 1);

  int i;

  for (i = 0; i< n; i++)
  {
    a -= x[i];

    if (a <= 0)
      return i;
  }

  if (a < WeightTH)
    return weighted_sample(x, n);
  else
    ::Rf_error("weighted vector is not properly normalized");
}


// min
double dmin(double a, double b)
{
  if (a <= b)
    return a;
  return b;
}

// max
double dmax(double a, double b)
{
  if (a>=b)
    return a;
  return b;
}

// min
int imin(int a, int b)
{
  if (a <= b)
    return a;
  return b;
}

// max
int imax(int a, int b)
{
  if (a>=b)
    return a;
  return b;
}


void swap_d(double* x, int i, int j)
{
  double temp;
  temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}

void swap_i(int* x, int i, int j)
{
  int temp;
  temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}

// sort by index

// quick sort function, for larger size

int partition_d(arma::vec& keys,
                int low,
                int high,
                arma::ivec& index)
{
  double pivot = keys[high];    // pivot
  double dtemp;
  size_t itemp;

  int i = (low - 1);            // Index of smaller element

  for (int j = low; j <= high- 1; j++)
  {
    // If current element is smaller than or equal to pivot
    if (keys[j] <= pivot)
    {
      i++;    // increment index of smaller element

      // swap index
      itemp = index[i];
      index[i] = index[j];
      index[j] = itemp;

      // swap keys
      dtemp = keys[i];
      keys[i] = keys[j];
      keys[j] = dtemp;
    }
  }


  // swap index
  itemp = index[i + 1];
  index[i + 1] = index[high];
  index[high] = itemp;

  // swap keys
  dtemp = keys[i + 1];
  keys[i + 1] = keys[high];
  keys[high] = dtemp;

  return (i + 1);
}

void qSort_dindex(arma::vec& keys,
                  int low,
                  int high,
                  arma::ivec& index)
{
  if (low < high)
  {
    /* par is partitioning index, keys[p] is now at right side */
    int par = partition_d(keys, low, high, index);

    // Separately sort left and right
    qSort_dindex(keys, low, par - 1, index);
    qSort_dindex(keys, par + 1, high, index);
  }
}

// for small size, use this

void iSort_index(arma::vec& keys,
                 int low,
                 int high,
                 arma::ivec& index)
{
  bool bool_sorted = false;
  size_t itemp;
  double dtemp;

  //check whether all keys are in the correct order

  while (bool_sorted == false)
  {
    bool_sorted = true;

    for (int i = low; i < high; i++)
    {
      //if next value is lower

      if (keys[i] > keys[i + 1])
      {
        //swap + key index
        itemp = index[i];
        index[i] = index[i + 1];
        index[i + 1] = itemp;

        dtemp = keys[i];
        keys[i] = keys[i + 1];
        keys[i + 1] = dtemp;

        bool_sorted = false;
      }
    }
  }
}


// quick sort function, for integer index

int partition_i(arma::ivec& keys,
                int low,
                int high,
                arma::ivec& index)
{
  int pivot = keys[high];    // pivot
  int temp;

  int i = (low - 1);            // Index of smaller element

  for (int j = low; j <= high- 1; j++)
  {
    // If current element is smaller than or equal to pivot
    if (keys[j] <= pivot)
    {
      i++;    // increment index of smaller element

      // swap index
      temp = index[i];
      index[i] = index[j];
      index[j] = temp;

      // swap keys
      temp = keys[i];
      keys[i] = keys[j];
      keys[j] = temp;
    }
  }


  // swap index
  temp = index[i + 1];
  index[i + 1] = index[high];
  index[high] = temp;

  // swap keys
  temp = keys[i + 1];
  keys[i + 1] = keys[high];
  keys[high] = temp;

  return (i + 1);
}

void qSort_iindex(arma::ivec& keys,
                  int low,
                  int high,
                  arma::ivec& index)
{
  if (low < high)
  {
    /* par is partitioning index, keys[p] is now at right side */
    int par = partition_i(keys, low, high, index);

    // Separately sort left and right
    qSort_iindex(keys, low, par - 1, index);
    qSort_iindex(keys, par + 1, high, index);
  }
}

// categorical variable swap

void swap_SURVCAT(SURVCAT* a, SURVCAT* b)
{
  SURVCAT temp = *a;
  *a = *b;
  *b = temp;
}

void swap_SURVCAT_w(SURVCAT_w* a, SURVCAT_w* b)
{
  SURVCAT_w temp = *a;
  *a = *b;
  *b = temp;
}

// permutation

void permute_i(ivec &x, int n)
{
  int i;
  int j;
  int temp;

  for (i = 0; i<n-1; i++)
  {
    j = random_in_range(i, n);
    temp = x[i];
    x[i] = x[j];
    x[j] = temp;
  }
}

// permutation (float)

void permute(vec &x, int n)
{
  int i;
  int j;
  int temp;
  
  for (i = 0; i<n-1; i++)
  {
    j = random_in_range(i, n);
    temp = x[i];
    x[i] = x[j];
    x[j] = temp;
  }
}

// variable pack

double pack(const int nBits, const ivec bits) // from Andy's rf package
{
  int i;
  double value = bits[nBits - 1];

  for (i = nBits - 2; i >= 0; i--)
    value = 2.0*value + bits[i];

  return(value);
}

void unpack(const double pack, const int nBits, ivec bits) // from Andy's rf package
{
  int i;
  double x = pack;
  for (i = 0; i < nBits; ++i)
  {
    bits[i] = ((unsigned long) x & 1) ? 1 : 0;
    x /= 2;
  }
}

int unpack_goright(double pack, const int cat)
{
  int i;

  for (i = 0; i < cat; i++) pack /= 2;

  return(((unsigned long) pack & 1) ? 1 : 0);
}



/*
// test the pack
int nc = 10;
int* testpack = (int *) malloc(nc * sizeof(int));
double packed;

for (int k=0; k < 10; k++)
{
memset(testpack, 0, nc*sizeof(int));

for (int i=0; i<nc; i++)
testpack[i] = 1;

memset(testpack, 0, random_in_range(1, nc)*sizeof(int));
permute_i(testpack, nc);

for (int i=0; i<nc; i++)
R_DBP(" %i", testpack[i]);

R_DBP(" \n");

packed = pack(nc, testpack);

for (int i=0; i<nc; i++)
R_DBP(" %i", unpack_goright(packed, i));

R_DBP(" \n\n\n");

}
*/

