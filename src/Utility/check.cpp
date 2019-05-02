//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Utility Functions: check
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// header file
# include "utility.h"

// check cores

void checkCores(int &use_cores, int verbose)
{
  use_cores = max(1, use_cores);

  if (use_cores > 0) OMPMSG(1);

  int haveCores = omp_get_max_threads();

  if(use_cores > haveCores)
  {
    if (verbose) Rprintf("Do not have %i cores, use maximum %i cores. \n", use_cores, haveCores);
    use_cores = haveCores;
  }
}

// standardize

void standardize(vec &x, int n)
{
  double sumx=0;
  int i;

  for (i=0; i<n; i++)
    sumx += x[i];

  for (i=0; i<n; i++)
    x[i] = x[i]/sumx;
}


// copy, check, print parameters

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
  myPara->nmin_control = INTEGER(getListElement(list, "nmin.control"))[0];
  myPara->nmin_failure = INTEGER(getListElement(list, "nmin.failure"))[0];
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
  Rprintf("Control terminal node size:                     nmin.control = %s \n", myPara->nmin_control ? "Yes" : "No");
  Rprintf("Control terminal node failure count:            nmin.failure = %s \n", myPara->nmin_failure ? "Yes" : "No");
  Rprintf("Minimum proportion of each child node:                 alpha = %1.2f%%\n", myPara->alpha*100);
  Rprintf("Sample with replacement:                         replacement = %s \n", myPara->replacement ? "Yes" : "No");
  Rprintf("Re-sampling proportion:                        resample.prob = %2.1f%% \n", myPara->resample_prob*100);
  Rprintf("Subject weights used:                         use.sub.weight = %s \n", myPara->use_sub_weight ? "Yes" : "No");
  Rprintf("Variable weights used:                        use.var.weight = %s \n", myPara->use_var_weight ? "Yes" : "No");
  Rprintf("Variable importance calculated:                   importance = %s \n", myPara->importance ? "Yes" : "No");
  Rprintf("Number of imputations for importance:                nimpute = %i \n", myPara->nimpute);
  Rprintf("------------------------------------------------------------------------\n");
}

// tree functions

// get number of nodes
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

