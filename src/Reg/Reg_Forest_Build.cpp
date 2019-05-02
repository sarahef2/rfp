//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Regression
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// my header file
# include "..//regForest.h"
# include "..//Utility//utility.h"

void regForestBuild(const std::vector< colvec > &X,
                    const vec &Y,
                    const ivec &Ncat,
                    const PARAMETERS* myPara,
                    const vec &subjectweight,
                    const ivec &subj_id,
                    const int &N,
                    vec &variableweight,
                    ivec &var_id,
                    const int &P,
                    TREENODE** Forest,
                    imat &ObsTrack,
                    imat &ObsTerminal,
                    std::vector< std::vector< ivec > > &NodeRegi,
                    vec &VarImp,
                    int use_cores)
{
  int ntrees = myPara->ntrees;
  int verbose = myPara->verbose;
  int replacement = myPara->replacement;
  int importance = myPara->importance;
  double resample_prob = myPara->resample_prob;
  int nimpute = myPara->nimpute;
  int size = (int) N*resample_prob;
  int nt;
  std::vector< mat > tree_matrix(ntrees);
  mat VarImp_store(ntrees, P);
  VarImp_store.fill(0);

  // start parallel trees

#pragma omp parallel for schedule(static) num_threads(use_cores)
  for (nt = 0; nt < ntrees; nt++) // fit all trees
  {
    R_DBP("Start tree %i \n", nt);
  }

}


