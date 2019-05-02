//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Regression
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// my header file
# include "Utility//utility.h"

#ifndef regForest_Fun
#define regForest_Fun

SEXP regForestFit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

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
                    int use_cores);
#endif
