//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Regression
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// my header file
# include "regForest.h"
# include "Utility//utility.h"

// main function
// [[Rcpp::export()]]
List regForestFit(arma::mat datasetX_R,
                  arma::vec datasetY_R,
                  arma::ivec ncat_R,
                  arma::vec subjectweight_R,
                  arma::vec variableweight_R,
                  List parameters_R,
                  int usecores_R)
{
  // copy parameters and check
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  if (myPara->verbose) printParameters(myPara);

  // check number of cores
  int use_cores = usecores_R;
  checkCores(use_cores, myPara->verbose);

  // create data objects
  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;

  // get X, Y
  std::vector< colvec > X(P, colvec(N));

  for (int j = 0; j < P; j++) {
    X[j] = conv_to<colvec>::from(datasetX_R.col(j));
  }

  const vec Y = datasetY_R;
  const ivec Ncat = ncat_R;

  vec subjectweight = subjectweight_R;
  standardize(subjectweight, N);	// this could change the input data due to precision loss...

  vec variableweight = variableweight_R;
  standardize(variableweight, P);	// this could change the input data due to precision loss...

  // initiate tree and other objects

  TREENODE **Forest = new TREENODE*[ntrees];

  imat ObsTrack(N,ntrees);
  ObsTrack.fill(0);

  imat ObsTerminal(N,ntrees);
  ObsTerminal.fill(-1);

  std::vector< std::vector< ivec > > NodeRegi(ntrees);

  vec VarImp(P);
  VarImp.fill(0);

  ivec subj_id(N);
  for (int i = 0; i  < N; i++) subj_id(i) = i;

  ivec var_id(P);
  for (int i = 0; i  < P; i++) var_id(i) = i;

  // start to fit the model
  regForestBuild((const std::vector< colvec >) X,
                 (const vec) Y,
                 (const ivec) Ncat,
                 (const PARAMETERS*) myPara,
                 (const vec) subjectweight,
                 (const ivec) subj_id,
                 (const int) N,
                 variableweight,
                 var_id,
                 (const int) P,
                 Forest,
                 ObsTrack,
                 ObsTerminal,
                 NodeRegi,
                 VarImp,
                 use_cores);

  int TreeWidth = 4;

  // return subjects to R

  mat FittedTree;
  List FittedForest(ntrees);

  int Node;
  int TreeLength;

  for (int nt = 0; nt<ntrees; nt++)
  {
    TreeLength = TreeSize(Forest[nt]);
    Node = 0;

    FittedTree.set_size(TreeLength, TreeWidth);

    //Convert tree structure nodes into matrix to send to R matrix
    //Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);

    FittedForest(nt)=FittedTree;
  }

  List ReturnList;

  ReturnList["FittedForest"] = FittedForest;
  ReturnList["ObsTrack"] = ObsTrack;
  ReturnList["NodeRegi"] = NodeRegi;
  ReturnList["Importance"] = VarImp;
  ReturnList["ObsTerminal"] = ObsTerminal;

  delete[] myPara;
  delete[] Forest;

  return ReturnList;
}
