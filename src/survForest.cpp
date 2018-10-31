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
using namespace Rcpp;
using namespace arma;


// my header file
# include "survForest.h"
# include "utilities.h"

// main function
// [[Rcpp::export()]]
List survForestFit(arma::mat datasetX_R,
                   arma::ivec datasetY_R,
                   arma::ivec datasetCensor_R,
                   arma::ivec ncat_R,
                   arma::vec interval_R,
                   arma::vec subjectweight_R,
                   arma::vec variableweight_R,
                   List parameters_R,
                   int usecores_R)
{
  // copy parameters and check
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  if (myPara->verbose) printParameters(myPara);

  int use_cores = usecores_R;
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  //// create data objects

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  //int nmin = myPara->nmin;
  int Nfail = myPara->Nfail; //Used in variable importance calculations
  
  int i;
  int j;
  int nt = 0;
  int counter = 0;

  // get X, Y, Censor and Treatment

  std::vector< colvec > X(P, colvec(N));
  
  for (j = 0; j < P; j++) {
    X[j] = conv_to<colvec>::from(datasetX_R.col(j));
  };
  
  const ivec Y = datasetY_R;
  const ivec Censor= datasetCensor_R;
  const ivec Ncat=ncat_R;
  const vec Interval=interval_R;

  vec subjectweight = subjectweight_R;
  standardize(subjectweight, N);	// this could change the input data due to precision loss...

  vec variableweight = variableweight_R;
  standardize(variableweight, P);	// this could change the input data due to precision loss...

  // initiate tree

  TREENODE **Forest = new TREENODE*[ntrees];

  imat ObsTrack(N,ntrees);
  ObsTrack.fill(0);
  imat ObsTerminal(N,ntrees);
  ObsTerminal.fill(-1);
  std::vector< std::vector< ivec > > NodeRegi(ntrees);
  vec VarImp(P);

  ivec subj_id(N);
  for (i = 0; i  < N; i++) subj_id(i) = i;
  
  ivec var_id(P);
  for (i = 0; i  < P; i++) var_id(i) = i;  
  
  mat oob_surv_matrix(N,Nfail+1);
  oob_surv_matrix.fill(0);
  
  vec oob_residuals(N);
  oob_residuals.fill(0);

  // start to fit the model
  survForestBuild((const std::vector<  colvec >) X,
                  (const ivec) Y,
                  (const ivec) Censor,
                  (const ivec) Ncat,
                  (const vec) Interval,
                  (const PARAMETERS*) myPara,
                  subjectweight,
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
                  use_cores,
                  oob_surv_matrix,
                  oob_residuals,
                  counter);
  
  int TreeWidth = 4;
  
  // return subjects to R 
  
  mat FittedTree;
  List FittedForest(ntrees);
  
  int Node;
  int TreeLength;

  for (nt = 0; nt<ntrees; nt++)
  {
    TreeLength = TreeSize(Forest[nt]);
    Node = 0;

    FittedTree.set_size(TreeLength,TreeWidth);

    //Convert tree structure nodes into matrix to send to R matrix
    Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);

    FittedForest(nt)=FittedTree;
  }

  List ReturnList;

  ReturnList["FittedForest"] = FittedForest;
  ReturnList["ObsTrack"] = ObsTrack;
  ReturnList["NodeRegi"] = NodeRegi;
  ReturnList["Importance"] = VarImp;
  ReturnList["counter"] = counter;
  ReturnList["ObsTerminal"] = ObsTerminal;
  ReturnList["OobSurvival"] = oob_surv_matrix;
  ReturnList["Residuals"] = oob_residuals;
  
  delete[] myPara;
  delete[] Forest;

  return ReturnList;
}


// print
// [[Rcpp::export()]]
void survForestPrint(List parameters_R)
{
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);
  printParameters(myPara);
  delete[] myPara;
}


// predict
// [[Rcpp::export()]]
 arma::mat survForestPredict(arma::mat testsetX_R,
                        List FittedForest_R,
                        arma::imat datasetY_R,
                        arma::ivec datasetCensor_R,
                        arma::ivec datasetNcat_R,
                        arma::vec subjectweight_R,
                        arma::imat ObsTrackMat_R,
                        List NodeRegiMat_R,
                        List parameters_R,
                        int usecores_R)
 {
//   
  // copy parameters and check
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  //// create data objects
  int use_cores = usecores_R;
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  int Nfail = myPara->Nfail;

  int testN = testsetX_R.n_rows;

  int i;
  int j;
  int nt = 0;

  // get X, Y, Censor

  std::vector< colvec > testX(P, colvec(testN));

  for (j = 0; j < P; j++) {
    testX[j] = conv_to<colvec>::from(testsetX_R.col(j));
  };

  const ivec Y = datasetY_R;
  const ivec Censor= datasetCensor_R;
  const ivec Ncat=datasetNcat_R;

  const vec subjectweight = subjectweight_R;

  // get tree matrix

  //int TreeWidth = 4;
  //int TreeLength;

  std::vector< mat > tree_matrix(ntrees);
  std::vector< std::vector< ivec > > NodeRegi(ntrees);

  for (nt=0; nt<ntrees; nt++)
  {
    tree_matrix[nt] = as<mat>(FittedForest_R[nt]);
    for(i = 0; i < Rcpp::as<List>(NodeRegiMat_R[nt]).size(); i++){
      NodeRegi[nt].push_back(as<ivec>(as<List>(NodeRegiMat_R[nt])[i]));
    }
  }

  // get ObsTrack

  imat ObsTrack(N,ntrees);
  ObsTrack=ObsTrackMat_R;

  // get NodeRegi
  mat surv_matrix(testN,Nfail+1);
  ivec tmp(testN);
  tmp.fill(0);
  
  ivec subj_id(testN);
  for (i = 0; i  < testN; i++) subj_id(i) = i;

  imat ObsTerminal(N,ntrees);
  ObsTerminal.fill(-1);
  PredictSurvivalKernel((const std::vector< colvec >) testX,
                        Y,
                        Censor,
                        Ncat,
                        subjectweight,
                        (const std::vector< mat >) tree_matrix,
                        (const imat) ObsTrack,
                        ObsTerminal,
                        (const std::vector< std::vector< ivec > >) NodeRegi,
                        surv_matrix,
                        (const PARAMETERS*) myPara,
                        testN,
                        (const ivec) subj_id,
                        (const int) -1,
                        (const ivec) tmp,
                        (int) 0,
                        false,
                        use_cores);


  return surv_matrix;//SurvMat;
 }
 
 
