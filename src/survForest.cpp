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
  //PARAMETERS *myPara = (PARAMETERS*) malloc(sizeof(PARAMETERS));
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);
  //Rcout << "Got to A " << std::endl;
  
  if (myPara->verbose) printParameters(myPara);

  //int use_cores = INTEGER(usecores_R)[0];
  int use_cores = usecores_R;
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  //Rcout << "Got to B " << std::endl;
  //// create data objects

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  int nmin = myPara->nmin;

  int i;
  int j;
  int nt = 0;

  // get X, Y, Censor and Treatment

  std::vector< colvec > X(P, colvec(N));
  
  for (j = 0; j < P; j++) {
    X[j] = conv_to<colvec>::from(datasetX_R.col(j));
  };
  
  //Rcout << "Got to C " << std::endl;
  const ivec Y = datasetY_R;
  const ivec Censor= datasetCensor_R;
  const ivec Ncat=ncat_R;
  const vec Interval=interval_R;

  vec subjectweight = subjectweight_R;
  standardize(subjectweight, N);	// this could change the input data due to precision loss...

  //Rcout << "Got to D " << std::endl;
  vec variableweight = variableweight_R;
  standardize(variableweight, P);	// this could change the input data due to precision loss...

  // initiate tree

  TREENODE **Forest = new TREENODE*[ntrees];

  imat ObsTrack(N,ntrees);
  ObsTrack.fill(0);
  //imat NodeRegi(N,ntrees);
  //std::vector< List > NodeRegi(ntrees);
  std::vector< std::vector< ivec > > NodeRegi(ntrees);
  //NodeRegi.fill(0);
  mat VarImp(ntrees,N);
  //Rcout << "Got to E " << std::endl;
  
  ivec subj_id(N);
  for (i = 0; i  < N; i++) subj_id(i) = i;
  
  ivec var_id(N);
  for (i = 0; i  < P; i++) var_id(i) = i;  
  
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
                  NodeRegi,
                  VarImp,
                  use_cores);
  
  //Rcout << "Got to F " << std::endl;
  
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

    //PROTECT(FittedTree = Rf_allocMatrix(REALSXP, TreeLength, TreeWidth));
    
    FittedTree.set_size(TreeLength,TreeWidth);

    //Convert tree structure nodes into matrix to send to R matrix
    Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);

    //THESE TWO LINES NEED TO BE FIXED
    //FittedTreeR = as< NumericMatrix >( &FittedTree );
    //Rf_setAttrib(FittedTreeR, R_DimNamesSymbol, TreeNames);

    //SET_VECTOR_ELT(FittedForest, nt, FittedTree);
    FittedForest(nt)=FittedTree;
  }
  //Rcout << "Got to G " << std::endl;
  
  List ReturnList;

  ReturnList["FittedForest"] = FittedForest;
  ReturnList["ObsTrack"] = ObsTrack;
  ReturnList["NodeRegi"] = NodeRegi;

  delete[] myPara;
  delete[] Forest;
  //Rcout << "Got to H " << std::endl;
  
  return ReturnList;
}


// print
// [[Rcpp::export()]]
void survForestPrint(List parameters_R)
{
  //PARAMETERS *myPara = malloc(sizeof(PARAMETERS));
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);
  printParameters(myPara);
  //free(myPara);
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
  //PARAMETERS *myPara = (PARAMETERS*) malloc(sizeof(PARAMETERS));
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  //// create data objects
  //int use_cores = INTEGER(usecores_R)[0];
  int use_cores = usecores_R;
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  int Nfail = myPara->Nfail;

  //SEXP dataX_dim = Rf_getAttrib(testsetX_R, R_DimSymbol);
  //int testN = INTEGER(dataX_dim)[0];
  int testN = testsetX_R.n_rows;

  int i;
  int j;
  int nt = 0;

  // get X, Y, Censor

  //double **testX = (double **) malloc(P * sizeof(double *));
  //double **testX = new double*[P];
  //if (testX == NULL) error("Unable to malloc testX");
  //if (testX == NULL) stop("Unable to malloc testX");
  //for (j = 0; j < P; j++)
  //  testX[j] = &REAL(testsetX_R)[j*testN];

  std::vector< colvec > testX(P, colvec(testN));

  for (j = 0; j < P; j++) {
    testX[j] = conv_to<colvec>::from(testsetX_R.col(j));
  };

  //const int *Y = INTEGER(datasetY_R);
  //const int *Censor = INTEGER(datasetCensor_R);
  //const int *Ncat = INTEGER(datasetNcat_R);
  const ivec Y = datasetY_R;
  const ivec Censor= datasetCensor_R;
  const ivec Ncat=datasetNcat_R;

  //const double *subjectweight = REAL(subjectweight_R);
  const vec subjectweight = subjectweight_R;

  // get tree matrix

  int TreeWidth = 4;
  int TreeLength;

  //double ***tree_matrix = (double ***) malloc(ntrees * sizeof(double **));
  //double ***tree_matrix = new double**[ntrees];
  std::vector< mat > tree_matrix(ntrees);
  std::vector< std::vector< ivec > > NodeRegi(ntrees);
  //if (tree_matrix == NULL) error("Unable to malloc for tree_matrix");
  //if (tree_matrix == NULL) stop("Unable to malloc for tree_matrix");

  for (nt=0; nt<ntrees; nt++)
  {
    //tree_matrix[nt] = (double **) malloc(TreeWidth  * sizeof(double*));
    //tree_matrix[nt] = new double*[TreeWidth];
    //if (tree_matrix[nt] == NULL) error("Unable to malloc for tree_matrix");
    //if (tree_matrix[nt] == NULL) stop("Unable to malloc for tree_matrix");

    //TreeLength = INTEGER(Rf_getAttrib(VECTOR_ELT(FittedForest_R, nt), R_DimSymbol))[0];
    //TreeLength=FittedForest_R[nt].size;
    
    tree_matrix[nt] = as<mat>(FittedForest_R[nt]);
    for(i = 0; i < Rcpp::as<List>(NodeRegiMat_R[nt]).size(); i++){
      NodeRegi[nt].push_back(as<ivec>(as<List>(NodeRegiMat_R[nt])[i]));
    }
    //NodeRegi[nt] = NodeRegiMat_R[nt];
    //for (i = 0; i < TreeWidth; i++)
    //  tree_matrix[nt][i] = &REAL(VECTOR_ELT(FittedForest_R, nt))[i*TreeLength];
  }

  // get ObsTrack

  //int **ObsTrack = (int **) malloc(ntrees * sizeof(int *));
  //int **ObsTrack = new int*[ntrees];
  imat ObsTrack(N,ntrees);
  //if (ObsTrack == NULL) error("Unable to malloc ObsTrack");
  //if (ObsTrack == NULL) stop("Unable to malloc ObsTrack");
  //for (nt = 0; nt < ntrees; nt++)
  //  ObsTrack[nt] = &INTEGER(ObsTrackMat_R)[nt*N];
  ObsTrack=ObsTrackMat_R;

  // get NodeRegi
  //int **NodeRegi = (int **) malloc(ntrees * sizeof(int *));
  //int **NodeRegi = new int*[ntrees];
  //imat NodeRegi(N,ntrees);
  //if (NodeRegi == NULL) error("Unable to malloc NodeRegi");
  //if (NodeRegi == NULL) stop("Unable to malloc NodeRegi");
  //for (nt = 0; nt < ntrees; nt++)
  //  NodeRegi[nt] = &INTEGER(NodeRegiMat_R)[nt*N];
  //NodeRegi=NodeRegiMat_R;

  //SEXP SurvMat;
  //mat SurvMat(Nfail+1,testN);
  //PROTECT(SurvMat = Rf_allocMatrix(REALSXP, Nfail+1, testN));

  //double **surv_matrix = (double **) malloc(testN * sizeof(double *));
  //double **surv_matrix = new double*[testN];
  mat surv_matrix(testN,Nfail+1);
  surv_matrix.fill(0);
  //if (surv_matrix == NULL) error("Unable to malloc surv_matrix");
  //if (surv_matrix == NULL) stop("Unable to malloc surv_matrix");
  
  //The following loop should be unnecessary after filling the matrix with 0's
  //for (i = 0; i < testN; i++)
  //{
    //surv_matrix[i] = &REAL(SurvMat)[i*(Nfail+1)];//What is this line supposed to do?
    //for(j=0; j<= Nfail; j++)
      //surv_matrix(i,j) = 0;
      //surv_matrix[i][j] = 0;
  //}

  PredictSurvivalKernel((const std::vector< colvec >) testX,
                        Y,
                        Censor,
                        Ncat,
                        subjectweight,
                        (const std::vector< mat >) tree_matrix,
                        (const imat) ObsTrack,
                        (const std::vector< std::vector< ivec > >) NodeRegi,
                        surv_matrix,
                        (const PARAMETERS*) myPara,
                        testN,
                        use_cores);


  //free(testX);
  //delete[] testX;

  //for (nt = 0; nt < ntrees; nt++)
    //free(tree_matrix[nt]);
    //delete[] tree_matrix[nt];
  //free(tree_matrix);
  //delete[] tree_matrix;

  //free(NodeRegi);
  //delete[] NodeRegi;
  //free(ObsTrack);
  //delete[] ObsTrack;
  //free(surv_matrix);
  //delete[] surv_matrix;

  //UNPROTECT(1);

  return surv_matrix;//SurvMat;
 }
 
 
