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

//# include <RcppArmadillo.h>
# include <Rcpp.h>
//# include <Rdefines.h>
//# include <Rinternals.h>
//# include <R.h>
using namespace Rcpp;


// my header file
# include "survForest.h"
# include "utilities.h"

// main function
// [[Rcpp::export()]]
SEXP survForestFit(SEXP datasetX_R,
                   SEXP datasetY_R,
                   SEXP datasetCensor_R,
                   SEXP ncat_R,
                   SEXP interval_R,
                   SEXP subjectweight_R,
                   SEXP variableweight_R,
                   SEXP parameters_R,
                   SEXP usecores_R)
{
  // copy parameters and check
  //PARAMETERS *myPara = (PARAMETERS*) malloc(sizeof(PARAMETERS));
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  if (myPara->verbose) printParameters(myPara);

  int use_cores = INTEGER(usecores_R)[0];
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  //// create data objects

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  int nmin = myPara->nmin;

  int i;
  int j;
  int nt = 0;

  // get X, Y, Censor and Treatment

  //double **X = (double **) malloc(P * sizeof(double *));
  double **X = new double*[P];
  //if (X == NULL) error("Unable to malloc X");
  if (X == NULL) stop("Unable to malloc X");
  for (j = 0; j < P; j++)
    X[j] = &REAL(datasetX_R)[j*N];

  const int *Y = INTEGER(datasetY_R);
  const int *Censor = INTEGER(datasetCensor_R);
  const int *Ncat = INTEGER(ncat_R);
  const double *Interval = REAL(interval_R);

  // this is a subject id indicator
  //int* subj_id = (int *) malloc (N * sizeof(int));
  int* subj_id = new int[N];
  //if (subj_id == NULL) error("Unable to malloc subj_id");
  if (subj_id == NULL) stop("Unable to malloc subj_id");
  for (i = 0; i < N; i++)
    subj_id[i] = i;

  // this is a variable id indicator
  //int* var_id = (int *) malloc (P * sizeof(int));
  int* var_id = new int[P];
  //if (var_id == NULL) error("Unable to malloc var_id");
  if (var_id == NULL) stop("Unable to malloc var_id");
  for (i = 0; i < P; i++)
    var_id[i] = i;

  // get weights (may not be used)
  //const double *subjectweight = REAL(subjectweight_R);
  double *subjectweight = REAL(subjectweight_R);
  standardize(subjectweight, N);	// this could change the input data due to precision loss...

  double *variableweight = REAL(variableweight_R);
  standardize(variableweight, P);	// this could change the input data due to precision loss...

  // initiate tree

  //TREENODE **Forest = (TREENODE **) malloc(ntrees * sizeof(TREENODE *));
  TREENODE **Forest = new TREENODE*[ntrees];

  // for tracking observations in each tree
  //int **ObsTrack = (int **) malloc(ntrees * sizeof(int *));
  int **ObsTrack = new int*[ntrees];
  //if (ObsTrack == NULL) error("Unable to malloc for ObsTrack");
  if (ObsTrack == NULL) stop("Unable to malloc for ObsTrack");

  for (nt = 0; nt<ntrees; nt++)
  {
    //ObsTrack[nt] = (int *) calloc(N, sizeof(int));
    ObsTrack[nt] = new int[N];
    //if (ObsTrack[nt] == NULL) error("Unable to calloc for ObsTrack");
    if (ObsTrack[nt] == NULL) stop("Unable to calloc for ObsTrack");
  }

  // for tracking the terminal nodes of each observation
  //int **NodeRegi = (int **) malloc(ntrees * sizeof(int *));
  int **NodeRegi = new int*[ntrees];
  //if (NodeRegi == NULL) error("Unable to malloc for NodeRegi");
  if (NodeRegi == NULL) stop("Unable to malloc for NodeRegi");

  for (nt = 0; nt<ntrees; nt++)
  {
    //NodeRegi[nt] = (int *) calloc(N, sizeof(int));
    NodeRegi[nt] = new int[N];
    //if (NodeRegi[nt] == NULL) error("Unable to calloc for NodeRegi");
    if (NodeRegi[nt] == NULL) stop("Unable to calloc for NodeRegi");
  }

  // for variable importance
  //double **VarImp = (double **) malloc(ntrees * sizeof(double *));
  double **VarImp = new double*[ntrees];
  //if (VarImp == NULL) error("Unable to calloc for mse recording");
  if (VarImp == NULL) stop("Unable to calloc for mse recording");

  for (nt=0; nt<ntrees; nt++)
  {
    //VarImp[nt] = (double *) calloc((P+1), sizeof(double));
    VarImp[nt] = new double[P+1];
    //if (VarImp[nt] == NULL) error("Unable to calloc for mse recording");
    if (VarImp[nt] == NULL) stop("Unable to calloc for mse recording");
  }

  // start to fit the model
  survForestBuild((const double**) X,
                  (const int*) Y,
                  (const int*) Censor,
                  (const int*) Ncat,
                  (const double*) Interval,
                  (const PARAMETERS*) myPara,
                  subjectweight,
                  (const int*) subj_id,
                  (const int) N,
                  variableweight,
                  var_id,
                  (const int) P,
                  Forest,
                  ObsTrack,
                  NodeRegi,
                  VarImp,
                  use_cores);
  
  delete[] var_id;

  // create objects to return

  int TreeWidth = 4;

  SEXP TreeNames; // both column and row names
  SEXP ColNames;  // just column names

  PROTECT(TreeNames = Rf_allocVector(VECSXP, 2));
  PROTECT(ColNames = Rf_allocVector(VECSXP, TreeWidth));

  //set column names

  SET_VECTOR_ELT(ColNames, 0, Rf_mkChar("SplitVar"));    // splitting variable
  SET_VECTOR_ELT(ColNames, 1, Rf_mkChar("SplitValue"));  // splitting value
  SET_VECTOR_ELT(ColNames, 2, Rf_mkChar("NextLeft"));    // left daughter node
  SET_VECTOR_ELT(ColNames, 3, Rf_mkChar("NextRight"));   // right daughter node
  //SET_VECTOR_ELT(ColNames, 4, Rf_mkChar("NodeID"));   // Node ID

  SET_VECTOR_ELT(TreeNames, 1, ColNames);

  // converting each tree into R subject
  SEXP FittedTree;
  SEXP FittedForest;
  PROTECT(FittedForest = Rf_allocVector(VECSXP, ntrees));
  int Node;
  int TreeLength;

  for (nt = 0; nt<ntrees; nt++)
  {
    TreeLength = TreeSize(Forest[nt]);
    Node = 0;

    PROTECT(FittedTree = Rf_allocMatrix(REALSXP, TreeLength, TreeWidth));

    Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);

    Rf_setAttrib(FittedTree, R_DimNamesSymbol, TreeNames);

    SET_VECTOR_ELT(FittedForest, nt, FittedTree);
  }

  delete[] Forest;
  
  SEXP ObsTrackMat;
  PROTECT(ObsTrackMat = Rf_allocMatrix(INTSXP, N, ntrees));

  for (nt = 0; nt<ntrees; nt++)
  {
    for (i = 0; i < N; i++)
      INTEGER(ObsTrackMat)[i + nt*N] = ObsTrack[nt][i];

    delete[] ObsTrack[nt];
  }
  delete[] ObsTrack;

  SEXP NodeRegiMat;
  PROTECT(NodeRegiMat = Rf_allocMatrix(INTSXP, N, ntrees));

  for (nt = 0; nt<ntrees; nt++)
  {
    for (i = 0; i < N; i++)
      INTEGER(NodeRegiMat)[i + nt*N] = NodeRegi[nt][i];

    delete[] NodeRegi[nt];
  }
  delete[] NodeRegi;

  // create output R object list
  SEXP ReturnList;
  PROTECT(ReturnList = Rf_allocVector(VECSXP, 3));

  // name the objects in the list
  SEXP list_names;
  PROTECT(list_names = Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(list_names, 0, Rf_mkChar("FittedForest"));
  SET_STRING_ELT(list_names, 1, Rf_mkChar("ObsTrack"));
  SET_STRING_ELT(list_names, 2, Rf_mkChar("NodeRegi"));

  // move to list
  SET_VECTOR_ELT(ReturnList, 0, FittedForest);
  SET_VECTOR_ELT(ReturnList, 1, ObsTrackMat);
  SET_VECTOR_ELT(ReturnList, 2, NodeRegiMat);

  // set names
  Rf_setAttrib(ReturnList, R_NamesSymbol, list_names);

  // free memory
  //free(myPara);
  delete[] myPara;
  //free(X);
  delete[] X;
  delete[] subj_id;

  for (nt=0; nt<ntrees; nt++)
    delete[] VarImp[nt];

  delete[] VarImp;

  // unprotect
  UNPROTECT(7+ntrees);

  return ReturnList;
}


// print
// [[Rcpp::export()]]
void survForestPrint(SEXP parameters_R)
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
SEXP survForestPredict(SEXP testsetX_R,
                       SEXP FittedForest_R,
                       SEXP datasetY_R,
                       SEXP datasetCensor_R,
                       SEXP datasetNcat_R,
                       SEXP subjectweight_R,
                       SEXP ObsTrackMat_R,
                       SEXP NodeRegiMat_R,
                       SEXP parameters_R,
                       SEXP usecores_R)
{

  // copy parameters and check
  //PARAMETERS *myPara = (PARAMETERS*) malloc(sizeof(PARAMETERS));
  PARAMETERS* myPara = new PARAMETERS();
  copyParameters(myPara, parameters_R);

  //// create data objects
  int use_cores = INTEGER(usecores_R)[0];
  if (use_cores <= 0) use_cores = imax(1, omp_get_max_threads() - 1);

  int N = myPara->N;
  int P = myPara->P;
  int ntrees = myPara->ntrees;
  int Nfail = myPara->Nfail;

  SEXP dataX_dim = Rf_getAttrib(testsetX_R, R_DimSymbol);
  int testN = INTEGER(dataX_dim)[0];

  int i;
  int j;
  int nt = 0;

  // get X, Y, Censor

  //double **testX = (double **) malloc(P * sizeof(double *));
  double **testX = new double*[P];
  //if (testX == NULL) error("Unable to malloc testX");
  if (testX == NULL) stop("Unable to malloc testX");
  for (j = 0; j < P; j++)
    testX[j] = &REAL(testsetX_R)[j*testN];

  const int *Y = INTEGER(datasetY_R);
  const int *Censor = INTEGER(datasetCensor_R);
  const int *Ncat = INTEGER(datasetNcat_R);

  const double *subjectweight = REAL(subjectweight_R);

  // get tree matrix

  int TreeWidth = 4;
  int TreeLength;

  //double ***tree_matrix = (double ***) malloc(ntrees * sizeof(double **));
  double ***tree_matrix = new double**[ntrees];
  //if (tree_matrix == NULL) error("Unable to malloc for tree_matrix");
  if (tree_matrix == NULL) stop("Unable to malloc for tree_matrix");

  for (nt=0; nt<ntrees; nt++)
  {
    //tree_matrix[nt] = (double **) malloc(TreeWidth  * sizeof(double*));
    tree_matrix[nt] = new double*[TreeWidth];
    //if (tree_matrix[nt] == NULL) error("Unable to malloc for tree_matrix");
    if (tree_matrix[nt] == NULL) stop("Unable to malloc for tree_matrix");

    TreeLength = INTEGER(Rf_getAttrib(VECTOR_ELT(FittedForest_R, nt), R_DimSymbol))[0];

    for (i = 0; i < TreeWidth; i++)
      tree_matrix[nt][i] = &REAL(VECTOR_ELT(FittedForest_R, nt))[i*TreeLength];
  }

  // get ObsTrack

  //int **ObsTrack = (int **) malloc(ntrees * sizeof(int *));
  int **ObsTrack = new int*[ntrees];
  //if (ObsTrack == NULL) error("Unable to malloc ObsTrack");
  if (ObsTrack == NULL) stop("Unable to malloc ObsTrack");
  for (nt = 0; nt < ntrees; nt++)
    ObsTrack[nt] = &INTEGER(ObsTrackMat_R)[nt*N];

  // get NodeRegi
  //int **NodeRegi = (int **) malloc(ntrees * sizeof(int *));
  int **NodeRegi = new int*[ntrees];
  //if (NodeRegi == NULL) error("Unable to malloc NodeRegi");
  if (NodeRegi == NULL) stop("Unable to malloc NodeRegi");
  for (nt = 0; nt < ntrees; nt++)
    NodeRegi[nt] = &INTEGER(NodeRegiMat_R)[nt*N];


  SEXP SurvMat;
  PROTECT(SurvMat = Rf_allocMatrix(REALSXP, Nfail+1, testN));

  //double **surv_matrix = (double **) malloc(testN * sizeof(double *));
  double **surv_matrix = new double*[testN];
  //if (surv_matrix == NULL) error("Unable to malloc surv_matrix");
  if (surv_matrix == NULL) stop("Unable to malloc surv_matrix");
  for (i = 0; i < testN; i++)
  {
    surv_matrix[i] = &REAL(SurvMat)[i*(Nfail+1)];
    for(j=0; j<= Nfail; j++)
      surv_matrix[i][j] = 0;
  }

  PredictSurvivalKernel((const double**) testX,
                        Y,
                        Censor,
                        Ncat,
                        subjectweight,
                        (const double***) tree_matrix,
                        (const int**) ObsTrack,
                        (const int**) NodeRegi,
                        surv_matrix,
                        (const PARAMETERS*) myPara,
                        testN,
                        use_cores);


  //free(testX);
  delete[] testX;

  for (nt = 0; nt < ntrees; nt++)
    //free(tree_matrix[nt]);
    delete[] tree_matrix[nt];
  //free(tree_matrix);
  delete[] tree_matrix;

  //free(NodeRegi);
  delete[] NodeRegi;
  //free(ObsTrack);
  delete[] ObsTrack;
  //free(surv_matrix);
  delete[] surv_matrix;

  UNPROTECT(1);

  return SurvMat;
}


