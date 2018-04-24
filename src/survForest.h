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
# include <R.h>

# include "utilities.h"


#ifndef survForest_Fun
#define survForest_Fun

SEXP survForestFit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void survForestBuild(//const double** X,
                     const std::vector<  colvec > X,
                     const ivec Y,
                     const ivec Censor,
                     const ivec Ncat,
                     const vec Interval,
                     const PARAMETERS* myPara,
                     const vec subjectweight,
                     const ivec subj_id,
                     const int N,
                     vec variableweight,
                     ivec var_id,
                     const int P,
                     TREENODE** Forest,
                     imat &ObsTrack,
                     imat &NodeRegi,
                     mat VarImp,
                     int use_cores);

void push_censor_front(int* inbagObs, int* Y, int* Censor, int size);

void Record_NodeRegi(int* Node, TREENODE* TreeRoot, imat& NodeRegi, int nt);

void Record_Tree(int* Node, TREENODE* TreeRoot, mat &FittedTree, int TreeLength);

void Surv_Split_A_Node(TREENODE* Node,
                  //const double** X,
                  const std::vector<  colvec > X,
                  const ivec Y,
                  const ivec Censor,
                  const ivec Ncat,
                  const vec Interval,
                  const PARAMETERS* myPara,
                  const vec subjectweight,
                  ivec useObs,
                  const int node_n,
                  vec variableweight,
                  ivec variableindex,
                  const int P);

void Surv_Find_A_Split(int* splitVar,
                  double* splitVal,
                  //const double** X,
                  std::vector<  colvec > X,
                  ivec Y,
                  const ivec Censor,
                  const ivec Ncat,
                  const vec Interval,
                  const PARAMETERS* myPara,
                  const vec subjectweight,
                  ivec useObs,
                  const int node_n,
                  vec variableweight,
                  ivec variableindex,
                  const int P);

void collapse(const ivec Y, const ivec Censor, ivec &Y_collapse, ivec &Censor_collapse, const ivec useObs, int node_n, int &nfail);



void Surv_One_Split_Cat_W(double* cut,
                           double* score,
                           const ivec useObs,
                           int node_n,
                           //const double* x,
                           const colvec x,
                           const ivec Y,
                           const ivec Censor,
                           const vec subjectweight,
                           int ncat,
                           int timepoints,
                           int split_gen,
                           int split_rule,
                           int nsplit,
                           int nmin,
                           int alpha);

void Surv_One_Split_Cat(double* cut,
                        double* score,
                        const ivec useObs,
                        int node_n,
                        //const double* x,
                        const vec x,
                        const ivec Y,
                        const ivec Censor,
                        int ncat,
                        int timepoints,
                        int split_gen,
                        int split_rule,
                        int nsplit,
                        int mincount);

void Surv_One_Split_Cont(double* cut,
                        double* score,
                        const ivec useObs,
                        int node_n,
                        //const double* x,
                        const colvec x,
                        const ivec Y,
                        const ivec Censor,
                        int timepoints,
                        int split_gen,
                        int split_rule,
                        int nsplit,
                        int mincount);

void Surv_One_Split_Cont_W(double* cut,
                        double* score,
                        const ivec useObs,
                        int node_n,
                        //const double* x,
                        const vec x,
                        const ivec Y,
                        const ivec Censor,
                        const vec subjectweight,
                        int timepoints,
                        int split_gen,
                        int split_rule,
                        int nsplit,
                        int nmin,
                        int alpha);

double logrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints);
double suplogrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints);

double logrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftN, double AllN, int timepoints);
double suplogrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftN, double AllN, int timepoints);


// prediction functions


List survForestPredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void PredictSurvivalKernel(const double** X,
                           const int* Y,
                           const int* Censor,
                           const int* Ncat,
                           const double* subjectweight,
                           const double*** tree_matrix,
                           const int** ObsTrack,
                           const int** NodeRegi,
                           double** surv_matrix,
                           const PARAMETERS* myPara,
                           int testN,
                           int use_cores);

void Get_Kernel_Weights(int subj,
                        const double** X,
                        const int* Ncat,
                        const double** tree_matrix_nt,
                        const int* ObsTrack_nt,
                        const int* NodeRegi_nt,
                        double* weights,
                        const int N);

void Get_Kernel_Weights_w(int subj,
                          const double** X,
                          const int* Ncat,
                          const double** tree_matrix_nt,
                          const int* ObsTrack_nt,
                          const int* NodeRegi_nt,
                          const double* subjectweight,
                          double* weights,
                          const int N);

int get_terminal(int node, int subj, const double ** X, const int* Ncat, const double** tree_matrix_nt);

#endif









