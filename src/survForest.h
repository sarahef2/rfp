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
# include <R.h>
using namespace Rcpp;

# include "utilities.h"


#ifndef survForest_Fun
#define survForest_Fun

SEXP survForestFit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

void survForestBuild(const std::vector<  colvec > &X,
                     const ivec &Y,
                     const ivec &Censor,
                     const ivec &Ncat,
                     const vec &Interval,
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
                     int &use_cores,
                     mat &oob_surv_matrix,
                     vec &oob_residuals,
                     int &counter);

void push_censor_front(int* inbagObs, int* Y, int* Censor, int size);

void Record_NodeRegi(int* Node, TREENODE* TreeRoot, std::vector< std::vector< ivec > > &NodeRegi, int nt, imat &ObsTerminal);

void Record_Tree(int* Node, TREENODE* TreeRoot, mat &FittedTree, int TreeLength);

void Surv_Split_A_Node(TREENODE* Node,
                  const std::vector<  colvec > X,
                  const ivec &Y,
                  const ivec &Censor,
                  const ivec &Ncat,
                  const vec &Interval,
                  const PARAMETERS* myPara,
                  const vec &subjectweight,
                  ivec &useObs,
                  const int &node_n,
                  vec &variableweight,
                  ivec &variableindex,
                  const int &P,
                  int &counter);

void Surv_Find_A_Split(int* splitVar,
                  double* splitVal,
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
                  const int P,
                  int &counter);

void collapse(const ivec Y, const ivec Censor, ivec &Y_collapse, ivec &Censor_collapse, const ivec useObs, int node_n, int &nfail);



void Surv_One_Split_Cat_W(double* cut,
                           double* score,
                           const ivec useObs,
                           int node_n,
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
                        const ivec &useObs,
                        int node_n,
                        const colvec &x,
                        const ivec &Y,
                        const ivec &Censor,
                        int &timepoints,
                        int &split_gen,
                        int &split_rule,
                        int &nsplit,
                        int &mincount);

void Surv_One_Split_Cont_W(double* cut,
                        double* score,
                        const ivec &useObs,
                        int node_n,
                        const vec &x,
                        const ivec &Y,
                        const ivec &Censor,
                        const vec &subjectweight,
                        int &timepoints,
                        int &split_gen,
                        int &split_rule,
                        int &nsplit,
                        int &nmin,
                        int &alpha);

double logrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints);
double suplogrank(ivec Left_Count_Fail, ivec Left_Count_Censor, ivec Right_Count_Fail, ivec Right_Count_Censor, double LeftN, double AllN, int timepoints);

double logrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftN, double AllN, int timepoints);
double suplogrank_w(vec Left_Count_Fail, vec Left_Count_Censor, vec Right_Count_Fail, vec Right_Count_Censor, double LeftN, double AllN, int timepoints);


// prediction functions


mat survForestPredict(mat, List, imat, ivec, ivec, vec, imat, List, List, int);

void PredictSurvivalKernel(const std::vector< colvec > &X,
                           const ivec &Y,
                           const ivec &Censor,
                           const ivec &Ncat,
                           const vec &subjectweight,
                           const std::vector< mat > &tree_matrix,
                           const imat &ObsTrack,
                           imat &ObsTerminal, //Save the terminal node number of the observations in each tree
                           const std::vector< std::vector< ivec > > &NodeRegi,
                           mat &surv_matrix,
                           const PARAMETERS* myPara,
                           int testN,
                           const ivec &use_obs, //The index of the obserations to predict
                           const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                           const ivec &perm_j, //The permuted variable (for variable importance)
                           int oob_only, //Should prediction be run on only oob obs
                           bool InTrainSet,
                           int use_cores);

void Get_Kernel_Weights(int subj,
                        const std::vector< vec > &X,
                        const ivec &Ncat,
                        const mat &tree_matrix_nt,
                        imat &ObsTerminal,
                        const std::vector< ivec > &NodeRegi_nt,
                        vec &weights,
                        const int &N,
                        bool &InTrainSet,
                        const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                        const ivec &perm_j, //The permuted variable (for variable importance)
                        int &subj_perm_loc,
                        int &nt
                          );

void Get_Kernel_Weights_w(int subj,
                          const std::vector< vec > &X,
                          const ivec &Ncat,
                          const mat &tree_matrix_nt,
                          imat &ObsTerminal,
                          const std::vector< ivec > &NodeRegi_nt,
                          const vec &subjectweight,
                          vec &weights,
                          const int &N,
                          bool &InTrainSet,
                          const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                          const ivec &perm_j, //The permuted variable (for variable importance)
                          int &subj_perm_loc,
                          int &nt
                            );

int get_terminal(int node, int &subj, const std::vector< vec > &X, const ivec &Ncat, const mat &tree_matrix_nt, const int &perm_ind, const ivec &perm_j, int &subj_perm_loc);

void martin_resid(const ivec &Censor, const ivec &Y, const ivec &obs, const int &Nb, mat &surv_matrix, vec &MResids);

void dev_resid(const ivec &Censor, const ivec &Y, const ivec &obs, const int &Nb, mat &surv_matrix, vec &DResid);
#endif