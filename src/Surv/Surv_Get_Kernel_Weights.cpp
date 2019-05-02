//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Survival
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// my header file
# include "..//survForest.h"
# include "..//Utility//utility.h"

void Get_Kernel_Weights(int subj,
                        const std::vector< colvec > &X,
                        const ivec &Ncat,
                        const mat &tree_matrix_nt,
                        imat &ObsTerminal,
                        const std::vector< ivec > &NodeRegi_nt,
                        vec& weights,
                        const int &N,
                        bool &InTrainSet,
                        const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                        const ivec &perm_j, //The permuted variable (for variable importance)
                        int &subj_perm_loc, //The location of the subject's permuted value
                        int &nt //tree index
                          )
{
  int node;
  bool check;
  if(InTrainSet){
    check = ObsTerminal(subj,nt) < 0 or perm_ind > -1 or (not InTrainSet);
  }else{
    check = perm_ind > -1 or (not InTrainSet);
  }
  if(check){
    node = get_terminal(0, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);// + 1;
    if(InTrainSet and perm_ind < 0){
      ObsTerminal(subj,nt) = node;
    }
  }else{
    node =  ObsTerminal(subj,nt);
  }

  for (int i = 0; i < NodeRegi_nt[node].n_elem; i++){
    weights[NodeRegi_nt[node][i]]++;
  }

  return;
}


void Get_Kernel_Weights_w(int subj,
                          const std::vector< vec > &X,
                          const ivec &Ncat,
                          const mat &tree_matrix_nt,
                          imat &ObsTerminal,
                          const std::vector< ivec > &NodeRegi_nt,
                          const vec &subjectweight,
                          vec& weights,
                          const int &N,
                          bool &InTrainSet,
                          const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                          const ivec &perm_j, //The permuted variable (for variable importance)
                          int &subj_perm_loc, //The location of the subject's permuted value
                          int &nt
                            )
{

  int node;
  if(ObsTerminal(subj,nt) < 0 or (not InTrainSet) or perm_ind > -1){
    node = get_terminal(0, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);// + 1;
    if(InTrainSet and perm_ind < 0){
      ObsTerminal(subj,nt) = node;
    }
  }else{
    node =  ObsTerminal(subj,nt) - 1;
  }

  for (int i = 0; i < NodeRegi_nt[node].n_elem; i++){
    weights[NodeRegi_nt[node][i]]+=subjectweight[NodeRegi_nt[node][i]];
  }

  return;
}


int get_terminal(int node, int &subj, const std::vector< colvec > &X, const ivec &Ncat, const mat &tree_matrix_nt,
                 const int &perm_ind, //The variable index, if any, permuted (for variable importance)
                 const ivec &perm_j, //The permuted variable (for variable importance)
                 int &subj_perm_loc //The location of the subject's permuted value
                   )
{
  if (tree_matrix_nt(node,0) < 0)
    return node;

  int splitvar = (int) tree_matrix_nt(node,0) - 1;

  double subj_val;

  if(splitvar == perm_ind){
    subj_val = X[splitvar][perm_j[subj_perm_loc]];
  }else{
    subj_val = X[splitvar][subj];
  }

  if (Ncat[splitvar] > 1)
  {
    if (unpack_goright(tree_matrix_nt(node,1), subj_val - 1) == 0)
      return get_terminal((int) tree_matrix_nt(node,2) - 1, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);
    else
      return get_terminal((int) tree_matrix_nt(node,3) - 1, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);
  }else{
    if (subj_val <= tree_matrix_nt(node,1))
      return get_terminal((int) tree_matrix_nt(node,2) - 1, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);
    else
      return get_terminal((int) tree_matrix_nt(node,3) - 1, subj, X, Ncat, tree_matrix_nt, perm_ind, perm_j, subj_perm_loc);
  }
}
