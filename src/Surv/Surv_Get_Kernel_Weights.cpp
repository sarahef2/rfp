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

//# include <Rdefines.h>
//# include <R.h>
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void Get_Kernel_Weights(int subj,
                        const std::vector< colvec > X,
                        const ivec Ncat,
                        const mat tree_matrix_nt,
                        const ivec ObsTrack_nt,
                        const std::vector< ivec > NodeRegi_nt,
                        vec& weights,
                        const int N)
{
  int node = get_terminal(0, subj, X, Ncat, tree_matrix_nt) + 1;

  //for (int i = 0; i < N; i++)
  //{
  //
  //  if (NodeRegi_nt[i] == node)
  //  {
  //
  //    weights[i] += ObsTrack_nt[i];
  //  }
  //}
  
  for (int i = 0; i < NodeRegi_nt[node].size(); i++){
    weights[NodeRegi_nt[node][i]]++;
  }

  return;
}


void Get_Kernel_Weights_w(int subj,
                          const std::vector< vec > X,
                          const ivec Ncat,
                          const mat tree_matrix_nt,
                          const ivec ObsTrack_nt,
                          const std::vector< ivec > NodeRegi_nt,
                          const vec subjectweight,
                          vec& weights,
                          const int N)
{
  int node = get_terminal(0, subj, X, Ncat, tree_matrix_nt); //+ 1;

//  for (int i = 0; i < N; i++)
//  {
//    if (NodeRegi_nt[i] == node)
//    {
//      weights[i] += ObsTrack_nt[i]*subjectweight[i];
//    }
//  }
  for (int i = 0; i < NodeRegi_nt[node].n_elem; i++){
    weights[NodeRegi_nt[node][i]]+=subjectweight[NodeRegi_nt[node][i]];
  }
  
  Rcout << "Finished" << std::endl;;
  
  return;
}


int get_terminal(int node, int subj, const std::vector< colvec > X, const ivec Ncat, const mat tree_matrix_nt)
{
  if (tree_matrix_nt(node,0) < 0)
    return node;

  int splitvar = (int) tree_matrix_nt(node,0) - 1;

  if (Ncat[splitvar] > 1)
  {
    if (unpack_goright(tree_matrix_nt(node,1), X[splitvar][subj] -1) == 0)
      return get_terminal((int) tree_matrix_nt(node,2) - 1, subj, X, Ncat, tree_matrix_nt);
    else
      return get_terminal((int) tree_matrix_nt(node,3) - 1, subj, X, Ncat, tree_matrix_nt);
  }else{
    if (X[splitvar][subj] <= tree_matrix_nt(node,1))
      return get_terminal((int) tree_matrix_nt(node,2) - 1, subj, X, Ncat, tree_matrix_nt);
    else
      return get_terminal((int) tree_matrix_nt(node,3) - 1, subj, X, Ncat, tree_matrix_nt);
  }
}
