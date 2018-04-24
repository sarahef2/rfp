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
                        const double** X,
                        const int* Ncat,
                        const double** tree_matrix_nt,
                        const int* ObsTrack_nt,
                        const int* NodeRegi_nt,
                        double* weights,
                        const int N)
{
  int node = get_terminal(0, subj, X, Ncat, tree_matrix_nt) + 1;

  for (int i = 0; i < N; i++)
  {

    if (NodeRegi_nt[i] == node)
    {

      weights[i] += ObsTrack_nt[i];
    }
  }

  return;
}


void Get_Kernel_Weights_w(int subj,
                          const double** X,
                          const int* Ncat,
                          const double** tree_matrix_nt,
                          const int* ObsTrack_nt,
                          const int* NodeRegi_nt,
                          const double* subjectweight,
                          double* weights,
                          const int N)
{
  int node = get_terminal(0, subj, X, Ncat, tree_matrix_nt) + 1;

  for (int i = 0; i < N; i++)
  {
    if (NodeRegi_nt[i] == node)
    {
      weights[i] += ObsTrack_nt[i]*subjectweight[i];
    }
  }

  return;
}


int get_terminal(int node, int subj, const double ** X, const int* Ncat, const double** tree_matrix_nt)
{
  if (tree_matrix_nt[0][node] < 0)
    return node;

  int splitvar = (int) tree_matrix_nt[0][node] - 1;

  if (Ncat[splitvar] > 1)
  {
    if (unpack_goright(tree_matrix_nt[1][node], X[splitvar][subj] -1) == 0)
      return get_terminal((int) tree_matrix_nt[2][node] - 1, subj, X, Ncat, tree_matrix_nt);
    else
      return get_terminal((int) tree_matrix_nt[3][node] - 1, subj, X, Ncat, tree_matrix_nt);
  }else{
    if (X[splitvar][subj] <= tree_matrix_nt[1][node])
      return get_terminal((int) tree_matrix_nt[2][node] - 1, subj, X, Ncat, tree_matrix_nt);
    else
      return get_terminal((int) tree_matrix_nt[3][node] - 1, subj, X, Ncat, tree_matrix_nt);
  }
}
