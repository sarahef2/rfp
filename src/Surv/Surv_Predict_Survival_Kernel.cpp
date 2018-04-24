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
//# include <R.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void PredictSurvivalKernel2(const double** X,
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
                           int use_cores)
{
  int N = myPara->N;
  int Nfail = myPara->Nfail;
  int ntrees = myPara->ntrees;
  int verbose = myPara->verbose;
  int use_sub_weight = myPara->use_sub_weight;
  int i;

  // parallel computing... set cores

  use_cores = imax(1, use_cores);

  int haveCores = omp_get_max_threads();

  if(use_cores > haveCores)
  {
    if (verbose) Rprintf("Do not have %i cores, use maximum %i cores. \n", use_cores, haveCores);
    use_cores = haveCores;
  }

  // matrix keep track of the sum of fail and censor at each time point
  //double **remove_matrix = (double **) malloc(testN * sizeof(double *));
  double **remove_matrix =new double*[testN];
  //if (remove_matrix == NULL) error("Unable to malloc remove_matrix");
  if (remove_matrix == NULL) stop("Unable to malloc remove_matrix");
  for (i = 0; i < testN; i++)
    //remove_matrix[i] = (double *) calloc(Nfail+1, sizeof(double));
    remove_matrix[i] = new double[Nfail+1];

#pragma omp parallel for schedule(guided) num_threads(use_cores)
  for (i = 0; i < testN; i++)
  {
    //double* weights = (double *) calloc(N, sizeof(double));
    double* weights = new double[N];
    int j, nt;

    if (use_sub_weight)
      for (nt = 0; nt < ntrees; nt++)
        Get_Kernel_Weights_w(i, X, Ncat, tree_matrix[nt], ObsTrack[nt], NodeRegi[nt], subjectweight, weights, N);
    else
      for (nt = 0; nt < ntrees; nt++)
        Get_Kernel_Weights(i, X, Ncat, tree_matrix[nt], ObsTrack[nt], NodeRegi[nt], weights, N);

    double weights_sum = 0;

    for (j = 0; j < N; j++)
    {
      remove_matrix[i][Y[j]] += weights[j];

      if (Censor[j] == 1)
        surv_matrix[i][Y[j]] += weights[j];

      weights_sum += weights[j];
    }

    surv_matrix[i][0] = 1;
    weights_sum -= remove_matrix[i][0];

    // KM survival function
    for (j = 1; j <= Nfail; j++)
    {
      surv_matrix[i][j] = surv_matrix[i][j-1] * (1 - surv_matrix[i][j]/weights_sum);
      weights_sum -= remove_matrix[i][j];
    }

    //free(weights);
    delete[] weights;
  }


return;
}

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
                           int use_cores){
  
  int N = myPara->N;
  int Nfail = myPara->Nfail;
  int ntrees = myPara->ntrees;
  int verbose = myPara->verbose;
  int use_sub_weight = myPara->use_sub_weight;
  int i;
  
  // parallel computing... set cores
  
  use_cores = imax(1, use_cores);
  
  int haveCores = omp_get_max_threads();
  
  if(use_cores > haveCores)
  {
    if (verbose) Rprintf("Do not have %i cores, use maximum %i cores. \n", use_cores, haveCores);
    use_cores = haveCores;
  }
  
  // matrix keep track of the sum of fail and censor at each time point
  double **remove_matrix = (double **) malloc(testN * sizeof(double *));
  //if (remove_matrix == NULL) error("Unable to malloc remove_matrix");
  for (i = 0; i < testN; i++)
    remove_matrix[i] = (double *) calloc(Nfail+1, sizeof(double));
  
#pragma omp parallel for schedule(guided) num_threads(use_cores)
  for (i = 0; i < testN; i++)
  {
    double* weights = (double *) calloc(N, sizeof(double));
    int j, nt;
    
    if (use_sub_weight)
      for (nt = 0; nt < ntrees; nt++)
        Get_Kernel_Weights_w(i, X, Ncat, tree_matrix[nt], ObsTrack[nt], NodeRegi[nt], subjectweight, weights, N);
    else
      for (nt = 0; nt < ntrees; nt++)
        Get_Kernel_Weights(i, X, Ncat, tree_matrix[nt], ObsTrack[nt], NodeRegi[nt], weights, N);
    
    double weights_sum = 0;
    
    for (j = 0; j < N; j++)
    {
      remove_matrix[i][Y[j]] += weights[j];
      
      if (Censor[j] == 1)
        surv_matrix[i][Y[j]] += weights[j];
      
      weights_sum += weights[j];
    }
    
    surv_matrix[i][0] = 1;
    weights_sum -= remove_matrix[i][0];
    
    // KM survival function
    for (j = 1; j <= Nfail; j++)
    {
      surv_matrix[i][j] = surv_matrix[i][j-1] * (1 - surv_matrix[i][j]/weights_sum);
      weights_sum -= remove_matrix[i][j];
    }
    
    free(weights);
  }


return;

}
  