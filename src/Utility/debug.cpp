//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Utility Functions: debug
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// my header file
# include "utility.h"

// debug function

void printLog(const char* mode, const char* x, const int n1, const double n2)
{
  FILE* pFile = fopen("RLT_Debug_log.txt", mode);

  if(pFile != NULL)
    fprintf(pFile, x, n1, n2);

  fclose(pFile);
  return;
}
