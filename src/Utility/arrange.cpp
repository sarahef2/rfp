//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Utility Functions: arrange
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

# include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// my header file
# include "utility.h"

// random sample an integer value within [min, max-1]
int random_in_range(int min, int max)
{
  if (min == max)
    return min;

  double u;
  do {u = R::runif((double) min, (double) max);} while (u <= min || u >= max);

  return (int) u; 
  
  // generates integers from min to max-1  
  // Solution from https://stackoverflow.com/questions/5008804/generating-random-integer-from-a-range
  // std::default_random_engine generator;
  // std::uniform_int_distribution<int> uni(min, max-1); // guaranteed unbiased
  // int u = uni(generator);
  // Rcout << "uniform_int_distribution ("<< min << ","<<max<<"): " << u << std::endl;
}

// sample with weight from the tail of a vector and rotate the sampled one to the front

int sample_rotate(ivec &index, vec &weights, int start, int end)
{
  // sample without replacement, will put sampled one to the front
  // will order weights vector in the same way
  int a = random_in_range(start, end);
  int temp = index[a];
  double tempw = weights[a];

  index[a] = index[start];
  index[start] = temp;

  weights[a] = weights[start];
  weights[start] = tempw;

  return temp;
}

// random sample with weight

int weighted_sample(const double* x, int n)
{
  double a = R::runif(0, 1);

  int i;

  for (i = 0; i< n; i++)
  {
    a -= x[i];

    if (a <= 0)
      return i;
  }

  if (a < WeightTH)
    return weighted_sample(x, n);
  else
    ::Rf_error("weighted vector is not properly normalized");
}

// variable pack

double pack(const int nBits, const ivec bits) // from Andy's rf package
{
  int i;
  double value = bits[nBits - 1];

  for (i = nBits - 2; i >= 0; i--)
    value = 2.0*value + bits[i];

  return(value);
}

void unpack(const double pack, const int nBits, ivec bits) // from Andy's rf package
{
  int i;
  double x = pack;
  for (i = 0; i < nBits; ++i)
  {
    bits[i] = ((unsigned long) x & 1) ? 1 : 0;
    x /= 2;
  }
}

int unpack_goright(double pack, const int cat)
{
  int i;

  for (i = 0; i < cat; i++) pack /= 2;

  return(((unsigned long) pack & 1) ? 1 : 0);
}

/*
// test the pack
int nc = 10;
int* testpack = (int *) malloc(nc * sizeof(int));
double packed;

for (int k=0; k < 10; k++)
{
memset(testpack, 0, nc*sizeof(int));

for (int i=0; i<nc; i++)
testpack[i] = 1;

memset(testpack, 0, random_in_range(1, nc)*sizeof(int));
permute_i(testpack, nc);

for (int i=0; i<nc; i++)
R_DBP(" %i", testpack[i]);

R_DBP(" \n");

packed = pack(nc, testpack);

for (int i=0; i<nc; i++)
R_DBP(" %i", unpack_goright(packed, i));

R_DBP(" \n\n\n");

}
*/

