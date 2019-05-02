//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Utility Functions and Definitions
//  **********************************
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <Rmath.h>

using namespace Rcpp;
using namespace arma;

#define WeightTH 1e-10

#ifndef RLT_DEFINITION
#define RLT_DEFINITION

// survival categories
typedef struct SURVCAT{
  int cat;
  int f;
  int c;
  ivec flist;
  ivec clist;
} SURVCAT;

typedef struct SURVCAT_w{
  int cat;
  double f;
  double c;
  vec flist;
  vec clist;
} SURVCAT_w;

// define tree node
typedef struct TREENODE {
  int Var;
  double Val;
  int NodeSize;
  ivec NodeObs;
  struct TREENODE *Left, *Right;
} TREENODE;

// parameters structure
typedef struct PARAMETERS{
  int N;
  int P;
  int Nfail;
  int ntrees;
  int mtry;
  int split_gen;
  int split_rule;
  int nsplit;
  int nmin;
  int nmin_control;
  int nmin_failure;
  double alpha;
  int replacement;
  double resample_prob;
  int use_sub_weight;
  int use_var_weight;
  int importance;
  int nimpute;
  int verbose;
} PARAMETERS;

#endif

// ****************//
// Check functions //
// ****************//

// check OMP
#ifdef _OPENMP
#include <omp.h>
#define OMPMSG(...)
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define OMPMSG(...) Rprintf("Package is not compiled with OpenMP (omp.h).\n")
#endif

#ifndef RLT_CHECK
#define RLT_CHECK

void checkCores(int &use_cores, int verbose);
void standardize(vec &x, int n);

SEXP getListElement(SEXP, const char *);
void copyParameters(PARAMETERS*, SEXP);
void printParameters(PARAMETERS*);

int TreeSize(TREENODE *root);
void CheckVar(mat tree_matrix_nt, int j, bool& Check);

#endif

// ******************//
// Arrange functions //
// ******************//

#ifndef RLT_ARRANGE
#define RLT_ARRANGE

int random_in_range(int, int);
int sample_rotate(ivec &index, vec &weights, int start, int end);
int weighted_sample(const double* x, int n);


// minimum

template <class T> const T& min(const T& a, const T& b) {
	return (a<b)?a:b;
};

// maximum
template <class T> const T& max(const T& a, const T& b) {
	return (a<b)?b:a;
};

// swap the values at two location

template <class T> void swap(T* x, const int i, const int j)
{
	T temp;
	temp = x[i];
	x[i] = x[j];
	x[j] = temp;
};

// swap index

template <class T> void swap(T* a, T* b)
{
	T temp = *a;
	*a = *b;
	*b = temp;
};

// permute

template <class T> void permute(T &x, int n)
{
  int i;
  int j;
  decltype(x[0]) temp = x[0];

  for (i = 0; i<n-1; i++)
  {
    j = random_in_range(i, n);
    temp = x[i];
    x[i] = x[j];
    x[j] = temp;
  }
};

// cat variables pack

double pack(const int nBits, const ivec bits);
void unpack(const double pack, const int nBits, ivec bits);
int unpack_goright(double pack, const int cat);

#endif

// ****************//
// Debug functions //
// ****************//

#ifndef RLT_DEBUG
#define RLT_DEBUG

// debug

#define DEBUG

// this debug function will output results directly to R
#ifdef DEBUG
#define R_DBP(...) Rprintf(__VA_ARGS__)
#else
#define R_DBP(...)
#endif

// this debug function will output results to a .txt file
#ifdef DEBUG
#define DEBUGPRINT(mode, x, n1, n2) printLog(mode, x, n1, n2)
#else
#define DEBUGPRINT(mode, x, n1, n2)
#endif

void printLog(const char*, const char*, const int, const double);

#endif










