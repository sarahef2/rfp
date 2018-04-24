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
# include <stdbool.h>
//# include <Rdefines.h>
//# include <R.h>
# include <Rmath.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#define OMPMSG(...)
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define OMPMSG(...) Rprintf("Package is not compiled with OpenMP (omp.h).\n")
#endif

#define WeightTH 1e-10

#ifndef survForest_utility
#define survForest_utility

// debug

#define DEBUG

#ifdef DEBUG
#define R_DBP(...) Rprintf(__VA_ARGS__)
#else
#define R_DBP(...)
#endif

#ifdef DEBUG
#define DEBUGPRINT(mode, x, n1, n2) printLog(mode, x, n1, n2)
#else
#define DEBUGPRINT(mode, x, n1, n2)
#endif

void printLog(const char*, const char*, const int, const double);

typedef std::vector<int> stdvecint;
typedef std::vector< std::vector<int> > stdvecvecint;
typedef std::vector<double> stdvec;
typedef std::vector< std::vector<double> > stdvecvec;

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
  double alpha;
  int replacement;
  double resample_prob;
  int use_sub_weight;
  int use_var_weight;
  //int honest;
  int importance;
  int nimpute;
  //int use_cores;
  int verbose;
} PARAMETERS;


SEXP getListElement(SEXP, const char *);
void copyParameters(PARAMETERS*, SEXP);
void printParameters(PARAMETERS*);


// survival categories
typedef struct SURVCAT{
  int cat;
  int f;
  int c;
  ivec flist;
  ivec clist;
} SURVCAT;

void swap_SURVCAT(SURVCAT* a, SURVCAT* b);

typedef struct SURVCAT_w{
  int cat;
  double f;
  double c;
  vec flist;
  vec clist;
} SURVCAT_w;

void swap_SURVCAT_w(SURVCAT_w* a, SURVCAT_w* b);

void swap_SURVCAT_w(SURVCAT_w* a, SURVCAT_w* b);


// define tree node

typedef struct TREENODE {
  int Var;
  double Val;
  int NodeSize;
  ivec NodeObs;//int* NodeObs;
  struct TREENODE *Left, *Right;
} TREENODE;


// tree functions
int TreeSize(TREENODE *root);



// other utility functions

void standardize(//double*
    vec&, int);
int imin(int, int);
int imax(int, int);
double dmin(double, double);
double dmax(double, double);

void swap_d(double*, int i, int j);
void swap_i(int*, int i, int j);


// random number

int random_in_range(int, int);
int sample_rotate(ivec index, vec weights, int start, int end);
int weighted_sample(const double* x, int n);

// sorting

int partition_d(vec& keys,
                int low,
                int high,
                ivec& index);

void qSort_dindex(vec& keys,
                 int low,
                 int high,
                 ivec& index);

void iSort_index(vec& keys,
                 int low,
                 int high,
                 ivec& index);

int partition_i(ivec& keys,
                int low,
                int high,
                ivec& index);

void qSort_iindex(ivec& keys,
                  int low,
                  int high,
                  ivec& index);


void permute_i(ivec &x, int n);

// cat variables pack

double pack(const int nBits, const ivec bits);
void unpack(const double pack, const int nBits, ivec bits);
int unpack_goright(double pack, const int cat);


#endif
