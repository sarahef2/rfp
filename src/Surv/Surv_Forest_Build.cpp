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
//# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

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
                     int use_cores)
{
  int ntrees = myPara->ntrees;
  int verbose = myPara->verbose;
  int replacement = myPara->replacement;
  int importance = myPara->importance;
  double resample_prob = myPara->resample_prob;
  int nimpute = myPara->nimpute;
  int size = (int) N*resample_prob;
  int nt;

  // normalize the variable weight

  standardize(variableweight, P);	// this cause precision loss...

  // parallel computing... set cores

  use_cores = imax(1, use_cores);

  if (use_cores > 0) OMPMSG(1);

  int haveCores = omp_get_max_threads();

  if(use_cores > haveCores)
  {
    if (verbose) Rprintf("Do not have %i cores, use maximum %i cores. \n", use_cores, haveCores);
    use_cores = haveCores;
  }

#pragma omp parallel for schedule(static) num_threads(use_cores)
  for (nt = 0; nt < ntrees; nt++) // fit all trees
  {
    //R_DBP("Start tree %i\n",nt);
    int i;
    // in-bag and out-of-bag data indicator
    //int *inbagObs = (int *) malloc(size * sizeof(int));
    //int *inbagObs = new int[size];
    uvec inbagObs(size);
    ////int *oobagObs = (int *) malloc(N * sizeof(int)); // initiate a longer one
    //int *oobagObs = new int[N]; // initiate a longer one
    ivec oobagObs(N);
    oobagObs.fill(0);
    

    int OneSub;
    int oobag_n;

    for (i=0; i < N; i++)
      oobagObs[i] = subj_id[i];

    // sample in-bag and out-of-bag observations
    if (replacement)
    {

      for (i = 0; i < size; i++)
      {
        OneSub = random_in_range(0, N);
        inbagObs[i] = subj_id[OneSub];
        oobagObs[OneSub] = -1;
      }

      oobag_n = N;

      for (i=0; i<oobag_n; i++)
        if (oobagObs[i] < 0)
          oobagObs[i--] = oobagObs[--oobag_n];

    }else{
      for (i = 0; i < size; i++)
      {
        OneSub = random_in_range(0, N-i);
        inbagObs[i] = oobagObs[OneSub];
        oobagObs[OneSub] = oobagObs[N-1-i];
      }
      oobag_n = N - size;
    }

    // record the observations
    // if the ObsTrack is positive, then its in the fitting set, can have multiple counts
    // if the ObsTrack is zero, then its in the out-of-bag data

    for (i=0; i< size; i++)
      ObsTrack(nt,inbagObs[i])++;

    //int * Ytemp = (int *) malloc(size * sizeof(int));
    //int * Ytemp = new int[size];
    ivec Ytemp(size);
    for (i = 0; i< size; i++) Ytemp[i] = Y[inbagObs[i]];

    //qSort_iindex(Ytemp, 0, size-1, inbagObs);
    inbagObs = sort_index(Ytemp);
    Ytemp = sort(Ytemp);

    //TREENODE *TreeRoot = (TREENODE*) malloc(sizeof(TREENODE));
    TREENODE *TreeRoot = new TREENODE;

    // index vector will get destroyed within the tree fitting process.
    //int *inbagObs_copy = (int *) malloc(size * sizeof(int));
    //int *inbagObs_copy = new int[size];
    ivec inbagObs_copy(size);
    for (i=0; i < size ; i++) inbagObs_copy[i] = inbagObs[i];
    
    // start to build the
    Surv_Split_A_Node(TreeRoot, X, Y, Censor, Ncat, Interval, myPara, subjectweight, inbagObs_copy, size, variableweight, var_id, P);
    
    // covert nodes to tree and NodeRegi.
    int Node = 0;

    Forest[nt] = TreeRoot;
    
    Record_NodeRegi(&Node, TreeRoot, NodeRegi.colptr(nt));//NodeRegi[nt]

    // summarize what observations are used in this tree

    //free(Ytemp);
    //delete[] Ytemp;
    //free(inbagObs);
    //delete[] inbagObs;
    //free(oobagObs);
    //delete[] oobagObs;

    if (importance)
    {
      Rprintf("Do importance here \n");
      bool Usedj = FALSE;

      //for (int j = 0; j < P; j++)
      //  Usedj = CheckVar(TreeRoot, j);


    }

  }





  return;
}




void Record_NodeRegi(int* Node, TREENODE* TreeRoot, int* NodeRegi_nt)
{
  *Node += 1;

  if (TreeRoot->Var == -1) // terminal node
  {
    for (int i = 0; i< TreeRoot->NodeSize; i++){
      NodeRegi_nt[TreeRoot->NodeObs[i]] = *Node;
    }
  }else{

    Record_NodeRegi(Node, TreeRoot->Left, NodeRegi_nt);

    Record_NodeRegi(Node, TreeRoot->Right, NodeRegi_nt);
  }
}


void Record_Tree(int* Node, TREENODE* TreeRoot, mat &FittedTree, int TreeLength)
{
  *Node += 1;
  
  if (TreeRoot->Var == -1) // terminal node
  {
    //REAL(FittedTree)[*Node - 1] = -1;
    FittedTree(*Node-1,0)=-1;
    //REAL(FittedTree)[*Node - 1 + TreeLength] = NAN;
    FittedTree(*Node - 1 + TreeLength)=NAN;
    //REAL(FittedTree)[*Node - 1 + 2*TreeLength] = NAN;
    FittedTree(*Node - 1 + 2*TreeLength) = NAN;
    //REAL(FittedTree)[*Node - 1 + 3*TreeLength] = NAN;
    FittedTree(*Node - 1 + 3*TreeLength) = NAN;
    //free(TreeRoot->NodeObs);
  }else{

    //REAL(FittedTree)[*Node - 1] = TreeRoot->Var + 1;
    FittedTree[*Node - 1] = TreeRoot->Var + 1;
    //REAL(FittedTree)[*Node - 1 + TreeLength] = TreeRoot->Val;
    FittedTree[*Node - 1 + TreeLength] = TreeRoot->Val;

    int currentNode = *Node - 1;

    //REAL(FittedTree)[currentNode + 2*TreeLength] = *Node + 1;
    FittedTree(currentNode + 2*TreeLength) = *Node + 1;
    Record_Tree(Node, TreeRoot->Left, FittedTree, TreeLength);

    //REAL(FittedTree)[currentNode + 3*TreeLength] = *Node + 1;
    FittedTree(currentNode + 3*TreeLength) = *Node + 1;
    Record_Tree(Node, TreeRoot->Right, FittedTree, TreeLength);

  }
}


