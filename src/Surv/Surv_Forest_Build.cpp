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
                     std::vector< std::vector< ivec > > &NodeRegi,
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

  auto start = std::chrono::system_clock::now();
    
  standardize(variableweight, P);	// this cause precision loss...

  auto t1 = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = t1-start;
  //Rcout << "Time to run standardize: " << diff.count() << std::endl;
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
    auto t2 = std::chrono::system_clock::now();
    //R_DBP("Start tree %i\n",nt);
    int i;
    // in-bag and out-of-bag data indicator
    //int *inbagObs = (int *) malloc(size * sizeof(int));
    //int *inbagObs = new int[size];
    ivec inbagObs(size);
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
      ObsTrack(inbagObs[i],nt)++;

    //int * Ytemp = (int *) malloc(size * sizeof(int));
    //int * Ytemp = new int[size];
    ivec Ytemp(size);
    for (i = 0; i< size; i++) Ytemp[i] = Y[inbagObs[i]];

    qSort_iindex(Ytemp, 0, size-1, inbagObs);
    //inbagObs = sort_index(Ytemp);
    //Ytemp = sort(Ytemp);

    //TREENODE *TreeRoot = (TREENODE*) malloc(sizeof(TREENODE));
    TREENODE *TreeRoot = new TREENODE;

    // index vector will get destroyed within the tree fitting process.
    //int *inbagObs_copy = (int *) malloc(size * sizeof(int));
    //int *inbagObs_copy = new int[size];
    ivec inbagObs_copy(size);
    for (i=0; i < size ; i++) inbagObs_copy[i] = inbagObs[i];

    auto t3 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff2 = t3-t2;
    //Rcout << "Time to run pre-split: " << diff2.count() << std::endl;
    
    // start to build the tree
    Surv_Split_A_Node(TreeRoot, X, Y, Censor, Ncat, Interval, myPara, subjectweight, inbagObs_copy, size, variableweight, var_id, P);
    
    auto t4 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff3 = t4-t3;
    //Rcout << "Time to run split: " << diff3.count() << std::endl;
    // covert nodes to tree and NodeRegi.
    int Node = 0;

    Forest[nt] = TreeRoot;
    
    Record_NodeRegi(&Node, TreeRoot, NodeRegi, nt);//NodeRegi[nt]

    auto t5 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff4 = t5-t4;
    //Rcout << "Time to record tree: " << diff4.count() << std::endl;
    // summarize what observations are used in this tree

    //free(Ytemp);
    //delete[] Ytemp;
    //free(inbagObs);
    //delete[] inbagObs;
    //free(oobagObs);
    //delete[] oobagObs;

    // if (importance)
    // {
    //   //Rprintf("Do importance here \n");
    //   bool Usedj = FALSE;
    // 
    //   //for (int j = 0; j < P; j++)
    //   //  Usedj = CheckVar(TreeRoot, j);
    //    
    //    int nsim = 5;
    //    
    //    vec imp_nt(P);
    //    for(int j=0; j < P; j++){
    //      bool Usedj = FALSE;
    //      CheckVar(TreeRoot, j, Usedj);
    //      if(Usedj){
    //        std::vector<  colvec > Xperm = X;
    //        vec sim(nsim); //Running five simulations with different permutations
    //        for(int k=0; k < nsim; k++){
    //          vec Xp(N);
    //          permute(Xp, N);
    //          Xperm[j] = Xp;
    //      }
    //        
    //      }else{
    //        imp_nt[P] = -1;
    //      }
    //   }
    // 
    // 
    // }
    auto t2b = std::chrono::system_clock::now();
    std::chrono::duration<double> difft = t2b-t2;
    //Rcout << "Time to run tree: " << difft.count() << std::endl;

  }





  return;
}

vec martin_resid(const ivec Censor, const ivec Y, ivec obs, mat surv_matrix){
  
  int Nb = obs.size(); //Find the number of observations
  vec MResid(Nb); //Create a vector in which to hold the residuals
  
  for(int i=0; i < Nb; i++){ //For each observation included
    MResid[i] = Censor[obs[i]] - log(surv_matrix(obs[i],Y[i])); //Find the Martingale residual
  }
  
  return(MResid);
}

vec dev_resid(const ivec Censor, ivec obs, vec MResid){
  
  int Nb = obs.size(); //Find the number of observations
  vec DResid(Nb); //Create a vector in which to hold the residuals
  
  for(int i=0; i < Nb; i++){ //For each observation included
    DResid[i] = ((MResid[i] > 0) - (MResid[i] < 0))*sqrt(-2*(MResid[i]+Censor[obs[i]]*log(Censor[obs[i]]-MResid[i]))); //Find the Deviance residual
  }
  
  return(DResid);
}


void Record_NodeRegi(int* Node, TREENODE* TreeRoot, std::vector< std::vector< ivec > > &NodeRegi, int nt)
{
  *Node += 1;

  if (TreeRoot->Var == -1) // terminal node
  {
    //for (int i = 0; i< TreeRoot->NodeSize; i++){
      //NodeRegi(TreeRoot->NodeObs[i],nt) = *Node;
      //NodeRegi[nt].push_back(TreeRoot->NodeObs[i]);//Adds all observations in this terminal node
    //}
    NodeRegi[nt].push_back(TreeRoot->NodeObs);
  }else{
    
    ivec empty(1);
    empty.fill(-1);
    
    NodeRegi[nt].push_back(empty);//Not terminal node, so pushes empty vector

    Record_NodeRegi(Node, TreeRoot->Left, NodeRegi, nt);

    Record_NodeRegi(Node, TreeRoot->Right, NodeRegi, nt);
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


