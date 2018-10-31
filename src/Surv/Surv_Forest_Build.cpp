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
# include <math.h>       /* isnan, sqrt */
// [[Rcpp::depends(RcppArmadillo)]]
//# include <Rcpp.h>
//# include <Rdefines.h>
//# include <R.h>
//# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void survForestBuild(const std::vector<  colvec > X,
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
                     imat &ObsTerminal,
                     std::vector< std::vector< ivec > > &NodeRegi,
                     vec &VarImp,
                     int use_cores,
                     mat &oob_surv_matrix,
                     vec &oob_residuals,
                     int &counter)
{
  int ntrees = myPara->ntrees;
  int verbose = myPara->verbose;
  int replacement = myPara->replacement;
  int importance = myPara->importance;
  double resample_prob = myPara->resample_prob;
  int nimpute = myPara->nimpute;
  int size = (int) N*resample_prob;
  int nt;
  int Nfail = myPara->Nfail; //Used in variable importance calculations
  std::vector< mat > tree_matrix(ntrees);
  mat VarImp_store(ntrees, P);
  VarImp_store.fill(0);

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
    //auto t2 = std::chrono::system_clock::now();
    //R_DBP("Start tree %i\n",nt);
    int i;
    // in-bag and out-of-bag data indicator
    ivec inbagObs(size);
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
      
      i=0;
      while(i<oobagObs.size())
        if(oobagObs(i) < 0)
          oobagObs.shed_row(i);
        else
          i++;

    }else{
      for (i = 0; i < size; i++)
      {
        OneSub = random_in_range(0, N-i);
        inbagObs[i] = oobagObs[OneSub];
        oobagObs[OneSub] = oobagObs[N-1-i];
      }
      i=1;
      while(i<oobagObs.size())
        if(oobagObs(i) == oobagObs(i-1))
          oobagObs.shed_row(i);
        else
          i++;
      oobag_n = N - size;
    }
    
    // record the observations
    // if the ObsTrack is positive, then its in the fitting set, can have multiple counts
    // if the ObsTrack is zero, then its in the out-of-bag data

    for (i=0; i< size; i++)
      ObsTrack(inbagObs[i],nt)++;

    ivec Ytemp(size);
    for (i = 0; i< size; i++) Ytemp[i] = Y[inbagObs[i]];

    qSort_iindex(Ytemp, 0, size-1, inbagObs);

    TREENODE *TreeRoot = new TREENODE;

    // index vector will get destroyed within the tree fitting process.
    ivec inbagObs_copy(size);
    for (i=0; i < size ; i++) inbagObs_copy[i] = inbagObs[i];

    // start to build the tree
    Surv_Split_A_Node(TreeRoot, X, Y, Censor, Ncat, Interval, myPara, subjectweight, inbagObs_copy, size, variableweight, var_id, P, counter);
    
    // covert nodes to tree and NodeRegi.
    int Node = 0;

    Forest[nt] = TreeRoot;
    
    Record_NodeRegi(&Node, TreeRoot, NodeRegi, nt, ObsTerminal);//NodeRegi[nt]

    mat FittedTree;
    int TreeLength;
    TreeLength = TreeSize(Forest[nt]);
    Node = 0;
    
    FittedTree.set_size(TreeLength,4);
    
    //Convert tree structure nodes into matrix
    Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);
    tree_matrix[nt] = FittedTree;
    // summarize what observations are used in this tree

  }

  if(importance){ //New version of importance
    ivec tmp(N);
    tmp.fill(0);
    PredictSurvivalKernel((const std::vector< colvec >) X,
                          Y,
                          Censor,
                          Ncat,
                          subjectweight,
                          (const std::vector< mat >) tree_matrix,
                          (const imat) ObsTrack,
                          ObsTerminal,
                          (const std::vector< std::vector< ivec > >) NodeRegi,
                          oob_surv_matrix,
                          (const PARAMETERS*) myPara,
                          N,
                          (const ivec) subj_id,
                          (const int) -1,
                          (const ivec) tmp,
                          1,
                          true,
                          use_cores);
    dev_resid(Censor, Y, subj_id, N, oob_surv_matrix, oob_residuals);
    double Dev_MSE = 0;
    for(int d=0; d < N; d++){
      Dev_MSE += oob_residuals[d]*oob_residuals[d];
    }

    //May want to make this variable?
    int nsim = 1;
#pragma omp parallel for schedule(static) num_threads(use_cores)
    for(int j=0; j < P; j++){
      vec imp_sim(nsim);
      imp_sim.fill(0);
      for(int k=0; k < nsim; k++){
          ivec perm_j = subj_id;//oobagObs;
          permute_i(perm_j, N);

          mat surv_matrix_perm(N,Nfail+1);
          surv_matrix_perm.fill(0);
          
          PredictSurvivalKernel((const std::vector< colvec >) X,
                                Y,
                                Censor,
                                Ncat,
                                subjectweight,
                                (const std::vector< mat >) tree_matrix, //tree_matrix, //Figure out
                                (const imat) ObsTrack,
                                ObsTerminal,
                                (const std::vector< std::vector< ivec > >) NodeRegi, //NodeRegi, //Figure out
                                surv_matrix_perm,
                                (const PARAMETERS*) myPara,
                                N,
                                (const ivec) subj_id,
                                (const int) j,
                                (const ivec) perm_j,
                                1,
                                true,
                                use_cores);
          vec resids_perm(N);
          resids_perm.fill(0);
          dev_resid(Censor, Y, subj_id, N, surv_matrix_perm, resids_perm);
          double Dev_MSE_perm = 0;
          for(int d=0; d < N; d++){
            Dev_MSE_perm += resids_perm[d]*resids_perm[d];
          }
          imp_sim[k] = Dev_MSE_perm;
        }
        for(int m = 0; m < nsim; m++){
          VarImp[j] += imp_sim[m] / nsim;
        }
        VarImp[j] = VarImp[j] / Dev_MSE - 1;
    }
    
    
  }
  
  return;
}

void martin_resid(const ivec Censor, const ivec Y, const ivec obs, int Nb, mat surv_matrix, vec &MResid){
  

  for(int i=0; i < Nb; i++){ //For each observation included
    MResid[i] = Censor[obs[i]] - (-log(surv_matrix(i,Y[obs[i]]))); //Find the Martingale residual
  }
}

void dev_resid(const ivec Censor, const ivec Y, const ivec obs, int Nb, mat surv_matrix, vec &DResid){
  vec MResid(Nb);
  MResid.fill(0);
  martin_resid(Censor, Y, obs, Nb, surv_matrix, MResid);

  double minMR = 1000;
  for(int i=0; i < Nb; i++){ //For each observation included
    if(MResid[i] < -100000000){
      if(minMR == 1000){
        for(int j=0; j < Nb; j++){
          if(MResid[j] < minMR &&  MResid[j] > -100000000){
            minMR = MResid[j];
          }
        }
      }
      //Rcout << "-infty res replace with " << minMR << std::endl;;
      MResid[i] = minMR;
    } 
    DResid[i] = ((MResid[i] > 0) - (MResid[i] < 0))*sqrt(-2*(MResid[i]+Censor[obs[i]]*log(Censor[obs[i]]-MResid[i]))); //Find the Deviance residual
  }
  
  double maxDR = -1000;
  double minDR = 1000;
  for(int i=0; i < Nb; i++){
    if(DResid[i] < -100000000 || DResid[i] != DResid[i]) {
      if(minDR == 1000){
        for(int j=0; j < Nb; j++){
          if(DResid[j] < minDR &&  DResid[j] > -100000000){
            minDR = DResid[j];
          }
        }
      }
      DResid[i] = minDR;
    }
    if(DResid[i] > 100000000) {
      if(maxDR == -1000){
        for(int j=0; j < Nb; j++){
          if(DResid[j] > maxDR &&  DResid[j] < 100000000){
            maxDR = DResid[j];
          }
        }
      }
      //Rcout << "+infty res replace with " << maxDR << std::endl;;
      DResid[i] = maxDR;
    }
  }

}


void Record_NodeRegi(int* Node, TREENODE* TreeRoot, std::vector< std::vector< ivec > > &NodeRegi, int nt, imat &ObsTerminal)
{
  *Node += 1;

  if (TreeRoot->Var == -1) // terminal node
  {
    NodeRegi[nt].push_back(TreeRoot->NodeObs);
    ivec inobs = TreeRoot->NodeObs;
    for(int i=0; i < inobs.size(); i++){
      int node_num = NodeRegi[nt].size();
      ObsTerminal(inobs(i),nt) = node_num-1;
    }
  }else{
    
    ivec empty(1);
    empty.fill(-1);
    
    NodeRegi[nt].push_back(empty);//Not terminal node, so pushes empty vector

    Record_NodeRegi(Node, TreeRoot->Left, NodeRegi, nt, ObsTerminal);

    Record_NodeRegi(Node, TreeRoot->Right, NodeRegi, nt, ObsTerminal);
  }
}


void Record_Tree(int* Node, TREENODE* TreeRoot, mat &FittedTree, int TreeLength)
{
  *Node += 1;
  
  if (TreeRoot->Var == -1) // terminal node
  {
    FittedTree(*Node-1,0)=-1;
    FittedTree(*Node - 1 + TreeLength)=NAN;
    FittedTree(*Node - 1 + 2*TreeLength) = NAN;
    FittedTree(*Node - 1 + 3*TreeLength) = NAN;
  }else{

    FittedTree[*Node - 1] = TreeRoot->Var + 1;
    FittedTree[*Node - 1 + TreeLength] = TreeRoot->Val;

    int currentNode = *Node - 1;

    FittedTree(currentNode + 2*TreeLength) = *Node + 1;
    Record_Tree(Node, TreeRoot->Left, FittedTree, TreeLength);

    FittedTree(currentNode + 3*TreeLength) = *Node + 1;
    Record_Tree(Node, TreeRoot->Right, FittedTree, TreeLength);

  }
}


