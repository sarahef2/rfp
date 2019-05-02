//  **********************************
//  Reinforcement Learning Trees (RLT)
//  Survival
//  **********************************

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// my header file
# include "..//survForest.h"
# include "..//Utility//utility.h"

void survForestBuild(const std::vector<  colvec > &X,
                     const ivec &Y,
                     const ivec &Censor,
                     const ivec &Ncat,
                     const vec &Interval,
                     const PARAMETERS* myPara,
                     const vec &subjectweight,
                     const ivec &subj_id,
                     const int &N,
                     vec &variableweight,
                     ivec &var_id,
                     const int &P,
                     TREENODE** Forest,
                     imat &ObsTrack,
                     imat &ObsTerminal,
                     std::vector< std::vector< ivec > > &NodeRegi,
                     vec &VarImp,
                     int &use_cores,
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
  std::vector< mat > tree_matrix(ntrees);
  mat VarImp_store(ntrees, P);
  VarImp_store.fill(0);

  // normalize the variable weight

  standardize(variableweight, P);	// this may cause precision loss...

  // check parallel computing cores

  checkCores(use_cores, verbose);

  // start trees

#pragma omp parallel for schedule(static) num_threads(use_cores)
  for (nt = 0; nt < ntrees; nt++) // fit all trees
  {
    //R_DBP("Start tree %i\n",nt);
    int i;
    // in-bag and out-of-bag data indicator
    ivec inbagObs(size);
    ivec oobagObs(N);
    oobagObs.fill(0);

    int OneSub;
    int oobag_n;

    //R_DBP("Bootstrap Sample");
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
    ivec index(size);
    for (i = 0; i< size; i++) {
      Ytemp[i] = Y[inbagObs[i]];
      index[i] = i;
    }

    //Rcout << "Ytemp= " << Ytemp << std::endl;;
    //R_DBP("Sort Y\n");
    //Rcout << "Y= " << Y << std::endl;;
    //Rcout << "useObs= " << inbagObs << std::endl;;
    std::sort(index.begin(), index.end(), [&Ytemp](size_t i, size_t j) {return (Ytemp[i] < Ytemp[j]);});

    TREENODE *TreeRoot = new TREENODE;

    // index vector will get destroyed within the tree fitting process.
    ivec inbagObs_copy(size);
    for (i=0; i < size ; i++) {
      inbagObs_copy[i] = inbagObs[index[i]];
      }
    
    for(i=0; i < size; i++) inbagObs[i] = inbagObs_copy[i];

    //Rcout << "useObs= " << inbagObs << std::endl;;
    //R_DBP("Split %i\n",nt);
    // start to build the tree
    Surv_Split_A_Node(TreeRoot, X, Y, Censor, Ncat, Interval, myPara, subjectweight, inbagObs_copy, size, variableweight, var_id, P, counter);

    // covert nodes to tree and NodeRegi.
    int Node = 0;

    //R_DBP("Forest[nt]\n");
    Forest[nt] = TreeRoot;

    //R_DBP("Record_NodeRegi\n");
    Record_NodeRegi(&Node, TreeRoot, NodeRegi, nt, ObsTerminal);//NodeRegi[nt]

    mat FittedTree;
    int TreeLength;
    TreeLength = TreeSize(Forest[nt]);
    Node = 0;

    FittedTree.set_size(TreeLength,4);

    //Convert tree structure nodes into matrix
    //R_DBP("Record_Tree\n");
    Record_Tree(&Node, Forest[nt], FittedTree, TreeLength);
    tree_matrix[nt] = FittedTree;
    // summarize what observations are used in this tree

  }

  if(importance){
    Variable_Importance(X, Y, Censor, Ncat, myPara, subjectweight, tree_matrix, subj_id, N, P, ObsTrack, ObsTerminal, NodeRegi, VarImp, use_cores, oob_surv_matrix, oob_residuals);
  }

  return;
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


