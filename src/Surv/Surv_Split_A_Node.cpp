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
# include <Rmath.h>
using namespace Rcpp;

// my header file
# include "..//survForest.h"
# include "..//utilities.h"

void Surv_Split_A_Node(TREENODE* Node,
                       //const double** X,
                       std::vector<  colvec > X,
                       const ivec Y,
                       const ivec Censor,
                       const ivec Ncat,
                       const vec Interval,
                       const PARAMETERS* myPara,
                       const vec subjectweight,
                       ivec useObs,
                       const int node_n,
                       vec variableweight,
                       ivec variableindex,
                       const int P)
{

  auto t1 = std::chrono::system_clock::now();
  //std::chrono::duration<double> diff2 = t3-t2;
  //Rcout << "Time to run pre-split: " << diff2.count() << std::endl;
  int nmin = myPara->nmin;
  int i;

  // calculate node information

  int node_fail = 0;

  for (i = 1; i<node_n; i++)
    node_fail += Censor[useObs[i]];

  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> diff1 = t2-t1;
  //Rcout << "Time to run pre-split: " << diff1.count() << std::endl;
  
  if (node_fail == 0 || node_n <= 2*nmin)
  {
    TERMINATE:;

    Node->Var = -1;			  // terminal node
    Node->NodeObs = useObs;
    Node->NodeSize = node_n;//

  }else{

    int splitVar = -1;
    double splitVal = 0;
    auto t3 = std::chrono::system_clock::now();
    //std::chrono::duration<double> diff2 = t3-t2;
    //Rcout << "Time to run pre-split: " << diff2.count() << std::endl;
    Surv_Find_A_Split(&splitVar, &splitVal, X, Y, Censor, Ncat, Interval, myPara, subjectweight, useObs, node_n, variableweight, variableindex, P);
    auto t4 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff2 = t4-t3;
    //Rcout << "Time to split: " << diff2.count() << std::endl;
    

    
    auto t5 = std::chrono::system_clock::now();
    if (splitVar == -1) // didnt find anything
      goto TERMINATE;
    

    // calculate left and right child node

    int LeftSize = 0;
    int RightSize = 0;

    //int* useObsLeft = (int *) malloc(node_n * sizeof(int));
    ivec useObsLeft(node_n);
    //int* useObsRight = (int *) malloc(node_n * sizeof(int));
    ivec useObsRight(node_n);
    
    if (Ncat[splitVar] > 1)
    {

      //int* goright = (int *) malloc(Ncat[splitVar]*sizeof(int));
      ivec goright(Ncat[splitVar]);
      unpack(splitVal, Ncat[splitVar], goright);

      for (i = 0; i<node_n; i++)
      {
        if ( goright[ (int) X[splitVar][useObs[i]] - 1 ] == 0 )
        {
          useObsLeft[LeftSize] = useObs[i];
          LeftSize ++;
        }else{
          useObsRight[RightSize] = useObs[i];
          RightSize ++;
        }
      }
      useObsLeft.resize(LeftSize);
      useObsRight.resize(RightSize);

      //free(goright);
      //delete[] goright;

    }else{

      for (i = 0; i<node_n; i++)
      {
        if (X[splitVar][useObs[i]] <= splitVal)
        {
          useObsLeft[LeftSize] = useObs[i];
          LeftSize ++;
        }else{
          useObsRight[RightSize] = useObs[i];
          RightSize ++;
        }
      }
    }
    useObsLeft.resize(LeftSize);
    useObsRight.resize(RightSize);

    if (LeftSize == 0 || RightSize == 0)
    {
      //free(useObsLeft);
      //delete[] useObsLeft;
      //free(useObsRight);
      //delete[] useObsRight;
      
      R_DBP("Did not produce a proper split at node %i, for variable %i, need to check node \n", Node, splitVar);
      goto TERMINATE;
    }

    // initiate left and right node

    Node->Var = splitVar;			// Splitting variable
    Node->Val = splitVal;		  // Splitting value

    //Node->Left = (TREENODE*) malloc(sizeof(TREENODE));
    Node->Left = new TREENODE();
    //Node->Right = (TREENODE*) malloc(sizeof(TREENODE));
    Node->Right = new TREENODE();
    auto t6 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff3 = t6-t5;
    //Rcout << "Test of time: " << diff3.count() << std::endl;
    
    //free(useObs);
    //delete[] useObs;

    Surv_Split_A_Node(Node->Left, X, Y, Censor, Ncat, Interval, myPara,
                      subjectweight, useObsLeft, LeftSize, variableweight, variableindex, P);

    Surv_Split_A_Node(Node->Right, X, Y, Censor, Ncat, Interval, myPara,
                      subjectweight, useObsRight, RightSize, variableweight, variableindex, P);

  }

  return;
}