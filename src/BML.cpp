// [[Rcpp::plugins("cpp11")]]

#include "Misra.h"

using namespace std;
using namespace Rcpp;

typedef vector<vector<bool> > tBoolMaxtrix;

tBoolMaxtrix transformMatrix(IntegerMatrix ma)
{
  tBoolMaxtrix data;
  
  int ncol = ma.ncol();
  int nrow = ma.nrow();
  
  for(int row_idx = 0; row_idx < nrow; ++row_idx)
  {
    vector<bool> aux_row(ncol, false);
    for(int col_idx = 0; col_idx < ncol; ++col_idx)
    {
      if(ma(row_idx, col_idx) == 1) aux_row[col_idx] = true;
    }
    data.push_back(aux_row);
  }
  return data;
}

// [[Rcpp::export]]
List BML(IntegerMatrix ma, int ntrees, double pthres) 
{
  // Copied from original code
  srand(12345);
  
  //Input data
  tBoolMaxtrix data = transformMatrix(ma);
  vector<string> names = as<vector<string> >(colnames(ma));
  vector<SimVar> SimDAG;
  
  List result = StructLearn(data, 
                            names, 
                            SimDAG, 
                            ntrees , 
                            pthres);
  
  // List simdag_data = List::create();
  // 
  // for(SimVar tmp : SimDAG)
  // {
  //   simdag_data[tmp.vname] = List::create(
  //     Named("vpi") = tmp.vpi,
  //     Named("vparam") = tmp.vparam,
  //     Named("indx") = tmp.indx
  //   );
  // };
  // 
  // result["Data"] = simdag_data;
  
  return(result);
}