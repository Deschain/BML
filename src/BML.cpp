/*
 BML is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 BML is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BML. If not, see <http://www.gnu.org/licenses/>.
 */

// [[Rcpp::plugins("cpp11")]]

#include <fstream>
#include <iostream>
#include <string>

#include "Misra.h"

using namespace std;
using namespace Rcpp;

typedef vector<vector<bool> > tBoolMaxtrix;

tBoolMaxtrix transformMatrix(IntegerMatrix ma)
{
 
  int ncol = ma.ncol();
  int nrow = ma.nrow();

  tBoolMaxtrix data(nrow, vector<bool>(ncol, false));
  
  for(int row_idx = 0; row_idx < nrow; ++row_idx)
  {
    for(int col_idx = 0; col_idx < ncol; ++col_idx)
    {
      if(ma(row_idx, col_idx) == 1) data[row_idx][col_idx] = true;
    }
  }
  return data;
}


List transformSimData(const vector<SimVar> &simDAG)
{
  List simdag_data = List::create();
  
  for(SimVar tmp : simDAG)
  {
    simdag_data[tmp.vname] = List::create(
      Named("vpi") = tmp.vpi,
      Named("vparam") = tmp.vparam,
      Named("indx") = tmp.indx
    );
  };
 return simdag_data;
}

// [[Rcpp::export]]
List BML(IntegerMatrix ma, int ntrees, double pthres, int nrep) 
{
  // Copied from original code
  srand(12345);
  
  //Input data
  tBoolMaxtrix data = transformMatrix(ma);
  vector<string> names = as<vector<string> >(colnames(ma));
  vector<SimVar> SimDAG;

  List result = StructLearn(data, names, SimDAG, ntrees, pthres);
  
  if(nrep > 0)
  {
    result["bootstrap"] = BootstrapAnalysis(data, SimDAG, nrep);
  }
  
  return(result);
}

// [[Rcpp::export]]
void writeDotFile(List bml, std::string file)
{

  List dag = bml["DAG"];
  
  vector<int> nodes = as<vector<int> >(dag["nodes"]);
  vector<double> probs = as<vector<double> >(dag["probs"]);
  vector<string> labels = as<vector<string> >(dag["labels"]);
  vector<int> edge_org = as<vector<int> >(dag["edges_1"]);
  vector<int> edge_dst = as<vector<int> >(dag["edges_2"]);
  
  fstream dotStream;
  dotStream.open(file.c_str(),ios::out);
  dotStream << "digraph G {" << endl << "\trankdir= LR;"<< endl;
  dotStream << "\t0 [label=\"Normal\" ,shape=ellipse, style=filled, color=\"0.000,1.0,0.85\",fontsize=\"36\", font=\"Helvetica\", fontcolor=white];" << endl;
  
  for(size_t node_idx = 1; node_idx < nodes.size(); ++node_idx)
  {
    
    dotStream << nodes[node_idx];
    dotStream << "\t[label=\"" << labels[node_idx] << "\" ,shape=ellipse, style=filled, color=\"";
    
    switch(std::count(labels[node_idx].begin(),  labels[node_idx].end(), ','))
    {
      case 0:
        dotStream << "0.07500," << probs[node_idx] << ",0.90\"";
        break;
      case 1:
        dotStream << "0.11000," << probs[node_idx] << ",0.90\"";
        break;
      case 2:
        dotStream << "0.13000," << probs[node_idx] << ",0.99\"";
        break;
    }
    dotStream << ",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];" << endl;
  }
  dotStream << endl;
  for(size_t edge_idx = 0; edge_idx < edge_org.size(); ++edge_idx)
  {
    dotStream << "\t" << edge_org[edge_idx] << "->" << edge_dst[edge_idx] << endl;
  }
  
  dotStream << "}";
  dotStream.close();
}