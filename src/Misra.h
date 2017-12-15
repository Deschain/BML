/*
 Copyright (c) 2012-2014 Navodit Misra.

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

/* Common includes */
#include <vector>
#include <string>
#include <Rcpp.h>

//A Tree is a vector of these NODES, it is also our data matrix.
struct NODE {
	int dad; 			      //	Parent node on the tree.
	std::vector<int> child;
	std::vector<int> mut; 	  //	# of mutations on the incoming edge.
	std::vector<bool> gtype; //	Binary genotype (mutated -> true).
	int NDes;			      //	Number of samples(taxa) that are descendants.
};

//Local structure of the Bayes Net for each family (parents + param for each gene)
struct FAM {
  std::vector<bool> pi;     //	Parents of a given gene.
  std::vector<bool> potpi;
  std::vector<double> param;
	double nmut;
	double ent;
	double score;
	int npa;
	int order; 		      //  Gene order during OBS
	bool npotpa;
};

struct SimVar {
  std::string vname;
  std::vector<int> vpi;
  std::vector<double> vparam;
	int indx;
};

/*
void InitNetRel(vector<vector<int> >& netrel1, vector<vector<int> >& netrel2, vector<vector<int> >& netrel3, vector<vector<string> >& EdList, vector<SimVar>& SimDAG);
void TPFPAnalysis( vector<vector<int> >& NetRel1,vector<vector<int> >& NetRel2, vector<vector<int> >& NetRel3, vector<SimVar>& SimDAG, int& NRep, string& job);
void ConfigParam(vector<SimVar>& SimDAG, vector<bool>& config, double& pratio);
void WriteLandscape(vector<SimVar>& SimDAG, vector<NODE>& Tree, vector<string>& GeneLabels, double& pthres, string& job);
void LandBoot(vector<SimVar>& SimDAG, vector< vector<double> >& landrel, vector<FAM>& DAG);
void EdBoot(vector< vector<string> >& EdList, vector<SimVar>& SimDAG, vector< vector<double> >& edrel, vector<FAM>& DAG);
void ProbEdBoot(vector< vector<string> >& EdList, vector< vector<double> >& edrel, string& job);
void ProbBootstrapAnalysis( vector<vector<double> >& NetRel1,vector<vector<double> >& NetRel2, vector<SimVar>& SimDAG, string& job);
void WriteDAG(vector<FAM>& DAG, vector<string>& GeneLabels, vector< vector<int> >& char_label, char op, string& job);
void ParseDAG(vector<FAM>& DAG, vector<SimVar>& SimDAG, vector<string>& GeneLabels, vector< vector<int> >& char_label);
void GenerateReplicate(vector<NODE>& Tree, vector<vector<bool> >& Data);
void  SimulateDAG( vector<NODE>& Tree, vector<SimVar>& SimDAG, vector< vector<bool> >& sdata);
void FreshDAG(vector<FAM>& DAG);
bool TestMI(int i, int j,  vector<NODE>& Tree );
void InitEnt(  vector<FAM>& BN, vector<NODE>& Phylogeny);
void LocalPruning( vector<FAM>& DAG, vector<NODE>& Tree );
void InitDAG(vector<FAM>& DAG, vector<FAM>& BN);
bool MScore(double& olap, double& Atot, double& Btot, double& NSam);
void GlobalPruning( vector<FAM>& BayesN, vector<NODE>& Tree);
void SanityCheckTree(vector<NODE>& Tree);
void InitTree( vector<NODE>& Tree, vector<FAM>& BN);
void Delete_Col(vector< vector<bool> >& Data, int col_no);
void EliminateRedundantCharacters(vector< vector<bool> >& Data, vector<vector<int> >& char_label);
bool IsInteger(string& inp, int& Num);
void ParseDataMatrix(vector<vector<bool> >& Data, vector<string>& TaxonLabels, vector<string>& GeneLabels, string& dname);
void DFS(vector<bool>& conf, vector<bool>& seen, int& label, vector< vector< double> >& cdata, bool& insta);
void K2Search( FAM& FamG, vector<FAM>& DAG, vector<NODE>& Tree, int g, char op);
void SearchDAGs( vector<FAM>& DAG, vector<NODE>& Tree, char op );
void RelabelOutlier(vector<NODE>& Tree,  int gn, int pdad, bool& tnew);
void swapnodes(vector<NODE>& Tree, vector<FAM>& DAG, int& i, int& gramps, int& spart, bool& change);
void SearchTrees( vector<FAM>& DAG, vector<NODE>& Tree);
 */
Rcpp::List StructLearn(std::vector<std::vector<bool> >& Data,
                       std::vector<std::string>& GeneLabels,
                       std::vector<SimVar>& SimDAG,
                       int& NTrees,
                       double& pthres);

Rcpp::List BootstrapAnalysis(std::vector<std::vector<bool> >& Data,
                             std::vector<SimVar>& SimDAG,
                             int& NRep);


