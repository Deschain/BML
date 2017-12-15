/*
 *  Algorithm.cpp
 *  
 *
 *  Created by Navodit Misra.
 *
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


/*	This is the main algorithm that controls and delegates specific tasks to other programs.		*/

/* Generic includes */
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <string>

/* Data Structures & Utils */
#include "Misra.h"

using namespace std;
using namespace Rcpp;

void InitNetRel(vector<vector<int> >& netrel1,vector<vector<int> >& netrel2,vector<vector<int> >& netrel3, vector< vector<string> >& EdList, vector<SimVar>& SimDAG)
{
  for(int i=0;i<SimDAG.size();i++)
  {
    vector<int> temprel;
    for(int j=0;j<SimDAG.size();j++)
    {
      temprel.push_back(0);
    }
    netrel1.push_back(temprel);
    netrel2.push_back(temprel);
    netrel3.push_back(temprel);
    
    string nm = SimDAG[i].vname;
    
    for(int j=0;j<SimDAG[i].vpi.size();j++)
    {
      int ind=SimDAG[i].vpi[j];
      netrel3[i][ind]=1;
      netrel3[ind][i]=1;
      vector<string> panm;
      panm.push_back(SimDAG[ind].vname);
      panm.push_back(nm);
      EdList.push_back(panm);
    }
  }
}

List TPFPAnalysis( vector<vector<int> >& NetRel1,vector<vector<int> >& NetRel2, vector<vector<int> >& NetRel3, vector<SimVar>& SimDAG, int& NRep)
{
  // fstream relFile;
  // string tp="./output/TruePositives"+job;
  // string fp= "./output/FalsePositives"+job;
  // relFile.open(tp.c_str(),ios::out);
  // relFile<<'\t'<<"Tree+DAG "<<'\t'<<"DAG"<<endl;
  // 
  // 
  // fstream TPconf;
  // string tpc="./output/ConfidenceTruePositives"+job;
  // TPconf.open(tpc.c_str(),ios::out);
  // 
  // fstream FPconf;
  // string fpc="./output/ConfidenceFalsePositives"+job;
  // FPconf.open(fpc.c_str(),ios::out);
  
  vector<int> Hist1, Hist2;
  
  IntegerVector truePositives_tree_dag;
  IntegerVector truePositives_dag;
  CharacterVector truePositives_names;
  
  for(int i=0;i<10;i++)
  {
    Hist1.push_back(0);
    Hist2.push_back(0);
    
  }
  for(int i=0;i<NetRel1.size()-1;i++)
  {
    
    for(int j=i+1;j<NetRel1.size();j++)
    {
      if(NetRel3[i][j]==1)
      {
        // relFile<<SimDAG[i].vname<<"-"<<SimDAG[j].vname<<'\t';
        // relFile<<NetRel1[i][j]<<'\t'<<NetRel2[i][j];
        // relFile<<endl;
        
        truePositives_names.push_back(SimDAG[i].vname+"-"+SimDAG[j].vname);
        truePositives_tree_dag.push_back(NetRel1[i][j]);
        truePositives_dag.push_back(NetRel2[i][j]);
        
        int bin1=int(10*NetRel1[i][j]/NRep);
        int bin2=int(10*NetRel2[i][j]/NRep);
        if(bin1==10){bin1=9;}
        if(bin2==10){bin2=9;}
        Hist1[bin1]++;
        Hist2[bin2]++;
      }
    }
  }
  // relFile<<endl<<endl;
  // relFile.close();
  IntegerMatrix truePositives(truePositives_tree_dag.size(), 2);
  truePositives(_, 0) = truePositives_tree_dag;
  truePositives(_, 1) = truePositives_dag;
  rownames(truePositives) = truePositives_names;
  colnames(truePositives) = CharacterVector({"Tree-DAG", "DAG"});
  
  // TPconf<<"Confidence"<<'\t'<<"Tree+DAG"<<'\t'<<"DAG"<<endl;
  IntegerMatrix ConfidenceTruePositives(9,3);
  
  int nted=0;
  int noed=0;
  for(int i=9;i>0;i--)
  {
    nted+=Hist1[i];
    noed+=Hist2[i];
    // TPconf<<10*i<<'\t'<<nted<<'\t'<<noed<<endl;
    ConfidenceTruePositives(9-i, 0) = 10*i;
    ConfidenceTruePositives(9-i, 1) = nted;
    ConfidenceTruePositives(9-i, 2) = noed;
    Hist1[i]=0;
    Hist2[i]=0;
  }
  Hist1[0]=0;
  Hist2[0]=0;
  // TPconf.close();
  colnames(ConfidenceTruePositives) = CharacterVector({"Confidence", "Tree+DAG", "DAG"});
  
  // fstream relFile2;
  // relFile2.open(fp.c_str(),ios::out);
  // relFile2<<'\t'<<"Tree+DAG "<<'\t'<<"DAG"<<endl;
  IntegerVector falsePositives_tree_dag;
  IntegerVector falsePositives_dag;
  CharacterVector falsePositives_names;
  
  for(int i=0;i<NetRel1.size()-1;i++)
  {
    
    for(int j=i+1;j<NetRel1.size();j++)
    {
      if((NetRel1[i][j]+NetRel2[i][j]>0)&&(NetRel3[i][j]==0))
      {
        // relFile2<<SimDAG[i].vname<<"-"<<SimDAG[j].vname<<'\t';
        // relFile2<<NetRel1[i][j]<<'\t'<<NetRel2[i][j];
        // relFile2<<endl;
        falsePositives_names.push_back(SimDAG[i].vname+"-"+SimDAG[j].vname);
        falsePositives_tree_dag.push_back(NetRel1[i][j]);
        falsePositives_dag.push_back(NetRel2[i][j]);
        
        int bin1=int(10*NetRel1[i][j]/NRep);
        int bin2=int(10*NetRel2[i][j]/NRep);
        if(bin1==10){bin1=9;}
        if(bin2==10){bin2=9;}
        Hist1[bin1]++;
        Hist2[bin2]++;
      }
    }
  }
  // relFile2.close();
  IntegerMatrix falsePositives(falsePositives_tree_dag.size(),2);
  falsePositives(_,0) = falsePositives_tree_dag;
  falsePositives(_,1) = falsePositives_dag;
  rownames(falsePositives) = falsePositives_names;
  colnames(falsePositives) = CharacterVector({"Tree-DAG", "DAG"});
  
  nted=0;
  noed=0;
  // FPconf<<"Confidence"<<'\t'<<"Tree+DAG"<<'\t'<<"DAG"<<endl;
  IntegerMatrix confidenceFalsePositives(9,3);
  for(int i=9;i>0;i--)
  {
    nted+=Hist1[i];
    noed+=Hist2[i];
    // FPconf<<10*i<<'\t'<<nted<<'\t'<<noed<<endl;
    confidenceFalsePositives(9-i,0) = 10*i;
    confidenceFalsePositives(9-i,1) = nted;
    confidenceFalsePositives(9-i,2) = noed;
    Hist1[i]=0;
    Hist2[i]=0;
  }
  Hist1[0]=0;
  Hist2[0]=0;
  // FPconf.close();
  colnames(confidenceFalsePositives) = CharacterVector({"Confidence", "Tree+DAG", "DAG"});
  
  List result = List::create();
  
  result["TruePositives"] = truePositives;
  result["ConfidenceTruePositives"] = ConfidenceTruePositives;
  result["FalsePositives"] = falsePositives;
  result["ConfidenceFalsePositives"] = confidenceFalsePositives;
  
  return result;
}

void ConfigParam(vector<SimVar>& SimDAG, vector<bool>& config, double& pratio)
{
  for(int l=0;l<SimDAG.size();l++)
  {
    int y=0;
    int y2=SimDAG[l].vpi.size();
    for(int pa=0;pa<y2;pa++)
    {
      y+=config[SimDAG[l].vpi[pa]]*(int(pow(2.0,pa)));
      
    }
    double z=SimDAG[l].vparam[y]/(SimDAG[l].vparam[y]+SimDAG[l].vparam[y+(int(pow(2.0,y2)))] );
    if(config[l]==true){pratio=pratio*(1.0-z);}
    else{pratio=pratio*z;}
  }
}


List WriteLandscape(vector<SimVar>& SimDAG, vector<NODE>& Tree, vector<string>& GeneLabels, vector< vector<int> >& char_label, double& pthres)
{
  double epsilon=0.0001;
  int NVars= SimDAG.size();
  vector<double> NormStat;
  vector<bool>config(SimDAG.size());
  
  // fstream pathFile;
  // string bpath="./output/Path"+job+".dot";
  // pathFile.open(bpath.c_str(),ios::out);
  for(int i=0;i<SimDAG.size();i++)
  {
    config[i]=false;
  }
  double pnor=1.0;
  for(int i=0;i<SimDAG.size();i++)
  {
    
    double tempsc=0.0;
    tempsc=SimDAG[i].vparam[0];
    NormStat.push_back(tempsc);
    pnor=pnor*tempsc;
  }
  
  
  vector<int> nodes;
  vector<string> labels;
  vector<double> probs;
  vector<int> edges_1;
  vector<int> edges_2;
  
  // pathFile<<"digraph G{"<<endl<<"rankdir= LR;"<<endl;
  // pathFile<<0<<"[label=\"Normal\" ,shape=ellipse, style=filled, color=\"0.000,1.0,0.85\",fontsize=\"36\", font=\"Helvetica\", fontcolor=white];"<<endl;
  nodes.push_back(0);
  labels.push_back("Normal");
  probs.push_back(1.0);
  
  double max=0.0;
  vector<double>  prob;
  vector<string> plabel;
  vector<bool> present;
  for(int i=0;i<SimDAG.size();i++)
  {
    config[i]=true;
    double pratio=1.0;
    ConfigParam( SimDAG,  config,  pratio);
    if(pratio>max){max=pratio;}
    prob.push_back(pratio);
    present.push_back(false);
    config[i]=false;
  }
  
  vector<double> prob2;
  vector< vector<int> > plabel2;
  vector<bool> present2;
  double max2=0.0;
  for(int i=0;i<SimDAG.size()-1;i++)
  {
    config[i]=true;
    for(int j=i+1;j<SimDAG.size();j++)
    {
      config[j]=true;
      
      double pratio=1.0;
      ConfigParam( SimDAG,  config,  pratio);
      
      
      prob2.push_back(pratio);
      vector<int> lab;
      lab.push_back(i);
      lab.push_back(j);
      plabel2.push_back(lab);
      present2.push_back(false);
      
      if(pratio>max2){max2=pratio;}
      
      config[j]=false;
    }
    config[i]=false;
  }
  
  vector<double> prob3;
  vector< vector<int> > plabel3;
  vector<bool> present3;
  double max3=0.0;
  for(int i=0;i<SimDAG.size()-2;i++)
  {
    config[i]=true;
    for(int j=i+1;j<SimDAG.size()-1;j++)
    {
      
      config[j]=true;
      
      
      for(int k=j+1;k<SimDAG.size();k++)
      {
        config[k]=true;
        bool pre=false;
        for(int pr=1;pr<Tree[0].NDes+1;pr++)
        {
          if(Tree[pr].gtype[SimDAG[i].indx]==true && Tree[pr].gtype[SimDAG[j].indx]==true && Tree[pr].gtype[SimDAG[k].indx]==true)
          {
            pre=true;
          }
        }
        present3.push_back(pre);
        double pratio=1.0;
        ConfigParam( SimDAG,  config,  pratio);
        
        
        prob3.push_back(pratio);
        vector<int> lab;
        lab.push_back(i);
        lab.push_back(j);
        lab.push_back(k);
        plabel3.push_back(lab);
        
        if(pratio>max3){max3=pratio;}
        config[k]=false;
        
      }
      
      config[j]=false;
    }
    config[i]=false;
  }
  for(int i=0;i<prob3.size();i++)
  {
    
    if( present3[i]==true && prob3[i]>=max3*(pthres-epsilon))
    {
      // pathFile<<i+1+NVars+NVars*(NVars-1)/2
      // <<"[label=\""<<SimDAG[plabel3[i][0]].vname<<","<<SimDAG[plabel3[i][1]].vname<<","<<SimDAG[plabel3[i][2]].vname
      // <<"\" ,shape=ellipse, style=filled, color=\"0.13000,"<<prob3[i]/max3<<",0.99\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
      nodes.push_back(i+1+NVars+NVars*(NVars-1)/2);
      labels.push_back(SimDAG[plabel3[i][0]].vname+","+SimDAG[plabel3[i][1]].vname+","+SimDAG[plabel3[i][2]].vname);
      probs.push_back(prob3[i]/max3);
      
      int pa;
      double papr=0;
      for(int j=0;j<prob2.size();j++)
      {
        if( ( (plabel2[j][0]==plabel3[i][0]) && ( (plabel2[j][1]==plabel3[i][1]) || (plabel2[j][1]==plabel3[i][2]) ) ) || ( (plabel2[j][0]==plabel3[i][1]) &&  (plabel2[j][1]==plabel3[i][2])  ) )
        {
          if(prob2[j]>papr)
          {
            pa=j;
            papr=prob2[j];
          }
        }
      }
      present2[pa]=true;
      // pathFile<<NVars+pa+1<<"->"<<i+1+NVars+NVars*(NVars-1)/2<<endl;
      edges_1.push_back(NVars+pa+1);
      edges_2.push_back(i+1+NVars+NVars*(NVars-1)/2);
    }
  }
  
  for(int i=0;i<prob2.size();i++)
  {
    int g1 = plabel2[i][0];
    int g2 = plabel2[i][1];
    if(prob2[i]>=max2*(pthres-epsilon))
    {
      for(int pr=1;pr<Tree[0].NDes+1;pr++)
      {
        bool pr2=false;
        
        if(Tree[pr].gtype[SimDAG[g1].indx]==true && Tree[pr].gtype[SimDAG[g2].indx]==true )
        {
          pr2=true;
        }
        int mcnt=0;
        for(int j=0;j<Tree[pr].gtype.size();j++)
        {
          if( Tree[pr].gtype[j]==true){mcnt++;}
        }
        if(pr2==true && mcnt ==2){present2[i]=true;}
      }
    }
    
    if(present2[i]==true )
    {
      // pathFile<<NVars+i+1<<"[label=\""<<SimDAG[g1].vname<<","<<SimDAG[g2].vname<<"\" ,shape=ellipse, style=filled, color=\"0.11000,"<<prob2[i]/max2<<",0.90\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
      nodes.push_back(NVars+i+1);
      labels.push_back(SimDAG[g1].vname+","+SimDAG[g2].vname);
      probs.push_back(prob2[i]/max2);
      
      if(prob[g1]>prob[g2])
      {
        // pathFile<<g1+1<<"->"<<NVars+i+1<<endl;
        edges_1.push_back(g1+1);
        present[g1]=true;
      }else
      {
        // pathFile<<g2+1<<"->"<<NVars+i+1<<endl;
        edges_1.push_back(g2+1);
        present[g2]=true;
      }
      edges_2.push_back(NVars+i+1);
    }
    
  }
  
  for(int i=0;i<prob.size();i++)
  {
    if(prob[i]>=max*(pthres-epsilon))
    {
      for(int pr=1;pr<Tree[0].NDes+1;pr++)
      {
        bool pr1=false;
        
        if(Tree[pr].gtype[SimDAG[i].indx]==true )
        {
          pr1=true;
        }
        int mcnt=0;
        for(int j=0;j<Tree[pr].gtype.size();j++)
        {
          if( Tree[pr].gtype[j]==true){mcnt++;}
        }
        if(pr1==true && mcnt ==1){present[i]=true;}
      }
    }
    
    if(present[i]==true)
    {
      // pathFile<<0<<"->"<<i+1<<endl;
      // pathFile<<i+1<<"[label=\""<<SimDAG[i].vname<<"\" ,shape=ellipse, style=filled, color=\"0.075000,"<<prob[i]/max <<",0.90\",fontsize=\"36\", font=\"Helvetica\", fontcolor=black];"<<endl;
      nodes.push_back(i+1);
      labels.push_back(SimDAG[i].vname);
      probs.push_back(prob[i]/max);
      edges_1.push_back(0);
      edges_2.push_back(i+1);
    }
    
  }
  
  // pathFile<<"}";
  // pathFile.close();
  return List::create(
    Named("nodes") = wrap(nodes),
    Named("probs") = wrap(probs),
    Named("labels") = wrap(labels),
    Named("edges_1") = wrap(edges_1),
    Named("edges_2") = wrap(edges_2)
  );
}

void LandBoot(vector<SimVar>& SimDAG, vector< vector<double> >& landrel, vector<FAM>& DAG)
{
  int NVars=SimDAG.size();
  vector<string>VarName(NVars);
  vector<vector<int> > VarPi(NVars);
  vector<vector<double> >VarParam(NVars);
  vector<bool>config(DAG.size());
  for(int i=0;i<NVars;i++)
  {
    string vname;
    int ord,npi;
    vname=SimDAG[i].vname;
    ord=DAG[i].order;
    npi=DAG[i].npa;
    VarName[ord]=vname;
    for(int j=0;j<NVars;j++)
    {
      if(DAG[i].pi[j]==true){
        VarPi[ord].push_back(DAG[j].order);}
    }
    int nparam=int(pow(2.0,npi+1));
    for(int j=0;j<nparam;j++)
    {
      double vparam=DAG[i].param[j];
      VarParam[ord].push_back(vparam);
    }
  }
  double pnor=1.0;
  for(int i=0;i<DAG.size();i++)
  {
    config[i]=false;
    pnor=pnor*DAG[i].param[0];
  }
  double max=0.0;
  vector<double>  prob;
  vector<string> plabel;
  for(int ng=0;ng<NVars;ng++)
  {
    
    for(int i=0;i<DAG.size();i++)
    {
      if(VarName[i]==SimDAG[ng].vname)
      {
        config[i]=true;
        double pratio=1.0;
        for(int l=0;l<DAG.size();l++)
        {
          int y=0;
          int y2=VarPi[l].size();
          for(int pa=0;pa<VarPi[l].size();pa++)
          {
            y+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));
            
          }
          double z=VarParam[l][y]/(VarParam[l][y]+VarParam[l][y+(int(pow(2.0,y2)))] );
          if(config[l]==true){pratio=pratio*(1.0-z);}
          else{pratio=pratio*z;}
          
          
        }
        if(pratio>max){max=pratio;}
        prob.push_back(pratio);
        config[i]=false;
      }
    }
  }
  
  for(int i=0;i<prob.size();i++)
  {
    landrel[i].push_back(prob[i]/(pnor));
  }
  
}

void EdBoot(vector< vector<string> >& EdList, vector<SimVar>& SimDAG, vector< vector<double> >& edrel, vector<FAM>& DAG)
{
  int NVars=DAG.size();
  vector<string>VarName(NVars);
  vector<vector<int> > VarPi(NVars);
  vector<vector<double> >VarParam(NVars);
  vector<bool>config(DAG.size());
  for(int i=0;i<NVars;i++)
  {
    string vname;
    int ord,npi;
    vname=SimDAG[i].vname;
    ord=DAG[i].order;
    npi=DAG[i].npa;
    VarName[ord]=vname;
    for(int j=0;j<NVars;j++)
    {
      if(DAG[i].pi[j]==true){
        VarPi[ord].push_back(DAG[j].order);}
    }
    int nparam=int(pow(2.0,npi+1));
    for(int j=0;j<nparam;j++)
    {
      double vparam=DAG[i].param[j];
      VarParam[ord].push_back(vparam);
    }
  }
  double pnor=1.0;
  for(int i=0;i<DAG.size();i++)
  {
    config[i]=false;
    pnor=pnor*DAG[i].param[0];
  }
  for(int ng=0;ng<EdList.size();ng++)
  {
    vector<double>  prob(3);
    
    vector<int> glab(2);
    
    for(int i=0;i<DAG.size();i++)
    {
      if(VarName[i]==EdList[ng][0] || VarName[i]==EdList[ng][1] )
      {
        
        config[i]=true;
        double pratio=1.0;
        for(int l=0;l<DAG.size();l++)
        {
          int y=0;
          int y2=VarPi[l].size();
          for(int pa=0;pa<VarPi[l].size();pa++)
          {
            y+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));
            
          }
          double z=VarParam[l][y]/(VarParam[l][y]+VarParam[l][y+(int(pow(2.0,y2)))]);
          if(config[l]==true){pratio=pratio*(1.0-z);}
          else{pratio=pratio*z;}
          
          
        }
        if(VarName[i]==EdList[ng][0])
        {
          glab[0]=i;
          prob[0]=pratio;
        }else
        {
          glab[1]=i;
          prob[1]=pratio;
        }
        config[i]=false;
      }
    }
    config[glab[0]]=true;
    config[glab[1]]=true;
    
    double pratio=1.0;
    for(int l=0;l<DAG.size();l++)
    {
      int y=0;
      int y2=VarPi[l].size();
      
      for(int pa=0;pa<VarPi[l].size();pa++)
      {
        y+=config[VarPi[l][pa]]*(int(pow(2.0,pa)));
      }
      double z=VarParam[l][y]/(VarParam[l][y]+VarParam[l][y+(int(pow(2.0,y2)))]);
      
      if(config[l]==true){pratio=pratio*(1.0-z);}
      else{pratio=pratio*z;}
    }
    prob[2]=pratio;
    
    config[glab[1]]=false;
    config[glab[0]]=false;
    for(int bp=0;bp<3;bp++)
    {
      edrel[3*ng+bp].push_back((prob[bp]/pnor));
    }
  }
}

NumericMatrix ProbEdBoot(vector< vector<string> >& EdList, vector< vector<double> >& edrel)
{
  // string edb= "./output/EdgeProbabilities"+job;
  // fstream edFile;
  // edFile.open(edb.c_str(),ios::out);
  CharacterVector edge_names;
  NumericMatrix   edge_probs(4*EdList.size(), edrel[0].size());
  
  for(int ng=0;ng<EdList.size();ng++)
  {
    // edFile<<EdList[ng][0];
    edge_names.push_back(EdList[ng][0]);
    for(int bp=0;bp< edrel[0].size();bp++)
    {
      // edFile<<'\t'<<edrel[3*ng][bp];
      edge_probs(ng*4 + 0, bp) = edrel[3*ng][bp];
    }
    // edFile<<endl;
    // edFile<<EdList[ng][1];
    
    edge_names.push_back(EdList[ng][1]);
    for(int bp=0;bp< edrel[0].size();bp++)
    {
      // edFile<<'\t'<<edrel[3*ng+1][bp];
      edge_probs(ng*4 + 1, bp) = edrel[3*ng+1][bp];
    }
    // edFile<<endl;
    // edFile<<EdList[ng][0]<<","<<EdList[ng][1];
    edge_names.push_back(EdList[ng][0]+","+EdList[ng][1]);
    for(int bp=0;bp< edrel[0].size();bp++)
    {
      // edFile<<'\t'<<edrel[3*ng+2][bp];
      edge_probs(ng*4 + 2, bp) = edrel[3*ng+2][bp];
    }
    // edFile<<endl;
    // edFile<<EdList[ng][0]<<","<<EdList[ng][1]<<"(indep)";
    edge_names.push_back(EdList[ng][0]+","+EdList[ng][1]+"(indep)");
    for(int bp=0;bp< edrel[0].size();bp++)
    {
      // edFile<<'\t'<<edrel[3*ng+1][bp]*edrel[3*ng][bp];
      edge_probs(ng*4 + 3, bp) = edrel[3*ng+1][bp]*edrel[3*ng][bp];
    }
    // edFile<<endl;
  }
  // edFile.close();
  rownames(edge_probs) = edge_names;
  return edge_probs;
}

List ProbBootstrapAnalysis( vector<vector<double> >& NetRel1,vector<vector<double> >& NetRel2, vector<SimVar>& SimDAG)
{
  // string lr1= "./output/OBS_Probabilities"+job;
  // string lr2= "./output/Tree_OBS_Probabilities"+job;
  // 
  // fstream relFile1;
  // fstream relFile2;
  // 
  // relFile1.open(lr1.c_str(),ios::out);
  // relFile2.open(lr2.c_str(),ios::out);
  NumericMatrix obsProbabilities(NetRel1.size(),NetRel1[0].size());
  NumericMatrix TreeObsProbabilities(NetRel1.size(),NetRel1[0].size());
  CharacterVector prob_names;
  
  for(int i=0;i<NetRel1.size();i++)
  {
    // relFile1<<SimDAG[i].vname<<'\t';
    // relFile2<<SimDAG[i].vname<<'\t';
    prob_names.push_back(SimDAG[i].vname);
    
    for(int j=0;j<NetRel1[0].size();j++)
    {
      
      // relFile1<<NetRel1[i][j]<<'\t';
      // relFile2<<NetRel2[i][j]<<'\t';
      obsProbabilities(i,j) = NetRel1[i][j];
      TreeObsProbabilities(i,j) = NetRel2[i][j];
      
    }
    // relFile1<<endl;
    // relFile2<<endl;
    
  }
  
  // relFile1.close();
  // relFile2.close();
  rownames(obsProbabilities) =  prob_names;   
  rownames(TreeObsProbabilities) =  prob_names;
  
  List result = List::create();
  result["OBS_Probabilities"] = obsProbabilities;
  result["Tree_OBS_Probabilities"] = TreeObsProbabilities;
  return result;
}
/*	Output the local search based MPE estimates for DAG structure		*
*	for each starting tree.												*/



void WriteDAG(vector<FAM>& DAG, vector<string>& GeneLabels, vector< vector<int> >& char_label, char op, string& job)
{
  ofstream DAGFile;
  fstream outFile;
  
  string st1= "./output/Tree_DAG";
  string dt1= st1 + ".dot";
  string st2= "./output/OBS_DAG";
  string dt2= st2 + ".dot";
  
  if(op=='y')
  {
    DAGFile.open(st1.c_str(),ios::out);
    outFile.open(dt1.c_str(),ios::out);//prepare the output dot file
    
  }else
  {
    DAGFile.open(st2.c_str(),ios::out);
    outFile.open(dt2.c_str(),ios::out);//prepare the output dot file
    
  }
  
  outFile<<"digraph G{"<<endl;
  DAGFile<<DAG.size()<<'\n';
  vector<bool> cflag(DAG.size());
  for(int i=0;i<cflag.size();i++)
  {
    cflag[i]=false;
  }
  for(int i=0;i<DAG.size();i++)
  {
    
    for(int j=0;j< DAG[i].pi.size();j++)
    {
      int wt=0;
      if(DAG[i].pi[j]==1)
      {
        outFile<<j<<"->"<<i<<endl;
        cflag[i]=true;
        cflag[j]=true;
      }
      
    }
    
  }
  for(int i=0; i<DAG.size();i++)
  {
    if(cflag[i]==true)
    {
      outFile<<i<< "[style=filled, color=\"skyblue2\",font=\"Helvetica\",label=\""<<GeneLabels[char_label[i][0]];
      for(int k=1;k<char_label[i].size();k++)
      {
        outFile<<","<<GeneLabels[char_label[i][k]];
      }
      outFile<< "\"]" << endl;
    }
    DAGFile<<GeneLabels[char_label[i][0]]<<'\t'<<DAG[i].order<<'\t'<<DAG[i].npa<<'\t';
    for(int j=0;j< DAG[i].pi.size();j++)
    {
      if(DAG[i].pi[j]==1)
      {
        DAGFile<<DAG[j].order<<'\t';
      }
      
    }
    
    for(int j=0;j< DAG[i].param.size();j++)
    {
      DAGFile<<DAG[i].param[j]<<'\t';
      
    }
    DAGFile<<endl;
    
    
  }
  
  outFile<<"}"<<endl;
  outFile.close();
  DAGFile.close();
  
}



/*	Exception handling and error messages								*/

void ExceptionHandler(int& ex)
{
  switch (ex)
  {
  case 0:
    cerr<<"Could not open the data file. Perhaps you entered an incorrect path or file name.\nBailing Out !\n";
    break;
  case 1:
    cerr<<"Please check the data matrix file.\nSee README file or look at the examples in the data folder for detailed formatting instructions.\nBailing Out !\n";
    break;
  case 2:
    cerr<<"child-parent relations awry. Failed the sanity check in InitTree.cpp";
    break;
  case 3:
    cerr<<"Number of random restarts you entered is not a positive integer.\nPlease restart the analysis. \nBailing out !\n";
    break;
  case 4:
    cerr<<"You must enter a number between 0 and 1 as the threshold for path reconstruction.\nPlease restart the analysis. \nBailing out !\n";
    break;
  case 5:
    cerr<<"You must enter a positive integer for the number of bootstrap replicates.\nPlease restart the analysis. \nBailing out !\n";
    break;
  case 6:
    cerr<<"You must enter 'y' or 'n' for the bootstrap option.\nPlease restart the analysis. \nBailing out !\n";
    break;
  }
}

void ParseDAG(vector<FAM>& DAG, vector<SimVar>& SimDAG, vector<string>& GeneLabels, vector< vector<int> >& char_label)
{
  SimDAG.resize(DAG.size());
  for(int i=0;i<DAG.size();i++)
  {
    SimVar temp;
    temp.vname=GeneLabels[char_label[i][0]];
    int ord=DAG[i].order;
    int npi=DAG[i].npa;
    for(int j=0;j<DAG.size();j++)
    {
      if(DAG[i].pi[j]==true){
        temp.vpi.push_back(DAG[j].order);}
    }
    int nparam=int(pow(2.0,npi+1));
    for(int j=0;j<nparam;j++)
    {
      temp.vparam.push_back(DAG[i].param[j]);
    }
    temp.indx=i;
    SimDAG[ord]=temp;
  }
  
  
}

/*	Read the data file, store taxon and gene names as a vector and partially initialize the variable "Tree".					*/

void GenerateReplicate(vector<NODE>& Tree, vector<vector<bool> >& Data)
{
  NODE v;
  v.dad=0;
  v.NDes=1;
  for(int i=0;i<Data[0].size();i++)
  {
    v.gtype.push_back(false);
  }
  
  Tree.push_back(v);                      //Root initialized    srand(time(0));
  
  for(int i=0;i<Data.size();i++)
  {
    v.gtype=Data[i];
    Tree.push_back(v);
  }
  Tree[0].NDes=Data.size();
}

void  SimulateDAG( vector<NODE>& Tree, vector<SimVar>& SimDAG, vector< vector<bool> >& sdata)
{
  
  int tsamp=0;
  NODE v;
  v.dad=0;
  v.NDes=0;
  
  for(int i=0;i<SimDAG.size();i++)
  {
    v.gtype.push_back(false);
    
  }
  Tree.push_back(v);
  int cof=1;
  while(tsamp<sdata.size())
  {
    vector<bool> data(SimDAG.size());
    
    int cot=0;
    
    for(int j=0;j<SimDAG.size();j++)
    {
      double x;
      x=(double) rand()/RAND_MAX;
      int y1=0;
      int y2=SimDAG[j].vpi.size();
      double z=0.0;
      for(int k=0;k<y2;k++)
      {
        y1+=data[SimDAG[j].vpi[k]]*(int(pow(2.0,k)));
      }
      z=SimDAG[j].vparam[y1]/(SimDAG[j].vparam[y1]+SimDAG[j].vparam[y1+(int(pow(2.0,y2)))] );
      if(SimDAG[j].vparam[y1]+SimDAG[j].vparam[y1+(int(pow(2.0,y2)))] <.999){Rcout<<"Prob norm \n";}
      if(SimDAG[j].vparam[y1]+SimDAG[j].vparam[y1+(int(pow(2.0,y2)))] >1.001){Rcout<<"Prob norm 2 \n";}
      if(z>x)
      {
        data[j]=0;
      }else
      {
        data[j]=1;
        cot++;
      }
    }
    bool ctrue=false;
    int sdt=0;
    
    while(ctrue==false && sdt<sdata.size())
    {
      ctrue=true;
      int snmut=0;
      for(int gch=0;gch<SimDAG.size();gch++)
      {
        if(data[gch]==false && sdata[sdt][SimDAG[gch].indx]==true){ctrue=false;}
        if(sdata[sdt][SimDAG[gch].indx]==true){snmut++;}
        
      }
      sdt++;
      if(snmut<cof){ctrue=false;}
    }
    if(ctrue ==true && cot>= cof)
    {
      tsamp++;
      for(int j=0;j<SimDAG.size();j++)
      {
        v.gtype=data;
      }
      Tree.push_back(v);
    }
  }
}

void FreshDAG(vector<FAM>& DAG)
{
  FAM famtemp;
  
  for(int i=0;i<DAG.size();i++)
  {
    famtemp.pi.push_back(false);
  }
  famtemp.npa=0;
  famtemp.nmut=0.0;
  famtemp.score=0.0;
  for(int i=0;i<DAG.size();i++)
  {
    
    DAG[i].pi=famtemp.pi;
    DAG[i].npa=famtemp.npa;
    DAG[i].score=famtemp.score;
  }
  
}

bool TestMI(int i, int j,  vector<NODE>& Tree )
{
  double cdat[2][2];
  cdat[0][0]=0;
  cdat[0][1]=0;
  cdat[1][0]=0;
  cdat[1][1]=0;
  double count=0;
  double MInf=0.0;
  double nmin1=0.0;
  double nmin2=0.0;
  for(int k=0;k<Tree.size();k++)
  {
    bool ni=Tree[k].gtype[i];
    bool nj=Tree[k].gtype[j];
    cdat[ni][nj]++;
    count++;
  }
  double Minf2=cdat[0][0]*log(((count-cdat[1][0])*cdat[0][0])/((cdat[0][0]+cdat[0][1])*(cdat[0][0])));
  double Minf3=cdat[0][0]*log(((count-cdat[0][1])*cdat[0][0])/((cdat[0][0])*(cdat[0][0]+cdat[1][0])));
  
  
  
  MInf=cdat[0][0]*log((count*cdat[0][0])/((cdat[0][0]+cdat[0][1])*(cdat[0][0]+cdat[1][0])));
  if(cdat[0][1]!=0)
  {
    MInf+= cdat[0][1]*log((count*cdat[0][1])/((cdat[0][0]+cdat[0][1])*(cdat[0][1]+cdat[1][1])));
    Minf2=cdat[0][1]*log(((count-cdat[1][0])*cdat[0][1])/((cdat[0][0]+cdat[0][1])*(cdat[0][1]+cdat[1][1])));
    
  }
  if(cdat[1][0]!=0)
  {
    MInf+= cdat[1][0]*log((count*cdat[1][0])/((cdat[1][0]+cdat[1][1])*(cdat[0][0]+cdat[1][0])));
    Minf3=cdat[1][0]*log(((count-cdat[0][1])*cdat[1][0])/((cdat[1][0]+cdat[1][1])*(cdat[0][0]+cdat[1][0])));
    
  }
  if(cdat[1][1]!=0)
  {
    MInf+=cdat[1][1]*log((count*cdat[1][1])/((cdat[1][0]+cdat[1][1])*(cdat[0][1]+cdat[1][1])));
    Minf2=cdat[0][1]*log(((count-cdat[1][0])*cdat[1][1])/((cdat[1][1]+cdat[0][1])*(cdat[1][1])));
    Minf3=cdat[1][1]*log(((count-cdat[0][1])*cdat[1][1])/((cdat[1][0]+cdat[1][1])*(cdat[1][1])));
    
    nmin1=cdat[1][1]/(cdat[0][1]+cdat[1][1]);
    nmin2=cdat[1][1]/(cdat[1][0]+cdat[1][1]);
  }else
  {
    nmin1=0.0;
    nmin2=0.0;
  }
  if((MInf>log(count)) && (cdat[0][1]/(cdat[0][0]+cdat[0][1])<=nmin2)  )    {
    return true;
  }else
  {
    return false;
    
  }
  //return true;
  
  
}
void InitEnt(  vector<FAM>& BN, vector<NODE>& Phylogeny)
{
  for(int i=0;i<BN.size();i++)
  {
    double nummut=0.0;
    double nvert=0.0;
    for(int k=0;k<Phylogeny.size();k++)
    {
      if(Phylogeny[k].gtype[i]==true){nummut+=1.0;}
      nvert+=1.0;
      
    }
    BN[i].nmut=nummut;
    double theta=nummut/nvert;
    if(theta>0.0){BN[i].ent= - nummut*log(theta)-(nvert-nummut)*log(1.0-theta);}else{BN[i].ent=0.0;}
    BN[i].param.erase(BN[i].param.begin(),BN[i].param.end());
    BN[i].param.push_back(1.0-theta);
    BN[i].param.push_back(theta);
    
  }
}


/*	Initialize the subtree structure and potential parents for each gene																	*
*	in variable "Genes" given "Tree"																										*
*	Identify the set of possible parents for gene "G" given "Tree"																			*/
void LocalPruning( vector<FAM>& DAG, vector<NODE>& Tree )
{
  FreshDAG(DAG);
  for(int i=0;i<DAG.size()-1;i++)
  {
    if(DAG[i].npotpa==true)
    {
      for(int j=i+1;j<DAG.size();j++)
      {
        if( DAG[i].potpi[j]==true)
        {
          if(TestMI(i,j,Tree)==true)
          {
            DAG[j].potpi[i]=true;
          }
          else
          {
            DAG[j].potpi[i]=false;
          }
        }
      }
    }
  }
  InitEnt(DAG,Tree);
}

void InitDAG(vector<FAM>& DAG, vector<FAM>& BN)
{
  for(int i=0;i<BN.size();i++)
  {
    DAG[i]=BN[i];
    
  }
}

bool MScore(double& olap, double& Atot, double& Btot, double& NSam)
{
  double MinfMax=0.0;
  
  double A1B0= (double) 1.0*(Atot-1.0*olap);
  double A0B1= (double) 1.0*(Btot-1.0*olap);
  double A1B1= olap;
  double TotS= (double) 2.0*NSam -2.0;
  double norm= (double) 1.0*(TotS-A1B0-A0B1-1.0*olap);
  double A0B0= norm;
  for(int i=0; i<olap-1; i++)
  {
    A0B0=norm-i;
    A1B1=olap+i;
    double MinfA=0.0;
    if(A0B0>0)
    {
      MinfA+= (double)A0B0*log((1.0*TotS*A0B0)/(1.0*(A0B0+A0B1)*(A0B0+A1B0)));
    }
    
    if(A1B0>0)
    {
      MinfA+= (double)A1B0*log((1.0*TotS*A1B0)/(1.0*(A0B0+A1B0)*(A1B1+A1B0)));
      
    }
    if(A0B1>0)
    {
      MinfA+= (double)A0B1*log((1.0*TotS*A0B1)/(1.0*(A0B0+A0B1)*(A1B1+A0B1)));
    }
    if(A1B1>0)
    {
      MinfA+= (double)A1B1*log((1.0*TotS*A1B1)/(1.0*(A1B1+A1B0)*(A1B1+A0B1)));
    }
    if(MinfA>MinfMax){MinfMax=MinfA;}
  }
  
  
  
  if(MinfMax>log(TotS) )
  {
    return true;
  }else{return false;}
}

void GlobalPruning( vector<FAM>& BayesN, vector<NODE>& Tree)
{
  FAM famtemp;
  vector<int> TempCO;
  vector<vector<int> > CoOccur;
  
  for(int i=0;i<Tree[0].gtype.size();i++)
  {
    TempCO.push_back(0);
    famtemp.pi.push_back(false);
  }
  famtemp.potpi=famtemp.pi;
  famtemp.npotpa=false;
  famtemp.npa=0;
  famtemp.nmut=0.0;
  famtemp.score=0.0;
  for(int i=0;i<Tree[0].gtype.size();i++)
  {
    
    BayesN.push_back(famtemp);
    CoOccur.push_back(TempCO);
  }
  
  vector<double> freq(BayesN.size());
  vector<int> Order;
  for(int i=0;i<BayesN.size();i++)
  {
    Order.push_back(i);
    freq[i]=0.0;
    for(int j=0;j<Tree.size();j++)
    {
      if(Tree[j].gtype[i]==true){freq[i]+=1.0;}
      for(int k=i+1;k<BayesN.size();k++)
      {
        if((Tree[j].gtype[i]==true) && (Tree[j].gtype[k]==true))
        {
          CoOccur[i][k]++;
          CoOccur[k][i]++;
        }
      }
      
    }
    
  }
  int c1=0;
  for(int i=0;i<BayesN.size()-1;i++)
  {
    for(int j=i+1;j<BayesN.size();j++)
    {
      if(freq[i]<freq[j])
      {
        int bufog=Order[i];
        Order[i]=Order[j];
        Order[j]=bufog;
        double ftem=freq[i];
        freq[i]=freq[j];
        freq[j]=ftem;
      }
    }
  }
  for(int i=0;i<BayesN.size();i++)
  {
    BayesN[Order[i]].order=i;
    
  }
  
  
  
  for(int i=0;i<BayesN.size()-1;i++)
  {
    
    bool indep=true;
    for(int j=i+1;j<BayesN.size();j++)
    {
      double olap=(double) CoOccur[i][j];
      double nsamp= (double) Tree.size();
      int a=BayesN[i].order;
      int b=BayesN[j].order;
      
      if(MScore(olap,freq[a],freq[b],nsamp)==true)
      {
        BayesN[i].potpi[j]=true;
        BayesN[j].potpi[i]=true;
        BayesN[i].npotpa=true;
        BayesN[j].npotpa=true;
        c1++;
        
      }
      
    }
  }

  InitEnt( BayesN, Tree);
}


void SanityCheckTree(vector<NODE>& Tree)
{
  if(Tree[0].child.size()!=1)
  {
    Rcout<<"PROBLEM in INITTREE.cpp. Zero has "<<Tree[0].child.size()<<" children";
    throw 2;
  }
  for(int i=1;i<Tree.size();i++)
  {
    int dad=Tree[i].dad;
    if((Tree[dad].child[0]!=i) && (Tree[dad].child[1]!=i))
    {
      Rcout<<"PROBLEM in INITTREE.cpp";
      throw 2;
    }
    
  }
}


void InitTree( vector<NODE>& Tree, vector<FAM>& BN)
{
  
  vector<int> Orphans;
  
  int Ntax=int(Tree.size());
  
  for(int i=1;i<Ntax;i++)
  {
    Orphans.push_back(i);
  }
  
  vector<int> Nmut;
  for(int i=0;i<BN.size();i++)
  {
    Nmut.push_back(BN[i].nmut);
  }
  
  for(int i=Ntax;i<2*Ntax-2;i++)
  {
    double y = (double) rand()/RAND_MAX;
    int x = int(y*(Orphans.size() - 1));
    int rchld=Orphans[x];
    
    Orphans.erase(Orphans.begin() +x);
    
    
    y = (double) rand()/RAND_MAX;
    
    x= int(y*(Orphans.size() - 1));
    int lchld =Orphans[x];
    
    
    NODE n;
    n.child.push_back(rchld);
    n.child.push_back(lchld);
    
    for(int k=0;k<Tree[0].gtype.size();k++)									/*	each time a new internal vertex is found			*/
    {
      bool gbit= Tree[rchld].gtype[k]*Tree[lchld].gtype[k];
      
      if(Nmut[k]>Ntax-2)
      {
        gbit=false;
      }
      n.gtype.push_back(gbit);
      if(gbit==true){Nmut[k]++;}
    }
    
    n.NDes=Tree[rchld].NDes + Tree[lchld].NDes;
    Tree.push_back(n);
    
    
    
    Tree[rchld].dad=i;
    Tree[lchld].dad=i;
    
    
    Orphans.erase(Orphans.begin() +x);
    Orphans.push_back(Tree.size()-1);
    
  }
  if(Orphans.size()!=1){Rcout<<"Problem initializing tree. "<<Orphans.size()<<'\n';}else
  {
    int rchld=Tree.size()-1;
    
    Tree[0].child.push_back(rchld);
    Tree[rchld].dad=0;
    
  }
  for(int i=1;i<Tree.size();i++)
  {
    for(int j=0;j<Tree[i].gtype.size();j++)
    {
      if(Tree[i].gtype[j]!=Tree[Tree[i].dad].gtype[j])
      {
        Tree[i].mut.push_back(j);
      }
    }
  }
  SanityCheckTree( Tree);
  
  
  
}

void Delete_Col(vector< vector<bool> >& Data, int col_no)
{
  for(int i =0;i<Data.size();i++)
  {
    vector<bool> temp =Data[i];
    temp.erase(temp.begin() + col_no);
    Data[i].swap(temp);
  }
  
}
void EliminateRedundantCharacters(vector< vector<bool> >& Data, vector<vector<int> >& char_label)//discard repeated columns
{
  vector<bool> RedChar;
  for (int i=0; i<Data[0].size(); i++)
  {
    RedChar.push_back(false);
  }
  
  for(int ch1=0;ch1<Data[0].size();ch1++)
  {
    
    if(RedChar[ch1]==false)
    {
      vector<int> templ;
      templ.push_back(ch1);
      for(int k=ch1+1;k<Data[0].size();k++)
      {
        
        
        if(RedChar[k]==false)
        {
          RedChar[k]=true;
          for(int i =1; i<Data.size();i++)
          {
            
            if(Data[i][ch1]!=Data[i][k])
            {
              
              RedChar[k]=false;
              
            }
          }
          if(RedChar[k]==true)
          {
            templ.push_back(k);
            //Rcout<<ch1<<" and "<<k<<'\n';
          }
          
        }
      }
      char_label.push_back(templ);
    }
  }
  int j=0;
  
  while(Data[0].size()>char_label.size())
  {
    
    if(RedChar[j]==true)
    {
      //Rcout<<"REDCHAR "<<j<<'\n';
      Delete_Col(Data,j);
      RedChar.erase(RedChar.begin() +j);
      j--;
    }
    j++;
  }
  Rcout<<Data[0].size()<<'\t'<<char_label.size()<<'\n';
}


bool IsInteger(string& inp, int& Num)
{
  Num=0;
  for(int i=0;i<inp.length();i++)
  {
    if( (int(inp[i])-int('0')>=0) && (int(inp[i])-int('0')<=9))
    {
      Num =Num*10 + int(inp[i])-int('0');
    }
    else
    {
      return false;
    }
  }
  return true;
}

/*	Read the data file, store taxon and gene names as a vector and partially initialize the variable "Tree".					*/

void ParseDataMatrix(vector<vector<bool> >& Data, vector<string>& TaxonLabels, vector<string>& GeneLabels, string& dname)
{
  ifstream DataFile(dname.c_str());
  if(!DataFile){throw 0;}else{Rcout<<"Reading "<<dname<<'\n';}
  int NTaxa, NGenes;
  string inp;
  DataFile>>inp;
  if(!IsInteger(inp,NTaxa))
  {
    cerr<<"\nNon integral entry "<<inp<<" for # taxa in file "<<dname<<'\n';
    throw 1;
  }
  DataFile>>NGenes;
  if(!IsInteger(inp,NTaxa))
  {
    cerr<<"\nNon integral entry "<<inp<<" for # genes/characters in file "<<dname<<'\n';
    throw 1;
  }
  
  for(int i=0;i<NGenes;i++)				//	Storing gene names, so that we can work with integers from
  {										//	now on and have peace of mind. Gene names are weird !
    DataFile>>inp;
    GeneLabels.push_back(inp);
  }
  vector<int> hist(11);
  for(int i=0;i<11;i++)
  {
    hist[i]=0;
  }
  for(int i=0;i<NTaxa;i++)
  {
    DataFile>>inp;
    TaxonLabels.push_back(inp);	//	Storing taxon or sample names. Good riddance! These are even weirder than gene names.
    vector<bool> v;
    int nm=0;
    for(int j=0;j<NGenes;j++)
    {
      string state;
      DataFile>>state;
      if(state=="0"){v.push_back(false);}
      else
      {
        if(state=="1")
        {
          v.push_back(true);
          nm++;
        }
        else
        {
          cerr<<"\nProblem parsing the data matrix entry for taxon (row) # "<<i<<" and gene (column) # "<<j<<'\n'<<"Either the entries are not 0/1 or Num Genes is incorrectly specified in the header as "<<NGenes<<'\n';
          throw 1;
        }
      }
    }
    if(nm<10)
    {
      hist[nm]++;
    }else{hist[10]++;}
    Data.push_back(v);
  }
  for(int i=0;i<11;i++)
  {
    Rcout<<i<<'\t'<<hist[i]<<endl;
  }
  Rcout<<"\nTaxa initialized\n"<<Data.size()<<" taxa and "<<NGenes<<" genes/characters found."<<'\n';
  DataFile.close();
  
}


void DFS(vector<bool>& conf, vector<bool>& seen, int& label, vector< vector< double> >& cdata, bool& insta)
{
  for(int nb=0;nb<conf.size();nb++)
  {
    if(conf[nb]==false)
    {
      conf[nb]=true;
      label+=int(pow(2.0,nb));
      double nmin1, nmin2;
      if(cdata[0][label]>0.0){nmin1=cdata[0][label]/(cdata[0][label]+cdata[1][label]);}else{nmin1=0.0;}
      if(cdata[0][label-int(pow(2.0,nb))]>0.0){nmin2=cdata[0][label-int(pow(2.0,nb))]/(cdata[0][label-int(pow(2.0,nb))]+cdata[1][label-int(pow(2.0,nb))]);}else{nmin2=0.0;}
      
      if(nmin1<=nmin2 && insta==true)
      {
        if(seen[label]==false){DFS(conf,seen,label, cdata, insta);}
        
      }else
      {
        insta=false;
      }
      seen[label]=true;
      conf[nb]=false;
      label-=int(pow(2.0,nb));
    }
    
  }
}



/*	Greedy K2 heuristic for updating the parents of gene "G" in "FamG" given "Tree"				*
*	cf. Cooper & Herskovits (1992).																*/

void K2Search( FAM& FamG, vector<FAM>& DAG, vector<NODE>& Tree, int g, char op)
{
  for(int j=0;j<DAG.size();j++)
  {
    FamG.pi[j]=false;
  }
  FamG.npa=0;
  double theta = FamG.nmut/Tree.size();
  FamG.param.erase(FamG.param.begin(),FamG.param.end());
  
  FamG.param.push_back(1.0 -theta);
  FamG.param.push_back(theta);
  FamG.score=0.0;
  bool repeat=true;
  while(repeat==true)
  {
    double besc=FamG.score;
    int add=0;
    double nvert=(double) Tree.size();
    
    int npi=FamG.npa;
    vector<double> famparam(int(pow(2.0,npi+2))), cst(2);
    repeat=false;
    for(int i=0;i<DAG.size();i++)
    {
      vector< vector<double> > cdata(2);
      
      if((DAG[g].potpi[i]==1 && DAG[i].potpi[g]==1)|| op=='n')
      {
        if( DAG[g].order>DAG[i].order && FamG.pi[i]==0 )
        {
          for(int m=0;m<int(pow(2.0,npi+1));m++)
          {
            cdata[0].push_back(0);
            cdata[1].push_back(0);
          }
          cst[0]=0;
          cst[1]=0;
          
          nvert=0.0;
          FamG.pi[i]=1;
          for(int j=0;j<Tree.size();j++)
          {
            int indxg= Tree[j].gtype[g];
            int indxp=0;
            int indxi=0;
            int par=0;
            for(int k=0;k<Tree[j].gtype.size();k++)
            {
              if(FamG.pi[k]==1)
              {
                
                indxp+=int(pow(2.0,par))* Tree[j].gtype[k];
                par++;
              }
            }
            cdata[indxg][indxp]++;
            cst[indxg]++;
            nvert++;
          }
          FamG.pi[i]=0;
          double Minf=0.0;
          vector<bool> seen;
          for(int l=0;l<int(pow(2.0,npi+1));l++)
          {
            if(cdata[0][l]>0.0){Minf+=cdata[0][l]*log((nvert*cdata[0][l])/((cdata[0][l]+cdata[1][l])*(cst[0])));}
            if(cdata[1][l]>0.0){Minf+=cdata[1][l]*log((nvert*cdata[1][l])/((cdata[0][l]+cdata[1][l])*(cst[1])));}
            seen.push_back(false);
          }
          vector<bool> conf;
          for(int l=0;l<npi+1;l++)
          {
            conf.push_back(false);
          }
          int label=0;
          bool insta=true;
          if(op=='y'){DFS(conf,seen,label, cdata, insta);}
          
          if( insta==true && Minf-pow(2.0,npi)*log(nvert)>besc)
          {
            besc=Minf-pow(2.0,npi)*log(nvert);
            add=i;
            for(int pin=0;pin<int(pow(2.0,npi+1));pin++)
            {
              if(cdata[0][pin]+cdata[1][pin] >0.0)
              {
                famparam[pin]=(cdata[0][pin]+ cst[0]/nvert)/(cdata[0][pin]+cdata[1][pin]+1.0);
                famparam[pin+int(pow(2.0,npi+1))]=(cdata[1][pin] + cst[1]/nvert)/((cdata[0][pin]+cdata[1][pin])+1.0);
              }
              else
              {
                famparam[pin]=1.0;
                famparam[pin+int(pow(2.0,npi+1))]=0.0;
              }
            }
          }
        }
      }
      
    }
    if((besc>FamG.score))
    {
      FamG.pi[add]=1;
      FamG.score=besc;
      FamG.npa++;
      FamG.param=famparam;
      repeat=true;
    }
    
  }
  
}


/*	Perform local search operations on "DAG" given "Tree" for structure learning.				*/

void SearchDAGs( vector<FAM>& DAG, vector<NODE>& Tree, char op )
{
  bool opt=false;
  
  vector<int> Ord(DAG.size());
  for(int i=0;i<DAG.size();i++)
  {
    Ord[DAG[i].order]=i;
  }
  
  
  
  for(int i=0;i<DAG.size();i++)
  {
    if(DAG[i].npotpa==true||op=='n'){K2Search(DAG[i],DAG,Tree,i,op);}
  }
  int scount=0;
  while(opt==false)
  {
    opt=true;
    for(int i=0;i<DAG.size()-1;i++)
    {
      if((DAG[Ord[i]].potpi[Ord[i+1]]==true) && (DAG[Ord[i+1]].potpi[Ord[i]]==true))
      {
        int tempord=Ord[i];
        DAG[Ord[i]].order+=1;
        DAG[Ord[i+1]].order-=1;
        Ord[i]=Ord[i+1];
        Ord[i+1]=tempord;
        FAM fg1, fg2;
        fg1=DAG[Ord[i]];
        fg2=DAG[Ord[i+1]];
        
        //if(Ord[DAG[i].order]!=i){Rcout<<"problem in reordering genes for gene"<<i<<'\t';}
        
        if(DAG[Ord[i]].pi[Ord[i+1]]==1)
        {
          K2Search(fg1,DAG,Tree,Ord[i],op);
        }
        K2Search(fg2,DAG,Tree,Ord[i+1],op);
        if((fg1.score+fg2.score)>(DAG[Ord[i]].score+DAG[Ord[i+1]].score))
        {
          scount ++;
          
          //Rcout<<"\n\n\n Reordering accepted ! Genes"<<Ord[i]<<" AND "<<Ord[i+1]<<'\n'<<"score increases by "<<fg1.score+fg2.score -(DAG[Ord[i]].score+DAG[Ord[i+1]].score)<<'\t';
          DAG[Ord[i]]=fg1;
          DAG[Ord[i+1]]=fg2;
          opt=false;
          
          double psc=0.0;
          for(int i=0;i<DAG.size();i++)
          {
            psc+=DAG[i].score;
          }
          // Rcout<<" to "<<psc<<" reordered "<<scount<<" times \n";
        }else
        {
          int tempord=Ord[i];
          DAG[Ord[i]].order+=1;
          DAG[Ord[i+1]].order-=1;
          Ord[i]=Ord[i+1];
          Ord[i+1]=tempord;
        }
      }
    }
  }
  
}

void RelabelOutlier(vector<NODE>& Tree,  int gn, int pdad, bool& tnew)
{
  int min=Tree.size();
  int lab=0;
  for(int i=Tree.size()/2+1;i<Tree.size();i++)
  {
    if(Tree[i].gtype[gn]==false &&(Tree[Tree[i].child[0]].gtype[gn]*Tree[Tree[i].child[1]].gtype[gn]==true))
    {
      if(min>Tree[i].NDes)
      {
        min=Tree[i].NDes;
        lab=i;
      }
    }
  }
  Tree[lab].gtype[gn]=true;
  
}

void swapnodes(vector<NODE>& Tree, vector<FAM>& DAG, int& i, int& gramps, int& spart, bool& change)
{
  int node1=0;
  int sib=Tree[Tree[i].dad].child[1];
  if(Tree[Tree[i].dad].child[1]==i)
  {
    node1=1;
    sib=Tree[Tree[i].dad].child[0];
  }
  int node2=0;
  if(Tree[gramps].child[0]==Tree[i].dad){node2=1;}
  
  spart = Tree[gramps].child[node2];
  
  change=false;
  for(int j=0;j<Tree[i].gtype.size();j++)
  {
    bool tempnew=Tree[spart].gtype[j]*Tree[sib].gtype[j];
    bool tempold=Tree[Tree[i].dad].gtype[j];
    bool CSOld = Tree[i].gtype[j]*Tree[sib].gtype[j];
    if((CSOld!=tempold))
    {
      if((DAG[j].nmut==Tree.size()/2) && tempold==false)
      {
        tempnew=false;
      }
      else
      {
        Rcout<<"incorrect label for node "<<i<<" sib "<<sib<<" dad "<<Tree[i].dad<<" gramps "<<gramps<<" spart "<<spart<<endl;
        Rcout<<" Nmut "<<DAG[j].nmut<<endl;
        Rcout<<"Dad's children "<<Tree[Tree[i].dad].child[0]<<" and "<<Tree[Tree[i].dad].child[1]<<endl;
        Rcout<<"gramps children "<<Tree[gramps].child[0]<<" and "<<Tree[gramps].child[1]<<endl;
        Rcout<<"spart's father "<<Tree[spart].dad<<endl;
        Rcout<<"sib's father "<<Tree[sib].dad<<endl;
        Rcout<<j<<" gene, node ="<<Tree[i].gtype[j]<<" sib = "<<Tree[sib].gtype[j]<<" dad = "<<Tree[Tree[i].dad].gtype[j]<<" gramps = "<<Tree[gramps].gtype[j]<<" spart = "<<Tree[spart].gtype[j]<<endl;
        throw 2;
      }
    }
    double theta=0.0;
    if(tempold!=tempnew)
    {
      change=true;
      if(tempold==true)
      {
        if(DAG[j].nmut==Tree.size()/2)
        {
          int canj=j;
          RelabelOutlier(Tree,canj,Tree[i].dad, tempnew);
        }else
        {
          DAG[j].nmut-=1.0;
        }
      }else
      {
        if(DAG[j].nmut<Tree.size()/2){DAG[j].nmut+=1.0;}else{tempnew=false;}
      }
      theta = DAG[j].nmut/Tree.size();
      if(theta>0.0)
      {
        DAG[j].ent=-DAG[j].nmut*log(theta) - (Tree.size()-DAG[j].nmut)*log(1.0-theta);
      }else
      {
        DAG[j].ent=0.0;
      }
      Tree[Tree[i].dad].gtype[j]=tempnew;
      for(int k=0;k<DAG.size();k++)
      {
        if(k<j && DAG[k].potpi[j]==true)
        {
          if(TestMI(k, j,Tree )==true)
          {
            
            DAG[j].potpi[k]=true;
          }
          else
          {
            DAG[j].potpi[k]=false;
          }
          
          
        }
        if(j<k && DAG[j].potpi[k]==true)
        {
          if(TestMI(j, k,Tree )==true)
          {
            
            DAG[k].potpi[j]=true;
          }
          else
          {
            DAG[k].potpi[j]=false;
          }
          
          
        }
        
      }
      
    }
  }
  
  Tree[gramps].child[node2]=i;
  Tree[Tree[i].dad].child[node1]=spart;
  
  Tree[spart].dad=Tree[i].dad;
  Tree[i].dad=gramps;
}

void SearchTrees( vector<FAM>& DAG, vector<NODE>& Tree)
{
  bool opt =false;
  double bscore=0.0;
  int deg=0;
  
  SearchDAGs(DAG, Tree,'y');
  for(int i=0; i<DAG.size();i++)
  {
    bscore+=DAG[i].score-DAG[i].ent;
    deg+=DAG[i].npa;
  }
  while(opt==false)
  {
    
    opt=true;
    Rcout<<" Tree size = "<<Tree.size()<<" present score "<<bscore<<'\n';
    Rcout<<"degree "<<deg<<endl;
    
    for(int i=1; i<Tree.size();i++)
    {
      
      SanityCheckTree(Tree);
      int spart=0;
      int cand=i;
      double nscore=0.0;
      vector<int> tempord;
      int gramps= Tree[Tree[i].dad].dad;
      bool change=false;
      if(gramps!=0)
      {
        swapnodes(Tree,DAG,cand,gramps,spart, change);
        if(change==true)
        {
          for(int j=0; j<DAG.size();j++)
          {
            tempord.push_back(DAG[j].order);
          }
          
          SearchDAGs(DAG,Tree,'y');
          
          for(int j=0; j<DAG.size();j++)
          {
            nscore+=DAG[j].score-DAG[j].ent;
          }
          if(nscore>bscore)
          {
            opt=false;
            bscore=nscore;
            deg=0;
            for(int j=0; j<DAG.size();j++)
            {
              deg+=DAG[j].npa;
            }
          }else
          {
            swapnodes(Tree,DAG,spart,gramps,cand, change);
            for(int j=0; j<DAG.size();j++)
            {
              DAG[j].order=tempord[j];
            }
          }
        }
      }
    }
  }
  Rcout<<"\n Search trees done \n";
  double finsc=0.0;
  
  SearchDAGs(DAG, Tree,'y');
  int nedges=.0;
  for(int j=0; j<DAG.size();j++)
  {
    finsc+=DAG[j].score-DAG[j].ent;
    nedges+=DAG[j].npa;
    
  }
  Rcout<<nedges<<" edges and log-likl = "<<finsc<<'\n';
  
  
}

List StructLearn(vector< vector<bool> >& Data, vector<string>& GeneLabels, vector<SimVar>& SimDAG, int& NTrees, double& pthres)
{
  
  vector<vector<int> > char_label;
  // ParseDataMatrix(Data, TaxonLabels, GeneLabels, dname);
  
  EliminateRedundantCharacters(Data,char_label);
  
  vector<NODE> Emtree;
  vector<FAM> BayesN;
  
  GenerateReplicate(Emtree,Data);
  GlobalPruning(BayesN,Emtree);
  
  List result = List::create();
  result["num_edges"] = BayesN.size()*(BayesN.size()-1)/2.0;
  result["num_genes"] = BayesN.size();
  
  int num_unpruned_edges = 0;
  int num_noparent_after_global_pruning = 0;
  
  for(size_t i=0;i<BayesN.size()-1;i++)
  {
    for(size_t j=i+1;j<BayesN.size();j++)
    {
      if(BayesN[i].potpi[j] && BayesN[j].potpi[i] && BayesN[i].npotpa && BayesN[j].npotpa)
      {
        num_unpruned_edges++;
      }
    }
    if(BayesN[i].npotpa==0) num_noparent_after_global_pruning++;
  }
  
  result["num_unpruned_edges"] = num_unpruned_edges;
  result["num_noparent_after_global_pruning"] = num_noparent_after_global_pruning;
  /*	Initiate a subroutine for randomly reseeding the search procedure. Each tree serves as	*
  *	a starting point for a local search heuristic through tree space						*/
  
  
  double bnscore=0.0;
  
  for(int Seed=0;Seed<NTrees;Seed++)
  {
    /*	Define the basic network object		*/
    vector<FAM> DAG;
    
    vector<NODE> Phy;
    GenerateReplicate(Phy,Data);
    GlobalPruning(DAG,Phy);
    /*	Initialize the missing data values in the                          *
    *	variable "Phy", cf. InitTree.cpp and perform local pruning													*/
    
    InitTree(Phy,DAG);
    Rcout<<"tree initialized\n Seed "<<Seed+1<<" out of "<<NTrees<<endl;
    LocalPruning(DAG,Phy);
    
    double dagscore=0.0;
    
    
    
    /*	Perform local search operations in tree space for imputing missing data values		*
    *	in variable "Tree" using the MPE (Most Probable Explanation) criterion. 			*/
    
    SearchTrees(DAG,Phy);
    
    for(int i=0;i<DAG.size();i++)
    {
      dagscore+=DAG[i].score-DAG[i].ent;
    }
    if(dagscore>bnscore || Seed==0)
    {
      for(int i=0;i<BayesN.size();i++)
      {
        BayesN[i]=DAG[i];
        bnscore=dagscore;
      }
      
    }
    
  }
  Rcout<<"writing Pathfile"<<endl;
  ParseDAG(BayesN, SimDAG, GeneLabels,char_label);
  result["DAG"] = WriteLandscape(SimDAG, Emtree, GeneLabels, char_label,pthres);
  return result;
}

List BootstrapAnalysis(vector< vector<bool> >& Data, vector<SimVar>& SimDAG, int& NRep)
{
  /* Edges present in Tree+OBSDAG, OBSDAG and SimDAG resp. */
  vector< vector<int> > netrel1, netrel2, netrel3;
  vector< vector<string> > EdList;
  InitNetRel(netrel1, netrel2, netrel3, EdList, SimDAG);
  
  /* Bootstrap probabilities for the top ordered genes for Tree+OBSDAG and OBSDAG (not to be confused with the ordering of mutations */
  vector<vector<double> > landrel1(SimDAG.size());
  vector<vector<double> > landrel2(SimDAG.size());
  
  /* Bootstrap probabilities for each pair of genes sharing an edge, 2 single mutatant states and the double mutant genotype using Tree+OBSDAG */
  vector< vector<double> > edrel(3*EdList.size());
  
  
  for(int Brep=1; Brep<NRep+1;Brep++)
  {
    vector<NODE> Tree;
    vector<FAM> BayesNP;
    
    SimulateDAG(Tree, SimDAG, Data);
    
    Rcout<<"\nGenerated bootstrap replicate "<<Brep<<" out of "<<NRep<<'\n';

    GlobalPruning(BayesNP,Tree);
    SearchDAGs(BayesNP,Tree,'n');
    for(int i=0;i<netrel2.size();i++)
    {
      for(int j=0;j<netrel2.size();j++)
      {
        netrel2[i][j]+=int(BayesNP[i].pi[j])+int(BayesNP[j].pi[i]);
      }
    }
    LandBoot(SimDAG, landrel1, BayesNP);
    
    
    double bnscore=0.0;
    
    
    /*	Parse each input tree and assign the missing data values in the						*
    *	variable "Tree", cf. InitTree.cpp.													*/
    
    vector<FAM> DAG(BayesNP.size());					/*	Define the basic network object	*/
    InitTree(Tree, BayesNP);
    Rcout<<"tree initialized\n";
    InitDAG(DAG, BayesNP);
    LocalPruning(DAG,Tree);
    
    double dagscore=0.0;
    
    /*	Perform local search operations in tree space for imputing missing data values		*
    *	in variable "Tree" using the MPE (Most Probable Explanation) criterion. 			*/
    
    SearchTrees(DAG,Tree);
    Rcout<<Brep<<" out of "<<NRep<<" bootstrap replicates analyzed\n";
    for(int i=0;i<netrel1.size();i++)
    {
      for(int j=0;j<netrel1.size();j++)
      {
        netrel1[i][j]+=int(DAG[i].pi[j])+int(DAG[j].pi[i]);
      }
    }
    LandBoot(SimDAG, landrel2, DAG);
    EdBoot(EdList,SimDAG,edrel, DAG);
  }
  // TPFPAnalysis(netrel1,netrel2,netrel3,SimDAG,NRep,job);
  List result = TPFPAnalysis(netrel1,netrel2,netrel3,SimDAG,NRep);
  List probs = ProbBootstrapAnalysis(landrel1,landrel2,SimDAG);
  result["OBS_Probabilities"] = probs["OBS_Probabilities"];
  result["Tree_OBS_Probabilities"] = probs["Tree_OBS_Probabilities"];
  result["EdgeProbabilities"] = ProbEdBoot(EdList, edrel);
  return result;
}