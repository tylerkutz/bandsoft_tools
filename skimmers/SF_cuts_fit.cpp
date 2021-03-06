#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "reader.h"
#include "bank.h"
#include "clas12fiducial.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;

double Epcal_max=1.2;

double FF(double E, double sf1, double sf2, double sf3, double sf4){
  return sf1 * (sf2 + (sf3/E) + (sf4/(E*E)) ); 
}

int main(int argc, char** argv) {

  if( argc != 4 ){
    cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [outPutTextFile] [inputFile] \n\n";
    return -1;
  }
  
  ofstream cutFile;
  cutFile.open(argv[2]);
  

  //Define Variables
  int Runno;
  double Ebeam, gated_charge, livetime, starttime, current;
  clashit * eHit = new clashit;
  char temp[100];

 //Creat input tree
  TFile * inFile = new TFile(argv[3]);
  TTree * inTree = (TTree*)inFile->Get("electrons");
  // 	Event branches:
  inTree->SetBranchAddress("Runno"               ,&Runno                 );
  inTree->SetBranchAddress("Ebeam"               ,&Ebeam                 );
  inTree->SetBranchAddress("gated_charge"        ,&gated_charge          );
  inTree->SetBranchAddress("livetime"            ,&livetime              );
  inTree->SetBranchAddress("starttime"           ,&starttime             );
  inTree->SetBranchAddress("current"             ,&current               );
  //	Electron branches:
  inTree->SetBranchAddress("eHit"		 ,&eHit			 );

  
  // Create output tree
  TFile * outFile = new TFile(argv[1],"RECREATE");
  /*
  TTree * outTree = new TTree("electrons","CLAS Electrons");
  // 	Event branches:
  outTree->Branch("Runno"		,&Runno			);
  outTree->Branch("Ebeam"		,&Ebeam			);
  outTree->Branch("gated_charge"	,&gated_charge		);
  outTree->Branch("livetime"	        ,&livetime		);
  outTree->Branch("starttime"	        ,&starttime		);
  outTree->Branch("current"	        ,&current		);
  //	Electron branches:
  outTree->Branch("eHit"		,&eHit			);
  */

  //Creating TH2D from tree
  vector<TH2*> hist_list;
  TH2D * h2[6];
  for(int i = 0; i<6; i++){
    sprintf(temp,"TH2D_sec%d",i);
    h2[i] = new TH2D(temp,"SF_v_Epcal;Epcal;SF;Events",250,0,1.7,150,0.05,0.40);
    hist_list.push_back(h2[i]);
  } 
  int fin = inTree->GetEntries();
  for(int i = 0; i < fin; i++){
    eHit->Clear();
    inTree->GetEntry(i);
    //Display completed  
    if((i%100000) == 0){
      cerr << (i*100.)/fin <<"% complete \n";
    }
    if(eHit->getPID() != 11){continue;}
    if(eHit->getCharge() != -1){continue;}
    if(eHit->getV() < 14){ continue; }
    if(eHit->getW() < 14){ continue; }
    h2[eHit->getSector() - 1]->Fill(eHit->getEpcal(),eHit->getEoP());
  }

  //Make Tgraphs from TH2D
  int ctr = 0;
  vector<TGraph*> graph_list;
  TGraph * tg_mu[6];
  TGraph * tg_sigma[6];
  for(int i = 0; i<6; i++){
    //Make the graphs
    sprintf(temp,"TGraph_Mean_sec%d",i);
    tg_mu[i] = new TGraph();
    tg_mu[i]->SetName(temp);
    graph_list.push_back(tg_mu[i]);
    sprintf(temp,"TGraph_Sigma_sec%d",i);
    tg_sigma[i] = new TGraph();
    tg_sigma[i]->SetName(temp);
    graph_list.push_back(tg_sigma[i]);
    //Now project the histogram    
    for(int j = 0; j < h2[i]->GetXaxis()->GetNbins(); j++){
      double x = h2[i]->GetXaxis()->GetBinCenter(j+1);
      ctr++;
      sprintf(temp,"Proj_num%d",ctr);
      TH1D * proj = h2[i]->ProjectionY(temp,j+1,j+2);
      //Now preform a guassian fit
      double A = proj->GetMaximum();
      double mode = proj->GetXaxis()->GetBinCenter(proj->GetMaximumBin());
      double mu = proj->GetMean();
      double sigma = proj->GetStdDev();
      TF1 * gFit = new TF1("GausFit","gaus",0,0.4);
      gFit->SetParameter(0,A);
      gFit->FixParameter(1,mode);
      gFit->SetParameter(2,sigma);
      TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",0,0.4);
      //Now determine if you want to add the points
      int check = gPoint;
      if(check == 0){
	mu = gPoint->Parameter(1);
	sigma = gPoint->Parameter(2);
	if(proj->GetEntries()>100){
	  if((mu>0.0) && (mu<0.4)){
	    //Write out x vs y and x vs s
	    tg_mu[i]->SetPoint(tg_mu[i]->GetN(),x,mu);
	    tg_sigma[i]->SetPoint(tg_sigma[i]->GetN(),x,sigma);
	  }
	}
      }
    }
  }
  
  cout<<"Now fit the SF\n";
  //The first two are the two functions I will fit to
  TF1 * f1_mu = new TF1("f1_mu",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]); },0.06,1.6,4);
  f1_mu->SetLineColor(2);
  TF1 * f1_sigma = new TF1("f1_sigma",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]); },0.06,1.6,4);
  f1_sigma->SetLineColor(2);
  
  //These are the functions that will show the maximum and minimum defined by the fit
  TF1 * f1_max = new TF1("f1_max",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]) + 5*FF(x[0],p[4],p[5],p[6],p[7]); },0.06,1.6,8);
  f1_max->SetLineColor(2);
  TF1 * f1_min = new TF1("f1_min",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]) - 5*FF(x[0],p[4],p[5],p[6],p[7]); },0.06,1.6,8);
  f1_min->SetLineColor(2);

  //These functions will be used to compare with the RGA values
  /*
  TF1 * f1_mu_rga = new TF1("f1_mu_rga",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]); },0.06,1.6,4);
  f1_mu_rga->SetLineColor(3);
  TF1 * f1_sigma_rga = new TF1("f1_sigma_rga",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]); },0.06,1.6,4);
  f1_sigma_rga->SetLineColor(3);
  TF1 * f1_max_rga = new TF1("f1_max_rga",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]) + 5*FF(x[0],p[4],p[5],p[6],p[7]); },0.06,1.6,8);
  f1_max_rga->SetLineColor(3);
  TF1 * f1_min_rga = new TF1("f1_min_rga",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2],p[3]) - 5*FF(x[0],p[4],p[5],p[6],p[7]); },0.06,1.6,8);
  f1_min_rga->SetLineColor(3);
  */
  //Give some estimates
  f1_mu->SetParameter(0,0.25);
  f1_mu->SetParameter(1,1.00);
  f1_mu->SetParameter(2,-0.02);
  f1_mu->SetParameter(3,-0.0001);
  f1_sigma->SetParameter(0,0.018);
  f1_sigma->SetParameter(1,1.00);
  f1_sigma->SetParameter(2,0.0);
  f1_sigma->SetParameter(3,0.0);

  for(int i = 0; i<6; i++){
    sprintf(temp,"Canvas_sec%d",i);
    TCanvas * c1 = new TCanvas(temp,temp,1200,1000);
    c1->cd();
    h2[i]->Draw("colz");
    tg_mu[i]->Draw("SAME");
    
    //Start with RGA
    /*
    f1_mu_rga->SetParameter(0,SF1[i]);
    f1_mu_rga->SetParameter(1,1.00);
    f1_mu_rga->SetParameter(2,SF3[i]);
    f1_mu_rga->SetParameter(3,SF4[i]);
    f1_sigma_rga->SetParameter(0,SFs1[i]);
    f1_sigma_rga->SetParameter(1,1.00);
    f1_sigma_rga->SetParameter(2,0.00);
    f1_sigma_rga->SetParameter(3,0.00);

    f1_max_rga->SetParameter(0,SF1[i]);
    f1_max_rga->SetParameter(1,1.00);
    f1_max_rga->SetParameter(2,SF3[i]);
    f1_max_rga->SetParameter(3,SF4[i]);
    f1_max_rga->SetParameter(4,SFs1[i]);
    f1_max_rga->SetParameter(5,1.00);
    f1_max_rga->SetParameter(6,0.00);
    f1_max_rga->SetParameter(7,0.00);

    f1_min_rga->SetParameter(0,SF1[i]);
    f1_min_rga->SetParameter(1,1.00);
    f1_min_rga->SetParameter(2,SF3[i]);
    f1_min_rga->SetParameter(3,SF4[i]);
    f1_min_rga->SetParameter(4,SFs1[i]);
    f1_min_rga->SetParameter(5,1.00);
    f1_min_rga->SetParameter(6,0.00);
    f1_min_rga->SetParameter(7,0.00);
    
    f1_mu_rga->Draw("SAME");
    f1_max_rga->Draw("SAME");
    f1_min_rga->Draw("SAME");
    */

    //Now do fit
    TFitResultPtr muPoint = tg_mu[i]->Fit(f1_mu,"qeSrn","",0.06,Epcal_max);
    TFitResultPtr sigmaPoint = tg_sigma[i]->Fit(f1_sigma,"qeSrn","",0.2,Epcal_max);
    f1_mu->Draw("SAME");

    f1_max->SetParameter(0,muPoint->Parameter(0));
    f1_max->SetParameter(1,muPoint->Parameter(1));
    f1_max->SetParameter(2,muPoint->Parameter(2));
    f1_max->SetParameter(3,muPoint->Parameter(3));
    f1_max->SetParameter(4,sigmaPoint->Parameter(0));
    f1_max->SetParameter(5,sigmaPoint->Parameter(1));
    f1_max->SetParameter(6,sigmaPoint->Parameter(2));
    f1_max->SetParameter(7,sigmaPoint->Parameter(3));
    f1_max->Draw("SAME");

    f1_min->SetParameter(0,muPoint->Parameter(0));
    f1_min->SetParameter(1,muPoint->Parameter(1));
    f1_min->SetParameter(2,muPoint->Parameter(2));
    f1_min->SetParameter(3,muPoint->Parameter(3));
    f1_min->SetParameter(4,sigmaPoint->Parameter(0));
    f1_min->SetParameter(5,sigmaPoint->Parameter(1));
    f1_min->SetParameter(6,sigmaPoint->Parameter(2));
    f1_min->SetParameter(7,sigmaPoint->Parameter(3));
    f1_min->Draw("SAME");

    sprintf(temp,"Canvas_sec%d.pdf",i);
    c1->SaveAs(temp);

    //Quickly write out
    cutFile << muPoint->Parameter(0) << " "
	    << muPoint->Parameter(1) << " "
	    << muPoint->Parameter(2) << " "
	    << muPoint->Parameter(3) << " "
	    << sigmaPoint->Parameter(0) << " "
	    << sigmaPoint->Parameter(1) << " "
	    << sigmaPoint->Parameter(2) << " "
	    << sigmaPoint->Parameter(3) << "\n";

    //Now do sigma pdf
    sprintf(temp,"Canvas_Sigma_sec%d",i);
    TCanvas * c2 = new TCanvas(temp,temp,1200,1000);
    tg_sigma[i]->Draw();
    //f1_sigma_rga->Draw("SAME");
    f1_sigma->Draw("SAME");
    sprintf(temp,"Canvas_Sigma_sec%d.pdf",i);
    c2->SaveAs(temp);
  }
  


  inFile->Close();
  outFile->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  for(int i=0; i<graph_list.size(); i++){
    graph_list[i]->Write();
  }
  //outTree->Write();
  outFile->Close();
  cout<<"Finished making file: "<<argv[1]<<"\n";
  return 0;
}
