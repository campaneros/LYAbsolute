#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <math.h>



#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"



#include "leggifile.cc"
#include "Fit.cc"




/*Double_t Fit ( Double_t *x, Double_t *par){

  Double_t poisson, gauss, sigma, ped, single, doub;


  poisson = exp(-par[0]);
  sigma=par[2];
  gauss=exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*sigma*sigma))/(sigma*sqrt(2*3.14));
  ped=poisson*gauss;


  poisson = exp(-par[0])*par[0];
  sigma=sqrt(par[2]*par[2]+par[4]*par[4]);
  gauss=exp(-((x[0]-(par[1]+par[3]))*(x[0]-(par[1]+par[3])))/(2*sigma*sigma))/(sigma*sqrt(2*3.14));
  single=poisson*gauss;


  poisson = exp(-par[0])*par[0]*par[0];
  sigma=sqrt(par[2]*par[2]+par[5]*par[5]);
  gauss=exp(-((x[0]-(par[1]+2*par[3]))*(x[0]-(par[1]+2*par[3])))/(2*sigma*sigma))/(sigma*sqrt(2*2*3.14));
  doub=poisson*gauss;
  

  return ped+single+doub;

  
  }*/

int main(){
   gROOT->Reset();
   gROOT->SetBatch(kTRUE);

   std::string File, DirPlot;


  TFile* DataFile;
  TTree* DataTree;
   
  //std::cout<<"Inserisci nome File: "<<std::endl;
  //getline(cin, File);
   //std::cout<<"Inserisci nome Directory: "<<std::endl;
    // getline(cin, DirPlot);
   // gSystem->Exec(("mkdir "+DirPlot+"/Plot").c_str());
   gSystem->Exec("mkdir Plot2");
   
   vector<string> FileList;
   string Filename="Filelist2.txt";
   FileList=ReadData(Filename);
   Float_t Charge;
   Double_t meanPed, meanSe, sigmaPed, sigmaSe;
   


   
   TH1F* HistoCharge=new TH1F ("histocharge","histocharge", 600, 0, 180);
   TCanvas* canvCharge=new TCanvas ("charge","charge",1200,600);
     

   /*   string Filename= DirData+"/Filelist.txt";
    std::cout << "Lista File Pedestal: "<< Filename<< std::endl;
    FileList=ReadData(Filename);*/

    bool printFileList= true;
  bool printOpenedFileNumber= true;



  if(printFileList){
    
    // for(int i=0;i<(int)FileList.size();i++){
      std::cout<<"File " << 0<< "   " << (FileList.at(0)).c_str() << std::endl;
      //}//chiudo for
    
   }//chiudo if

  //size=(int)FileList.size();



   DataFile=TFile::Open((FileList.at(0)).c_str());
            if(printOpenedFileNumber) std::cout << "FileOpened: " << 0 << std::endl;
	    DataTree= (TTree*)DataFile->Get("digi");

	    DataTree->SetBranchAddress("charge_tot",&Charge); 

	    
	    
	    for(int i=0;i<DataTree->GetEntries();i++){
	      DataTree->GetEntry(i);
	      HistoCharge->Fill(Charge);
	      }


	    
	    TF1 *fitcharge= new TF1("fitcharge",Fit, 0, 160, 7);

	     TF1 *gausPed= new TF1("gausPed","gaus", 0, 35);
	     TF1 *gausSe= new TF1("gausSE","gaus", 40, 80);
	     //  TF1 *gausDe= new TF1("gausSE","gaus", 80, 120);


	     HistoCharge->Fit("gausPed", "RN");
	    

	     meanPed=gausPed->GetParameter(1);
	     sigmaPed=gausPed->GetParameter(2);

	      HistoCharge->Fit("gausSE", "RN");

	       meanSe=gausSe->GetParameter(1);
	       sigmaSe=gausSe->GetParameter(2);


	     

	    fitcharge->SetParLimits(0,0,1);
	    fitcharge->SetParameter(1,meanPed);
	    fitcharge->SetParameter(2,sigmaPed);
	    fitcharge->SetParameter(3,(meanSe-meanPed));
	    fitcharge->SetParameter(4,(sigmaSe-sigmaPed));
	    fitcharge->SetParameter(5,(sigmaSe-sigmaPed));

	    
	    //  fitcharge->SetParLimits(0,-6,-4);    //questo -4.16
	    // fitcharge->SetParLimits(1,60,120);  //questo 47
	    // fitcharge->SetParLimits(2,4,40); //questo 19.29
	    // fitcharge->SetParLimits(3,0,20);
	    //fitcharge->SetParLimits(4,0,20);
	    // fitcharge->SetParLimits(5,0,20);

	     gStyle->SetOptStat(1111);
	     gStyle->SetOptFit(0011);

	    canvCharge->cd();
	    //cavCharge->GetXaxis()->SetLimits(0,37);
	    HistoCharge->Fit("fitcharge", "R");
	    
	    HistoCharge->Draw();

	    canvCharge->SaveAs(("Plot2/Chargedistri"+Filename+".png").c_str());
	    

	    
 }
