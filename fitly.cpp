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
#include "TLatex.h"
#include "leggifile.cc"
#include "Fit.cc"
#include "fitsomma.cc"



#include "TGraphErrors.h"
#include "TMath.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TRandom.h"









int main(){
   gROOT->Reset();
   gROOT->SetBatch(kFALSE);

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
   Double_t meanPed, meanSe, sigmaPed, sigmaSe, meanDe, sigmaDe, FitSomma;
   


   
   TH1F* HistoCharge=new TH1F ("histocharge","histocharge", 600, 0, 200);
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
	    DataTree= (TTree*)DataFile->Get("h4");
	    DataTree->SetBranchAddress("charge_tot",&Charge); 

	    
	    
	    for(int i=0;i<DataTree->GetEntries();i++){
	      // std::cout<<i<<std::endl;
	      DataTree->GetEntry(i);
	      HistoCharge->Fill(Charge);
	      }



	   
	    
	    TF1 *fitcharge= new TF1("fitcharge",Fit, 15, 90, 6);
	    TF1 *fitchargee= new TF1("fitchargee", Fitn, 40, 160, 6);
	    TF1 *termzero= new TF1("termzero", Pedestal, 0, 60,4);
	    TF1 *termone= new TF1("termone", singlePe, 0, 90, 6);
	    TF1 *termtwo= new TF1("termtwo", DoublePe, 15, 120, 6);
	    TF1 *termthree= new TF1("termthree", TriplePe, 15, 160, 6);
	    TF1 *gausPed= new TF1("gausPed","gaus", 0, 35);
	    TF1 *gausSe= new TF1("gausSE","gaus", 45,60);
	    TF1 *gausDe= new TF1("gausDe","gaus", 80, 160);
	    TF1 *fitsomma= new TF1("fitsomma",PMTFunction, 15, 90, 6);
	    
	     

	     HistoCharge->Fit("gausPed", "RN");
	    

	     meanPed=gausPed->GetParameter(1);
	     sigmaPed=gausPed->GetParameter(2);
	     
	     // gausSe->SetParLimits(1,0,60);
	     //gausSe->SetParameter(1,42);
	     HistoCharge->Fit("gausSE", "RN");

	       meanSe=gausSe->GetParameter(1);
	       sigmaSe=gausSe->GetParameter(2);

	      HistoCharge->Fit("gausDe", "RN");

	       meanDe=gausDe->GetParameter(1);
	       sigmaDe=gausDe->GetParameter(2);
	     
	     
	     
	    fitcharge->SetParameter(0, 0.3);
	    fitcharge->SetParameter(1,meanPed);
	    fitcharge->SetParameter(2,sigmaPed);
	    fitcharge->SetParameter( 3, HistoCharge->Integral()/5. );
	    fitcharge->SetParameter(4,(meanSe-meanPed));
	    fitcharge->SetParameter(5,(sigmaSe-sigmaPed));

	     fitcharge->SetParLimits( 0, 0.1, 5.  ); //poiss mu
	     fitcharge->SetParLimits( 4, 10., 40. ); //gauss step
	     fitcharge->SetParLimits( 5, 8., 30. ); //gauss sigma
	     fitcharge->SetParLimits( 1, 10, 30.); //offset
	     fitcharge->SetParLimits( 2, 2., 10. ); //gauss sigma

	     fitsomma->SetParameter(0, 0.3);
	    fitsomma->SetParameter(1,meanPed);
	    fitsomma->SetParameter(2,sigmaPed);
	    fitsomma->SetParameter( 3, HistoCharge->Integral()/5. );
	    fitsomma->SetParameter(4,(meanSe-meanPed));
	    fitsomma->SetParameter(5,(sigmaSe-sigmaPed));

	     fitsomma->SetParLimits( 0, 0.1, 5.  ); //poiss mu
	     fitsomma->SetParLimits( 4, 10., 40. ); //gauss step
	     fitsomma->SetParLimits( 5, 8., 30. ); //gauss sigma
	     fitsomma->SetParLimits( 1, 10, 30.); //offset
	     fitsomma->SetParLimits( 2, 2., 10. ); //gauss sigma




	       

	   
	    /* fitcharge->SetParLimits(6,0,100);
	    fitcharge->SetParameter(6,(sigmaDe-sigmaPed));
	    fitcharge->SetParLimits(7,0,100);
	    fitcharge->SetParameter(7,(sigmaSe-sigmaPed));*/
	    fitcharge->SetParNames("#mu","mean_{ped}","#sigma_{ped}","constant","SinglePhotoElectron","#sigma_{1}");
	 

	    
	    
	     fitsomma->SetParNames("#mu","mean_{ped}","#sigma_{ped}","constant","SinglePhotoElectron","#sigma_{1}");
	    //   fitchargee->SetParLimits(6,0,100);
	    //fitchargee->SetParameter(6,(sigmaDe-sigmaPed));
	    //fitchargee->SetParameter(7,(sigmaSe-sigmaPed));
	    // fitchargee->SetParameter(6,24);
	    
	    /* fitchargee->SetParLimits(0,0,1);
	    fitchargee->SetParameter(1,meanSe);
	    fitchargee->SetParameter(2,sigmaSe);
	    fitchargee->SetParLimits(1,45,70);  
	    fitchargee->SetParLimits(4,45,80);
	    fitchargee->SetParLimits(5,0,100);;
	    fitchargee->SetParameter(5,sigmaSe*sqrt(2));*/


	    

	     gStyle->SetOptStat(1111);
	     gStyle->SetOptFit(0011);

	    canvCharge->cd();
	    fitsomma->SetLineWidth(3);
	    fitcharge->SetLineWidth(3);
	    //cavCharge->GetXaxis()->SetLimits(0,37);
	     HistoCharge->Fit("fitcharge", "RLN");
	     HistoCharge->Fit("fitsomma", "R");
	     // fitchargee->SetParLimits(0,0.7,2);
	     // fitchargee->FixParameter(1,meanPed);
	     // fisomm->FixParameter(2,sigmaPed);
	    //fitchargee->SetParLimits(3,3100,5000);
	    fitchargee->SetParLimits(4,20,30);  
	    fitchargee->SetParameter(4,(meanSe-meanPed));
	    fitchargee->SetParameter(5,(sigmaSe-sigmaPed));
	    HistoCharge->Draw("PE");
	    fitcharge->SetLineColor(kGreen);
	    fitcharge->Draw("SAME");
	    //  HistoCharge->Fit("fitcharge", "RN");
	    termzero->FixParameter(0,fitcharge->GetParameter(0));
	    termzero->FixParameter(1,fitcharge->GetParameter(1));
	    termzero->FixParameter(2,fitcharge->GetParameter(2));
	    termzero->FixParameter(3,fitcharge->GetParameter(3));
	    termzero->SetLineColor(kBlue);
	    termzero->SetLineStyle(2);
	    //  termzero->Draw("SAME");
	    termone->FixParameter(0,fitcharge->GetParameter(0));
	    termone->FixParameter(1,fitcharge->GetParameter(1));
	    termone->FixParameter(2,fitcharge->GetParameter(2));
	    termone->FixParameter(3,fitcharge->GetParameter(3));
	    termone->FixParameter(4,fitcharge->GetParameter(4));
	    termone->FixParameter(5,fitcharge->GetParameter(5));
	    termone->SetLineColor(kBlue);
	    termone->SetLineStyle(2);
	    termone->Draw("SAME");
	    termtwo->FixParameter(0,fitcharge->GetParameter(0));
	    termtwo->FixParameter(1,fitcharge->GetParameter(1));
	    termtwo->FixParameter(2,fitcharge->GetParameter(2));
	    termtwo->FixParameter(3,fitcharge->GetParameter(3));
	    termtwo->FixParameter(4,fitcharge->GetParameter(4));
	    termtwo->FixParameter(5,fitcharge->GetParameter(5));
	    termtwo->SetLineColor(kBlue);
	    termtwo->SetLineStyle(2);
	    termtwo->Draw("SAME");
	    termthree->FixParameter(0,fitcharge->GetParameter(0));
	    termthree->FixParameter(1,fitcharge->GetParameter(1));
	    termthree->FixParameter(2,fitcharge->GetParameter(2));
	    termthree->FixParameter(3,fitcharge->GetParameter(3));
	    termthree->FixParameter(4,fitcharge->GetParameter(4));
	    termthree->FixParameter(5,fitcharge->GetParameter(5));
	    termthree->SetLineColor(kGreen);
	    termthree->SetLineStyle(2);
	    TText *t=new TText();
	    t->SetTextFont(42);
	    // t->DrawText(40,1400,"Pedestal");
	    
	    // termthree->Draw("SAME");
	    //fitcharge->Draw("SAME");
	    //fitchargee->Draw("SAME");
	    

	     canvCharge->SaveAs(("Plot2/Chargedistri"+Filename+".png").c_str());
	    
	    return 0;

	    
 }
