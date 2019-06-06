#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TRandom.h"
#include "leggifile.cc"

#include <vector>
#include <iostream>

//#include "Style.C"

#define START_ORDER 0
#define NORDERS 8

#define OFFSET 755
struct FitResults {

  float ped_mu;
  float ped_mu_err;
  float ped_sigma;
  float ped_sigma_err;

  float mu;
  float mu_err;
  float offset;
  float offset_err;
  float Q1;
  float Q1_err;
  float sigma;
  float sigma_err;

};

Double_t PMTFunction(Double_t *x, Double_t *par)
{

  float N = par[0];
  float mu = par[1];
  float Q1 = par[2];
  float sigma = par[3];
  float offset = par[4];
  float sigmaoffset = par[5];
//   float alpha = par[6];
//   float w = par[7];
//   float frac = par[8];

  float xx = x[0];
  double value = 0.;

  for( unsigned i=START_ORDER; i<NORDERS; ++i ) {

    //double Qn = offset + (double)(i)*Q1;
    double sigma_n = sqrt( (double)(i)*sigma*sigma + sigmaoffset*sigmaoffset);

    //double poisson = TMath::Poisson( i, mu );
    //double gauss   = TMath::Gaus( xx, Qn, sigma_n );

    //double xxp = xx     - Qn - alpha*sigma_n*sigma_n;
    //double Q0p = offset - Qn - alpha*sigma_n*sigma_n;
    //double bg      = 0.5*alpha * TMath::Exp(-alpha*xxp)* (
    //                     TMath::Erf( abs(Q0p)/(sigma_n*sqrt(2) ) ) +  xxp/abs(xxp) * TMath::Erf( abs(xxp)/(sigma_n*sqrt(2)) ) );
    //value = value + N*( poisson * ( (1.-w)*gauss + w*bg ) );
    float norm;
//     if (i==0)
//       norm=N*(1+frac);
//     else
    norm=N;
    value += norm*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, ((double)i*Q1 + offset), sigma_n, kTRUE) );
    //value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sqrt((double)i)*sigma ));
  }

  return value;

}

FitResults fitSingleHisto( TH1F* histo, double xMin, double xMax ) 
{
  FitResults fr;
  TF1* f1 = new TF1( "fPMT", PMTFunction, xMin, xMax, 6 );
  f1->SetParameter( 0, histo->Integral()/5. ); //normalization
  f1->SetParameter( 1, 0.3); //poiss mu
  f1->SetParameter( 2, 22.74 ); //gauss step
  f1->SetParameter( 3, 12.13 ); //gauss sigma
  f1->SetParameter( 4, 29.16 ); //offset
  f1->SetParameter( 5, 4.21 ); //sigmaoffset
//   f1->SetParameter( 6, 0.03 ); //alpha
//   f1->SetParameter( 7, 0.4 ); //w
//   f1->SetParameter( 8, 0.3 ); //w

  f1->SetParName(0,"Norm");
  f1->SetParName(1,"#mu");
  f1->SetParName(2,"PE charge");
  f1->SetParName(3,"PE resolution");
  f1->SetParName(4,"Pedestal");
  f1->SetParName(5,"Noise");
  // f1->FixParameter( 1, 1. ); //mu
  //   f1->FixParameter( 2, 29.5 ); //Q1
  //   //  f1->FixParameter( 3, 14.3 ); //sigmaQ1
  //   f1->FixParameter( 4, 0. ); //offset
  //   f1->FixParameter( 5, 0. ); //sigmaoffset
  //    f1->FixParameter( 6, 0. ); //alpha
  //    f1->FixParameter( 7, 0. ); //w
  //    f1->FixParameter( 8, 0. ); //w
      
   f1->SetParLimits( 1, 0.1, 20  ); //poiss mu
   f1->SetParLimits( 2, 21., 27. ); //gauss step
   f1->SetParLimits( 3, 10., 14. ); //gauss sigma
   f1->SetParLimits( 4, 28, 33.); //offset
   f1->SetParLimits( 5, 3., 5. ); //gauss sigma


   
   // f1->FixParameter( 2, 22.74 ); //gauss step
   //f1->FixParameter( 3, 12.13 ); //gauss sigma
   //f1->FixParameter( 4, 29.16 ); //offset
   //f1->FixParameter( 5, 4.21 ); //sigmaoffset

   

    

  f1->SetLineColor(kRed+4);
  f1->SetLineWidth(6);
  
  
  histo->Fit( f1, "EMLR+" );
  TString histoName(histo->GetName());
  
  fr.mu = f1->GetParameter(1);
  fr.mu_err = f1->GetParError(1);
  delete f1;
  return fr;
}


void SinglePEAnalysis()
{


   gSystem->Exec("mkdir FitSingle");
  gSystem->Exec("mkdir Fileroot");
  gROOT->SetBatch("kTRUE");


  
 
  
  vector<string> FileList;
  string Filename="Filelist.txt";
  FileList=ReadData(Filename);
  bool printFileList= true;
  bool printOpenedFileNumber= true;
  if(printFileList){
    for(int i=0;i<(int)FileList.size();i++){
      std::cout<<"File " << i<< "   " << (FileList.at(i)).c_str() << std::endl;
      }
   }//chiudo if


  
  TH1F* adcData[FileList.size()];
  TFile *f[FileList.size()];
  TTree* tree[FileList.size()];
  TFile* out[FileList.size()];
  
  TCanvas *c[FileList.size()];
  //TFile *f=TFile::Open("h4Reco_calib1325V.root");
  TFile *v=TFile::Open("Voltaggi.root");
  TTree *volt=(TTree*)v->Get("volt");
  Int_t voltag;
  

  volt->SetBranchAddress("voltaggio", &voltag);
  TGraphErrors *calib= new TGraphErrors(FileList.size());
  TGraphErrors *mugr= new TGraphErrors(FileList.size());
  TGraphErrors *sigmape= new TGraphErrors(FileList.size());
  TGraphErrors *ped= new TGraphErrors(FileList.size());
  
  TGraphErrors *calib2= new TGraphErrors(FileList.size());
  TCanvas *canvCalib= new TCanvas("Calibration", "calibration", 1200,700);
  TCanvas *canvmugr= new TCanvas("Calibration", "calibration", 1200,700);
  TCanvas *canvsigmape= new TCanvas("Calibration", "calibration", 1200,700);
  TCanvas *canvped= new TCanvas("Calibration", "calibration", 1200,700);

  TF1* calibfit= new TF1("fitcalib","[0]+x*[1]",1220, 1370);

  
  for(int i=0; i<(int)FileList.size(); i++){
    volt->GetEntry(i);
    c[i]=new TCanvas(("c"+to_string(i)).c_str(),("c"+to_string(i)).c_str(),800,700);
    f[i]=TFile::Open((FileList.at(i).c_str()));
   if(printOpenedFileNumber) std::cout << "FileOpened: " << i << std::endl;
   out[i]=TFile::Open(("Fileroot/SinglePEAnalysis"+to_string(i)+".root").c_str(),"RECREATE");
   

   adcData[i]= new TH1F(("ledData"+to_string(i)).c_str(),("ledData"+to_string(i)).c_str(),600,0,400);
  tree[i]=(TTree*)f[i]->Get("digi");
  //TTree* tree=(TTree*)f->Get("digi");
  // gSystem->Exec("mkdir calib1350V");


 
  
  tree[i]->Project(("ledData"+to_string(i)).c_str(),"charge_tot[C0]");
  // adcData[i]->Print();
  //  std::cout << "FIT RANGE " << adcData[i]->GetMean()-3*adcData[i]->GetRMS() << "," << adcData[i]->GetMean()+3*adcData[i]->GetRMS() << std::endl;
   FitResults fr=fitSingleHisto(adcData[i],0,300);
  //  c->SetLogy(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(10111);
  adcData[i]->SetMarkerStyle(20);
  adcData[i]->SetMarkerSize(0.6);
  adcData[i]->SetMarkerColor(kBlack);
  adcData[i]->SetLineColor(kBlack);
  adcData[i]->Draw("PE");
  adcData[i]->GetXaxis()->SetTitle("Charge [ADC Counts]");
  //gPad->SetLogy();
  adcData[i]->GetYaxis()->SetTitle("N of events");

   for (int ipe=0; ipe<20;++ipe)
    {
      TF1* peFunc=new TF1(Form("peFunc_%d",ipe),"gausn",0,800);
      peFunc->SetLineColor(1+ipe);
      peFunc->SetLineWidth(2);
      float mu_pe=ipe*adcData[i]->GetFunction("fPMT")->GetParameter(2)+adcData[i]->GetFunction("fPMT")->GetParameter(4);
      float sigma_pe=sqrt(ipe*adcData[i]->GetFunction("fPMT")->GetParameter(3)*adcData[i]->GetFunction("fPMT")->GetParameter(3)+adcData[i]->GetFunction("fPMT")->GetParameter(5)*adcData[i]->GetFunction("fPMT")->GetParameter(5));
      peFunc->SetParameter(0,adcData[i]->GetFunction("fPMT")->GetParameter(0)*TMath::Poisson(ipe,adcData[i]->GetFunction("fPMT")->GetParameter(1)));
      peFunc->SetParameter(1,mu_pe);
      peFunc->SetParameter(2,sigma_pe);
      peFunc->Draw("SAME");
      }

  out[i]->cd();
  adcData[i]->Write();
  //c[i]->SaveAs(("FitSingle/singlePEfit"+to_string(i)+".png").c_str());
  c[i]->SaveAs(("FitSingle/singlePEfit"+to_string(i)+".pdf").c_str());
    out[i]->Write();
    //  if(voltag!=1350){
   
    mugr->SetPoint(i, voltag, adcData[i]->GetFunction("fPMT")->GetParameter(1));
    mugr->SetPointError(i, 0., adcData[i]->GetFunction("fPMT")->GetParError(1));

    calib->SetPoint(i, voltag, adcData[i]->GetFunction("fPMT")->GetParameter(2));
    calib->SetPointError(i, 0., adcData[i]->GetFunction("fPMT")->GetParError(2));


    ped->SetPoint(i, voltag, adcData[i]->GetFunction("fPMT")->GetParameter(4));
    ped->SetPointError(i, 0., adcData[i]->GetFunction("fPMT")->GetParError(4));
    //} else {
    //calib2->SetPoint(i, voltag, adcData[i]->GetFunction("fPMT")->GetParameter(2));
    //calib2->SetPointError(i, voltag, adcData[i]->GetFunction("fPMT")->GetParError(2));
    //}
    
  }
  

  canvCalib->cd();
   gPad->SetLogy(0);
  calib->SetMarkerStyle(20);
  calib->SetMarkerSize(0.6);
  calib->Draw("AP");
  canvCalib->SaveAs("singlepe.pdf");

  canvmugr->cd();
   gPad->SetLogy(0);
  mugr->SetMarkerStyle(20);
  mugr->SetMarkerSize(0.6);
  mugr->Draw("AP");
  canvmugr->SaveAs("mupmt.pdf");

  canvped->cd();
   gPad->SetLogy(0);
  ped->SetMarkerStyle(20);
  ped->SetMarkerSize(0.6);
  ped->Draw("AP");
  canvped->SaveAs("ped.pdf");

 

 
  
  
}
