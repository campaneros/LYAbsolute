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


struct FitResults {
  float Ngaus;
  float meanpeak;
  float sigmagaus;
  float Nfrac;
  float boh;
};






Double_t Spectrum (Double_t *x, Double_t *par)
{
  float xx=x[0];
  float Ngaus=par[0];
  float meanpeak= par[1];
  float sigmagaus=par[2];
  float Nfrac=par[3];
  float expo=par[4];
  float Ngaus2=par[5];
  float meanpeak2=par[6];
  float sigmagaus2=par[7];
  float mu=par[8];


  double value;

  value = Ngaus*TMath::Gaus(xx,meanpeak,sigmagaus);//+Nfrac/(TMath::Exp(xx*expo-(2*meanpeak*meanpeak/(meanpeak+2*meanpeak)))+1);//+Ngaus2*TMath::Gaus(xx, meanpeak2, sigmagaus2);//TMath::crystalball_function(par[8],par[9],par[10],par[11],par[12]);//+ROOT::Math::crystalball_function(xx, alpha, n, sigma, mu);
  return value;
}

FitResults fitSingleHisto( TH1F* histo, double xMax)
{
  FitResults fr;
  Double_t peak, peak2;
  Double_t *peaak;
  Double_t RMS, RMS2;
  Double_t min, min2;
  Double_t max, max2;
  Int_t npeaks, EndPlot;


 
  TSpectrum *s = new TSpectrum();
  /* peak = histo->GetBinCenter(histo->GetMaximumBin());
  std::cout<<"il picco è a ="<<peak<<std::endl;
  peaak =s->GetPositionX();
  std::cout<<"il picco è a ="<<*peaak<<std::endl;*/
  // f1->SetParName(5,"Noise");


  // Int_t nfound = s->Search(histo,4,"nobackground""goff",0.50);
  //peaak =s->GetPositionX();
  histo->GetXaxis()->SetRangeUser(6000,13000);
  peak=histo->GetBinCenter(histo->GetMaximumBin());
  histo->GetXaxis()->UnZoom();
  /* while( peak<20000){
    int i=1;
    i+=1;
    peak=peaak[nfound-i];
    }*/
    std::cout<<"il picco è a ="<<peak<<std::endl;
  // std::cout<<"il picco è a ="<<nfound<<std::endl;
   std::cout<<"il picco è a ="<<((histo->GetMaximum()))<<std::endl;



  
  histo->GetXaxis()->SetRangeUser(peak/2,peak);
  min = histo->GetBinCenter(histo->GetMinimumBin());
  histo->GetXaxis()->UnZoom();

   histo->GetXaxis()->SetRangeUser( peak - (peak - min) , peak+ (peak - min) );
  TF1* gaus=new TF1("gaus", "gaus",(peak -RMS), (peak+ RMS));
  histo->Fit(gaus,"N");
  max=histo->GetMaximum();
  RMS = gaus->GetParameter(2);
  histo->GetXaxis()->UnZoom();
 

  RMS=gaus->GetParameter(2);
  peak=gaus->GetParameter(1);


   Int_t Last= histo->GetXaxis()->GetLast();
  std::cout << Last << std::endl;

  for(int i=Last;i>0;i--){
    if( histo->GetBinContent(i) > 100) {
      EndPlot=i;
      break;
    }
  }
  
  std::cout << EndPlot << std::endl;

  peak2=EndPlot/1.09;

  std::cout << peak2 << std::endl;
  
  histo->GetXaxis()->SetRangeUser( peak2 - ( EndPlot - peak2 ) , EndPlot );
  peak2 = histo->GetBinCenter(histo->GetMaximumBin());
  RMS2 = histo->GetRMS();
  histo->GetXaxis()->UnZoom();
  //peak2=histo->GetBinCenter(histo->GetMaximumBin());
  /* std::cout<<"il picco2 è a ="<<peak2<<std::endl;
  histo->GetXaxis()->UnZoom();
  histo->GetXaxis()->SetRangeUser(peak2/2,peak2);
  min2 = histo->GetBinCenter(histo->GetMinimumBin());
  histo->GetXaxis()->UnZoom();
  histo->GetXaxis()->SetRangeUser( peak2 - (peak2 - min2) , peak2+ (peak2 - min2) );
  //TF1* gaus1=new TF1("gaus1", "gaus",(peak -RMS), (peak+ RMS));
  //histo->Fit(gaus1,"N");
  //max2=histo->GetMaximum();
  //RMS2 = gaus1->GetParameter(2);
  histo->GetXaxis()->UnZoom();*/
  


  TF1* f1 = new TF1( "fSpectrum", Spectrum, peak-RMS, peak+RMS, 3);

  f1->SetParName(0,"Norm");
  f1->SetParName(1,"meangaus");
  f1->SetParName(2,"sigmagaus");
  f1->SetParName(3,"Nfrac");
  f1->SetParName(4,"expo");
  f1->SetParName(5,"Norm2");
  f1->SetParName(6,"meangaus2");
  f1->SetParName(7,"sigmagaus2");

  f1->SetParameter(0,max);
  f1->SetParameter(1,peak);
  f1->SetParameter(2,RMS);
  f1->SetParameter(3, 40);
  f1->SetParameter(5,max2);
  f1->SetParameter(6,peak2);
  f1->SetParameter(7,RMS2);
  
  
  
  


  histo->Fit( f1, "LR+" );
  TString histoName(histo->GetName());
  fr.meanpeak= f1->GetParameter(1);
  fr.Nfrac=f1->GetParameter(3);
  fr.Ngaus=f1->GetParameter(0);
  fr.sigmagaus=f1->GetParameter(2);
  fr.boh=f1->GetParameter(4);

  
  delete f1;
  return fr;

}


void AnalysisBarWrap()
{


   gSystem->Exec("mkdir FitSingle");
  gSystem->Exec("mkdir Fileroot");
  //gROOT->SetBatch("kTRUE");


  
 
  
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


  /*Float_t lowbin[5];

  lowbin[0]=7000;
  lowbin[1]=7300;
  lowbin[2]=7600;
  lowbin[3]=7900;
  lowbin[4]=8200;*/
  
  
 

  TH1F* HistoMean= new TH1F("Mean","Mean",10,7000,8200);

  
  for(int i=0; i<(int)FileList.size(); i++){
    c[i]=new TCanvas(("c"+to_string(i)).c_str(),("c"+to_string(i)).c_str(),800,700);
    f[i]=TFile::Open((FileList.at(i).c_str()));
   if(printOpenedFileNumber) std::cout << "FileOpened: " << i << std::endl;
   out[i]=TFile::Open(("Fileroot/SinglePEAnalysis"+to_string(i)+".root").c_str(),"RECREATE");
   

   adcData[i]= new TH1F(FileList.at(i).c_str(),FileList.at(i).c_str(),300,0,30000);
  tree[i]=(TTree*)f[i]->Get("digi");
  //TTree* tree=(TTree*)f->Get("digi");
  // gSystem->Exec("mkdir calib1350V");


 
  
  tree[i]->Project(FileList.at(i).c_str(),"charge_tot[C0]");
  // adcData[i]->Print();
  //  std::cout << "FIT RANGE " << adcData[i]->GetMean()-3*adcData[i]->GetRMS() << "," << adcData[i]->GetMean()+3*adcData[i]->GetRMS() << std::endl;
   FitResults fr=fitSingleHisto(adcData[i],30000);
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


  out[i]->cd();
  adcData[i]->Write();
  //c[i]->SaveAs(("FitSingle/singlePEfit"+to_string(i)+".png").c_str());
  c[i]->SaveAs(("FitSingle/singlePEfit"+to_string(i+1)+".pdf").c_str());
    out[i]->Write();
    //  if(voltag!=1350){


     HistoMean->Fill(fr.meanpeak);
    
  }
  


  gStyle->SetOptStat(1);
  HistoMean->GetXaxis()->SetRangeUser(6700,8500);
  HistoMean->Draw();
  
 
  
  
}
