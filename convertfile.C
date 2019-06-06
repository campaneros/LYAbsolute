#include "TString.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TStyle.h"



void Fileconvert(){

  TTree *volt=new TTree("volt","volt");
  volt->ReadFile("Voltaggi.dat","voltaggio/I");
  TFile *out=TFile::Open("Voltaggi.root","RECREATE");
  //volt->Branch("voltaggio", 
  volt->Print();
  out->cd();
  volt->Write();
  out->Close();
  TFile *Read=TFile::Open("Voltaggi.root");
  TTree *t=(TTree*)Read->Get("volt");
  Int_t prova;
  t->SetBranchAddress("voltaggio",&prova);
  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    std::cout<<prova<<std::endl;
  }

  
}
