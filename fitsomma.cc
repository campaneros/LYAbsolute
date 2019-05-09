#define START_ORDER 0
#define NORDERS 6

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


Double_t PMTFunction(Double_t *x, Double_t*par)
{

  float N = par[3];
  float mu = par[0];
  float Q1 = par[4];
  float sigma = par[5];
  float offset = par[1];
  float sigmaoffset = par[2];
//   float alpha = par[6];
//   float w = par[7];
//   float frac = par[8];

  float xx = x[0];
  double value = 0.;

  for( unsigned i=START_ORDER; i<NORDERS; ++i ) {

   
    double sigma_n = sqrt( (double)(i)*sigma*sigma + sigmaoffset*sigmaoffset);
    float norm;
     norm=N;
     value += norm*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, ((double)i*Q1 + offset), sigma_n, kTRUE) );
     }

  return value;

}


/*FitResults fitSingleHisto( TH1F* histo, double xMin, double xMax ) 
{
  FitResults fr;
  TF1* f1 = new TF1( "fPMT", PMTFunction, xMin, xMax, 6 );
  TF1 *gausPed= new TF1("gausPed","gaus", 0, 35);
  TF1 *gausSe= new TF1("gausSE","gaus", 45,60);
  histo->Fit("gausPed", "RN");
  histo->Fit("gausSE", "RN");
  

  f1->SetParameter( 0, histo->Integral()/5. ); //normalization
  f1->SetParameter( 1, 0.3); //poiss mu
  f1->SetParameter( 2, gausPed->GetParameter(1)); //gauss step
  f1->SetParameter( 3, (gausSe->GetParameter(2)-gausPed->GetParameter(2)) ); //gauss sigma
  f1->SetParameter( 4, (gausSe->GetParameter(1)-gausPed->GetParameter(1)) ); //offset
  f1->SetParameter( 5, gausPed->GetParameter(2) ); //sigmaoffset


  f1->SetParName(0,"Norm");
  f1->SetParName(1,"#mu");
  f1->SetParName(2,"PE charge");
  f1->SetParName(3,"PE resolution");
  f1->SetParName(4,"Pedestal");
  f1->SetParName(5,"Noise");

      
   f1->SetParLimits( 1, 0.1, 5.  ); //poiss mu
   f1->SetParLimits( 2, 20., 30. ); //gauss step
   f1->SetParLimits( 3, 10., 15. ); //gauss sigma
   f1->SetParLimits( 4, 25., 35.); //offset
   f1->SetParLimits( 5, 2., 10. ); //gauss sigma

  f1->SetLineColor(kCyan+2);
  f1->SetLineWidth(2);
  
  
  histo->Fit( f1, "LR+" );
  histo->Fit( f1, "LR+" );
  
  //f1->SetLineColor(kRed);
  TString histoName(histo->GetName());
  
  fr.mu = f1->GetParameter(1);
  fr.mu_err = f1->GetParError(1);
  delete f1, gausPed, gausSe;
  return fr;
  
  }*/
