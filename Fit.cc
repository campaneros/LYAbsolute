


Double_t Pedestal ( Double_t *x, Double_t *par){

  Double_t poisson, gauss, sigma, fit;
 

  //poisson = exp(-par[0]);
  poisson=TMath::Poisson(0,par[0]);
  sigma=par[2];
  gauss=par[3]*exp(-pow((x[0]-par[1]),2)/(2*sigma*sigma));//(sigma*sqrt(2*M_PI));
   //gauss=par[3]*TMath::Gaus(x[0], par[1], sigma);
  fit=poisson*gauss;
  

  return fit;

  
   }




Double_t singlePe ( Double_t *x, Double_t *par){

  Double_t poisson, gauss, sigma, fit;


  //poisson = exp(-par[0])*par[0];
  poisson=TMath::Poisson(1,par[0]);
  sigma=sqrt(par[2]*par[2]+par[5]*par[5]);
  gauss=par[3]*exp(-pow((x[0]-(par[1]+par[4])),2)/(2*sigma*sigma));//(sigma*sqrt(2*M_PI));
  // gauss=par[3]*TMath::Gaus(x[0], (par[1]+par[4]), sigma);
  fit=poisson*gauss;
  return fit;
  
  //return exp(-par[0])*par[0]*exp(-((x[0]-(par[1]+par[3]))*(x[0]-(par[1]+par[3])))/(2*(sqrt(par[2]*par[2]+par[3]*par[3])*(sqrt(par[2]*par[2]+par[3]*par[3])))))/((sqrt(par[2]^2+par[3]^2))*sqrt(2*pi));
   }

  Double_t DoublePe ( Double_t *x, Double_t *par){

    
 Double_t poisson, gauss, sigma, fit;


 //poisson = exp(-par[0])*(par[0]*par[0])/2;
 poisson=TMath::Poisson(2,par[0]);
  sigma=sqrt(par[2]*par[2]+2*par[5]*par[5]);
  gauss=par[3]*exp(-pow((x[0]-(par[1]+2*par[4])),2)/(2*sigma*sigma));//(sigma*sqrt(2*2*M_PI));
  // gauss=par[3]*TMath::Gaus(x[0], (par[1]+2*par[4]), sigma);
  fit=poisson*gauss;
  return fit;



    //return exp(-par[0])*par[0]^2*exp(-((x[0]-(par[1]+2*par[3]))^2)/(2*(sqrt(par[2]^2+par[4]^2)))^2)/((sqrt(par[2]^2+par[4]^2))*sqrt(2*2*pi));
   }


Double_t TriplePe ( Double_t *x, Double_t *par){

    
 Double_t poisson, gauss, sigma, fit;


 poisson = exp(-par[0])*pow(par[0],3)/6;
  sigma=sqrt(par[2]*par[2]+3*par[5]*par[5]);
  gauss=par[3]*exp(-pow((x[0]-(par[1]+3*par[4])),2)/(2*sigma*sigma))/(sigma*sqrt(2*3*M_PI));
  fit=poisson*gauss;
  return fit;



    //return exp(-par[0])*par[0]^2*exp(-((x[0]-(par[1]+2*par[3]))^2)/(2*(sqrt(par[2]^2+par[4]^2)))^2)/((sqrt(par[2]^2+par[4]^2))*sqrt(2*2*pi));
   }




  Double_t Fit ( Double_t *x, Double_t *par){
    return Pedestal(x,par)+singlePe(x,par)+DoublePe(x,par);//TriplePe(x,par);
    
    }

 Double_t Fitn (Double_t *x, Double_t*par){
  return singlePe(x,par)+DoublePe(x,par);//+TriplePe(x,par);
 }



Double_t singlePen ( Double_t *x, Double_t *par){

 Double_t poisson, gauss, sigma, fit;


  poisson = exp(-par[0])*par[0];
  gauss=par[3]*exp(-((pow((x[0]-par[1]),2))/(2*pow(par[2],2))))/(par[2]*sqrt(2*M_PI));
  fit=poisson*gauss;
  return fit;
  }


Double_t DoublePen ( Double_t *x, Double_t *par){

 Double_t poisson, gauss, sigma, fit;


  poisson = exp(-par[0])*par[0]*par[0]/2;
  gauss=par[3]*exp(-(pow((x[0]-par[4]),2))/(2*pow(par[5],2)))/(par[5]*sqrt(2*2*M_PI));
  fit=poisson*gauss;
  return fit;
  }

//Double_t Fitn (Double_t *x, Double_t*par){
//return singlePen(x,par)+DoublePen(x,par);
//}




