


Double_t Pedestal ( Double_t *x, Double_t *par){

  Double_t poisson, gauss, sigma, fit;


  poisson = exp(-par[0]);
  sigma=par[2];
  gauss=par[6]*exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
  fit=poisson*gauss;
  

  return fit;

  
   }



Double_t singlePe ( Double_t *x, Double_t *par){

 Double_t poisson, gauss, sigma, fit;


  poisson = exp(-par[0])*par[0];
  sigma=sqrt(par[2]*par[2]+par[4]*par[4]);
  gauss=par[6]*exp(-((x[0]-(par[1]+par[3]))*(x[0]-(par[1]+par[3])))/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
  fit=poisson*gauss;
  return fit;
  
  //return exp(-par[0])*par[0]*exp(-((x[0]-(par[1]+par[3]))*(x[0]-(par[1]+par[3])))/(2*(sqrt(par[2]*par[2]+par[3]*par[3])*(sqrt(par[2]*par[2]+par[3]*par[3])))))/((sqrt(par[2]^2+par[3]^2))*sqrt(2*pi));
   }

  Double_t DoublePe ( Double_t *x, Double_t *par){

    
 Double_t poisson, gauss, sigma, fit;


  poisson = exp(-par[0])*par[0]*par[0];
  sigma=sqrt(par[2]*par[2]+par[5]*par[5]);
  gauss=par[6]*exp(-((x[0]-(par[1]+2*par[3]))*(x[0]-(par[1]+2*par[3])))/(2*sigma*sigma))/(sigma*sqrt(2*2*M_PI));
  fit=poisson*gauss;
  return fit;



    //return exp(-par[0])*par[0]^2*exp(-((x[0]-(par[1]+2*par[3]))^2)/(2*(sqrt(par[2]^2+par[4]^2)))^2)/((sqrt(par[2]^2+par[4]^2))*sqrt(2*2*pi));
   }


  Double_t Fit ( Double_t *x, Double_t *par){
    return Pedestal(x,par)+singlePe(x,par)+DoublePe(x,par);
    
    }
    




