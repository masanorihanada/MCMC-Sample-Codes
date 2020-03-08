#include <iostream>
#include <cmath>
#include<fstream>
const int niter=1000;
const int ntau=10;
const double dtau=0.1e0;
/******************************************************************/
/*** Gaussian Random Number Generator with Box Muller Algorithm ***/
/******************************************************************/
int BoxMuller(double& p, double& q){
  
  double pi=2e0*asin(1e0);
  //uniform random numbers between 0 and 1
  double r = (double)rand()/RAND_MAX;
  double s = (double)rand()/RAND_MAX;
  //Gaussian random numbers, 
  //with weights proportional to e^{-p^2/2} and e^{-q^2/2}
  p=sqrt(-2e0*log(r))*sin(2e0*pi*s);
  q=sqrt(-2e0*log(r))*cos(2e0*pi*s);
  
  return 0;
}
/*********************************/
/*** Calculation of the action ***/
/*********************************/
// When you change the action, you should also change dH/dx,
// specified in "calc_delh". 
double calc_action(const double x){
  
  double action=0.5e0*x*x;
  
  return action;
}
/**************************************/
/*** Calculation of the Hamiltonian ***/
/**************************************/
double calc_hamiltonian(const double x,const double p){
  
  double ham=calc_action(x);
  
  ham=ham+0.5e0*p*p;
  
  return ham;
}
/****************************/
/*** Calculation of dH/dx ***/
/****************************/
// Derivative of the Hamiltonian with respect to x, 
// which is equivalent to the derivative of the action.
// When you change "calc_action", you have to change this part as well. 
double calc_delh(const double x){
  
  double delh=x;
  
  return delh;
}
/***************************/
/*** Molecular evolution ***/
/***************************/
int Molecular_Dynamics(double& x,double& ham_init,double& ham_fin){

  double r1,r2;
  BoxMuller(r1,r2);
  double p=r1;
  
  //*** calculate Hamiltonian ***
  ham_init=calc_hamiltonian(x,p);
  //*** first step of leap frog ***
  x=x+p*0.5e0*dtau;
  //*** 2nd, ..., Ntau-th steps ***
  for(int step=1; step!=ntau; step++){    
    double delh=calc_delh(x);
    p=p-delh*dtau;
    x=x+p*dtau;
  }
  //*** last step of leap frog ***
  double delh=calc_delh(x);
  p=p-delh*dtau;
  x=x+p*0.5e0*dtau;
  //*** calculate Hamiltonian again ***
  ham_fin=calc_hamiltonian(x,p);
  
  return 0;
}

int main()
{ 
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  double x=0e0;
  /*****************/
  /*** Main part ***/
  /*****************/
  std::ofstream outputfile("output.txt");
  int naccept=0;//counter for the number of acceptance
  double sum_xx=0e0;//sum of x^2, useed for <x^2>
  
  for(int iter=0; iter!=niter; iter++){
    
    double backup_x=x;
    double ham_init,ham_fin;
    Molecular_Dynamics(x,ham_init,ham_fin);
    double metropolis = (double)rand()/RAND_MAX;
    if(exp(ham_init-ham_fin) > metropolis){
      //accept
      naccept=naccept+1;
    }else{
      //reject
      x=backup_x;
    }
    
    /*******************/
    /*** data output ***/
    /*******************/
    sum_xx=sum_xx+x*x;

    // output x, <x^2>, acceptance
    
    std::cout << std::fixed << std::setprecision(6) 
      << x << "   " 
      << sum_xx/((double)(iter+1)) <<  "   " 
      << ((double)naccept)/((double)iter+1) 
      << std::endl;
     
    outputfile << std::fixed << std::setprecision(6) 
      << x << "   " 
      << sum_xx/((double)(iter+1)) <<  "   " 
      << ((double)naccept)/((double)iter+1) 
      << std::endl;
  }
    
  outputfile.close();
  
  return 0;
}
  


