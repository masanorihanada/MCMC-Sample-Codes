#include <iostream>
#include <cmath>
#include<fstream>
const int niter=10001000;
const double av_x=-0.0930181e0;
const double av_y=0.0475899e0;
const double av_xx=1.06614e0;
const double av_yy=1.28152e0;
const double av_xy=-0.504944e0;
const int nsample=100;
const int ndim=2;//don't chage this parameter.
const double step_A=0.1e0;
const int nskip=10;
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
double calc_action(const double A[ndim][ndim],const double mu[ndim]){
  
  double action=
    A[0][0]*((mu[0]-av_x)*(mu[0]-av_x)+av_xx-av_x*av_x)
    +A[1][1]*((mu[1]-av_y)*(mu[1]-av_y)+av_yy-av_y*av_y)
    +A[0][1]*((mu[0]-av_x)*(mu[1]-av_y)+av_xy-av_x*av_y)*2e0
    -log(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  action=0.5e0*action*nsample; 
  //add the contribution from the prior distribution
  action=action+0.5e0*A[0][0]*A[0][0]+0.5e0*A[0][0]*A[1][1]+A[0][1]*A[0][1]+0.5e0*mu[0]*mu[0]+0.5e0*mu[1]*mu[1];
  
  return action;
}
int main()
{
  double mu[ndim];
  double A[ndim][ndim];
  double dx;
  int naccept=0;
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  //NEVER choose a singular initial condition.
  A[0][0]=1e0;A[1][1]=1e0;
  A[0][1]=0e0;
  A[1][0]=A[0][1];
  mu[0]=0e0;mu[1]=0e0;
  /*****************/
  /*** Main part ***/
  /*****************/
  for(int iter=0; iter!=niter; iter++){
    /*************************/
    /* Gibbs sampling for mu */
    /*************************/
    double r1,r2;
    BoxMuller(r1,r2);
    mu[0]=r1/sqrt(1e0+nsample*A[0][0]);
    mu[0]=mu[0]-nsample/(1e0+nsample*A[0][0])*(A[0][1]*(mu[1]-av_y)-A[0][0]*av_x);

    mu[1]=r2/sqrt(1e0+nsample*A[1][1]);
    mu[1]=mu[1]-nsample/(1e0+nsample*A[1][1])*(A[1][0]*(mu[0]-av_x)-A[1][1]*av_y);
    /********************/
    /* Metropolis for A */
    /********************/
    double backup_A[ndim][ndim];
    backup_A[0][0]=A[0][0];
    backup_A[1][1]=A[1][1];
    backup_A[0][1]=A[0][1];
    backup_A[1][0]=A[1][0];

    double action_init=calc_action(A,mu);
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[0][0]=A[0][0]+dx*step_A*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[1][1]=A[1][1]+dx*step_A*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[1][0]=A[1][0]+dx*step_A*2e0;
    A[0][1]=A[0][1]+dx*step_A*2e0;
    double action_fin=calc_action(A,mu);
    /*******************/
    /* Metropolis test */
    /*******************/
    double metropolis = (double)rand()/RAND_MAX;    
    if(exp(action_init-action_fin) > metropolis){
      /* accept */
      naccept=naccept+1;
    }else{ 
      /* reject */
      A[0][0]=backup_A[0][0];
      A[1][1]=backup_A[1][1];
      A[0][1]=backup_A[0][1];
      A[1][0]=backup_A[1][0];
    }

    // output 
    if((iter+1)%nskip==0){
      std::cout << std::fixed << std::setprecision(6) 
        << A[0][0] << "  "
        << A[1][1] << "  "
        << A[0][1] << "  "
        << mu[0] << "  "
        << mu[1] << "  "
        << (double) naccept/(iter+1) << std::endl;
    }
  } 
  
  return 0;
}
  


