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
const double step_mu=0.1e0;
const int nskip=10;
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
  
  action=action+0.5e0*A[0][0]*A[0][0]+0.5e0*A[0][0]*A[1][1]+A[0][1]*A[0][1]+0.5e0*mu[0]*mu[0]+0.5e0*mu[1]*mu[1];
  
  return action;
}
int main()
{
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  //don't choose a singular initial condition.
  double A[ndim][ndim];
  A[0][0]=1e0;A[1][1]=1e0;
  A[0][1]=0e0;
  A[1][0]=A[0][1];
  double mu[ndim];
  mu[0]=0e0;mu[1]=0e0;
  /*****************/
  /*** Main part ***/
  /*****************/
  int naccept=0;
  for(int iter=0; iter!=niter; iter++){
    double backup_A[ndim][ndim],backup_mu[ndim];
    backup_A[0][0]=A[0][0];
    backup_A[1][1]=A[1][1];
    backup_A[0][1]=A[0][1];
    backup_A[1][0]=A[1][0];
    backup_mu[0]=mu[0];
    backup_mu[1]=mu[1];

    double action_init=calc_action(A,mu);
    double dx;
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[0][0]=A[0][0]+dx*step_A*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[1][1]=A[1][1]+dx*step_A*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    A[1][0]=A[1][0]+dx*step_A*2e0;
    A[0][1]=A[0][1]+dx*step_A*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    mu[0]=mu[0]+dx*step_mu*2e0;
    dx = (double)rand()/RAND_MAX-0.5e0;
    mu[1]=mu[1]+dx*step_mu*2e0;
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
      mu[0]=backup_mu[0];
      mu[1]=backup_mu[1];
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
  


