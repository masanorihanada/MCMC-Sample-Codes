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
const int ntau=50;
const double dtau_A=0.05e0;
const double dtau_mu=0.05e0;
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
  action=0.5e0*action*(double)nsample;
  
  action=action+0.5e0*A[0][0]*A[0][0]+0.5e0*A[1][1]*A[1][1]+A[0][1]*A[0][1]+0.5e0*mu[0]*mu[0]+0.5e0*mu[1]*mu[1];
  
  return action;
}
/**************************************/
/*** Calculation of the Hamiltonian ***/
/**************************************/
double calc_hamiltonian(const double p_A[ndim][ndim], const double p_mu[ndim],
                        const double A[ndim][ndim], const double mu[ndim]){
  
  double ham=calc_action(A,mu);
   
  for(int idim=0; idim!=ndim; idim++){
    ham=ham+0.5e0*p_mu[idim]*p_mu[idim];
    for(int jdim=0; jdim!=ndim; jdim++){
      ham=ham+0.5e0*p_A[idim][jdim]*p_A[jdim][idim];
    }
  }

  return ham;
}
/****************************/
/*** Calculation of dH/Dx ***/
/****************************/
// Derivatives of the Hamiltonian with respect to A and mu, 
// which is equivalent to the derivative of the action.
// When you change "calc_action", you have to change this part as well. 
int calc_delh(const double A[ndim][ndim],const double mu[ndim],
              double (&delh_A)[ndim][ndim], double (&delh_mu)[ndim]){

  double A_inv[ndim][ndim];
  double det_A=A[0][0]*A[1][1]-A[0][1]*A[1][0];
  A_inv[0][0]=A[1][1]/det_A;
  A_inv[1][1]=A[0][0]/det_A;
  A_inv[0][1]=-A[0][1]/det_A;
  A_inv[1][0]=-A[1][0]/det_A;

  delh_A[0][0]=A[0][0]+(mu[0]*mu[0]-2e0*av_x*mu[0]+av_xx-A_inv[0][0])*0.5e0*nsample;
  delh_A[1][1]=A[1][1]+(mu[1]*mu[1]-2e0*av_y*mu[1]+av_yy-A_inv[1][1])*0.5e0*nsample;
  delh_A[0][1]=A[0][1]+(mu[0]*mu[1]-av_x*mu[1]-av_y*mu[0]+av_xy-A_inv[0][1])*0.5e0*nsample;
  delh_A[1][0]=delh_A[0][1];

  delh_mu[0]=mu[0]+(A[0][0]*(mu[0]-av_x)+A[0][1]*(mu[1]-av_y))*nsample;
  delh_mu[1]=mu[1]+(A[1][0]*(mu[0]-av_x)+A[1][1]*(mu[1]-av_y))*nsample;

  return 0;
}
/***************************/
/*** Molecular evolution ***/
/***************************/
int Molecular_Dynamics(double (&A)[ndim][ndim],double (&mu)[ndim],
                       double& ham_init,double& ham_fin){
  double p_A[ndim][ndim];
  double p_mu[ndim];
  double delh_A[ndim][ndim];
  double delh_mu[ndim];
  double r1,r2;
  
  BoxMuller(r1,r2);
  p_mu[0]=r1;p_mu[1]=r2;
  BoxMuller(r1,r2);
  p_A[0][0]=r1;p_A[1][1]=r2;
  BoxMuller(r1,r2);
  p_A[0][1]=r1/sqrt(2e0);p_A[1][0]=p_A[0][1];
  
  //*** calculate Hamiltonian ***
  ham_init=calc_hamiltonian(p_A,p_mu,A,mu);
  //*** first step of leap frog ***
  for(int idim=0; idim!=ndim; idim++){
    mu[idim]=mu[idim]+p_mu[idim]*0.5e0*dtau_mu;
    for(int jdim=0; jdim!=ndim; jdim++){
      A[idim][jdim]=A[idim][jdim]+p_A[idim][jdim]*0.5e0*dtau_A;
    }
  }
  //*** 2nd, ..., Ntau-th steps ***
  for(int step=1; step!=ntau; step++){    
    calc_delh(A,mu,delh_A,delh_mu);
    for(int idim=0; idim!=ndim; idim++){
      p_mu[idim]=p_mu[idim]-delh_mu[idim]*dtau_mu;
      for(int jdim=0; jdim!=ndim; jdim++){    
        p_A[idim][jdim]=p_A[idim][jdim]-delh_A[idim][jdim]*dtau_A;     
      }
    }
    for(int idim=0; idim!=ndim; idim++){
      mu[idim]=mu[idim]+p_mu[idim]*dtau_mu;
      for(int jdim=0; jdim!=ndim; jdim++){
        A[idim][jdim]=A[idim][jdim]+p_A[idim][jdim]*dtau_A;
      }
    }
  }
  //*** last step of leap frog ***
  calc_delh(A,mu,delh_A,delh_mu);
  for(int idim=0; idim!=ndim; idim++){
    p_mu[idim]=p_mu[idim]-delh_mu[idim]*dtau_mu;
    for(int jdim=0; jdim!=ndim; jdim++){
      p_A[idim][jdim]=p_A[idim][jdim]-delh_A[idim][jdim]*dtau_A;
    }
  }
  for(int idim=0; idim!=ndim; idim++){
    mu[idim]=mu[idim]+p_mu[idim]*0.5e0*dtau_mu;
    for(int jdim=0; jdim!=ndim; jdim++){
      A[idim][jdim]=A[idim][jdim]+p_A[idim][jdim]*0.5e0*dtau_A;
    }
  }
  //*** calculate Hamiltonian again ***
  ham_fin=calc_hamiltonian(p_A,p_mu,A,mu);
  
  return 0;
}

int main(){
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  //NEVER choose a singular initial condition.
  double A[ndim][ndim];
  A[0][0]=1e0;A[1][1]=1e0;
  A[0][1]=0e0;
  A[1][0]=A[0][1];

  double mu[ndim];
  mu[0]=0e0;mu[1]=0e0;

  /*****************/
  /*** Main part ***/
  /*****************/
  int naccept=0;//counter for the number of acceptance
  for(int iter=0; iter!=niter; iter++){
    double backup_A[ndim][ndim],backup_mu[ndim];

    for(int idim=0; idim!=ndim; idim++){
      backup_mu[idim]=mu[idim];    
      for(int jdim=0; jdim!=ndim; jdim++){
        backup_A[idim][jdim]=A[idim][jdim];
      }
    }
    double ham_init,ham_fin;
    Molecular_Dynamics(A,mu,ham_init,ham_fin);
    double metropolis = (double)rand()/RAND_MAX;
    if(exp(ham_init-ham_fin) > metropolis){
      //accept
      naccept=naccept+1;
    }else{
      //reject
      for(int idim=0; idim!=ndim; idim++){
        mu[idim]=backup_mu[idim];
        for(int jdim=0; jdim!=ndim; jdim++){
          A[idim][jdim]=backup_A[idim][jdim];
        }
      }
    }
    /*******************/
    /*** data output ***/
    /*******************/
    // output 
    if((iter+1)%nskip==0){
      std::cout << std::fixed << std::setprecision(6) 
        << A[0][0] << "   " 
        << A[1][1] << "   " 
        << A[0][1] << "   " 
        << mu[0] << "   " 
        << mu[1] << "   " 
        << (double) naccept/(iter+1) 
        << std::endl;
    }
  } 
  return 0;
}
  


