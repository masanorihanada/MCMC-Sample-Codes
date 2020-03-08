#include <iostream>
#include <cmath>
#include<fstream>
const int niter=10000;
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
int main()
{
  double A[3][3];
  A[0][0]=1e0;A[1][1]=2e0;A[2][2]=2e0;
  A[0][1]=1e0;A[0][2]=1e0;A[1][2]=1e0;
  A[1][0]=A[0][1];A[2][0]=A[0][2];A[2][1]=A[1][2];

  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  double x=0e0;double y=0e0;double z=0e0;
  /*****************/
  /*** Main part ***/
  /*****************/
  for(int iter=0; iter!=niter; iter++){
    double sigma,mu;
    double r1,r2;    
    //update x
    sigma=1e0/sqrt(A[0][0]);
    mu=-A[0][1]/A[0][0]*y-A[0][2]/A[0][0]*z;
    BoxMuller(r1,r2);
    x=sigma*r1+mu;
    //update y
    sigma=1e0/sqrt(A[1][1]);
    mu=-A[1][0]/A[1][1]*x-A[1][2]/A[1][1]*z;
    BoxMuller(r1,r2);
    y=sigma*r1+mu;
    //update z
    sigma=1e0/sqrt(A[2][2]);
    mu=-A[2][0]/A[2][2]*x-A[2][1]/A[2][2]*y;
    BoxMuller(r1,r2);
    z=sigma*r1+mu;
    // output x,y,z
    if((iter+1)%10==0){
      std::cout << std::fixed << std::setprecision(6) 
        << x << "   " 
        << y << "   " 
        << z << std::endl;
    }
  } 
  
  return 0;
}
  


