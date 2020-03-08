#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void){
  int niter=100000;
  srand((unsigned)time(NULL)); 
  
  double pi=asin(1e0)*2e0;
  double sum_z=0e0;
  int n_in=0;
  /*************/
  /* Main loop */
  /*************/
  for(int iter=1;iter<niter+1;iter++){

    double x = (double)rand()/RAND_MAX;
    double y = (double)rand()/RAND_MAX;
    if(x*x+y*y < 1e0){
      n_in=n_in+1;
      double z=sqrt(1e0-x*x-y*y);
      sum_z=sum_z+z;
    }
    printf("%d   %.10f\n",iter,sum_z/n_in*2e0*pi);}
}
