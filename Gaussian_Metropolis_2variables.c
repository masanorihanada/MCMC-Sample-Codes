#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void){
  int niter=10000;
  double step_size_x=0.5e0;
  double step_size_y=0.5e0;

  srand((unsigned)time(NULL)); 

  /*********************************/
  /* Set the initial configuration */
  /*********************************/      
  double x=0e0;
  double y=0e0;
  int naccept=0;
  /*************/
  /* Main loop */
  /*************/
  for(int iter=1;iter<niter+1;iter++){
    double backup_x=x;
    double backup_y=y;
    double action_init=0.5e0*(x*x+y*y+x*y);
    
    double dx = (double)rand()/RAND_MAX;
    double dy = (double)rand()/RAND_MAX;
    dx=(dx-0.5e0)*step_size_x*2e0;
    dy=(dy-0.5e0)*step_size_y*2e0;
    x=x+dx;
    y=y+dy;
    double action_fin=0.5e0*(x*x+y*y+x*y);
    /*******************/
    /* Metropolis test */
    /*******************/
    double metropolis = (double)rand()/RAND_MAX;
    if(exp(action_init-action_fin) > metropolis){
      /* accept */
      naccept=naccept+1;
    }else{ 
      /* reject */
      x=backup_x;
      y=backup_y;}
    /***************/
    /* data output */
    /***************/
    // output the results every ten steps.	
    if(iter%10==0){
      printf("%.10f   %.10f   %f\n",x,y,(double)naccept/iter);}
  }
}
