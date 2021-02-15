#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void){
  int niter=1000;
  double step_size=0.03;
  srand((unsigned)time(NULL)); 
  /* Set the initial configuration */      
  double p=0.5e0;
  /* Main part */
  int naccept=0;
  for(int iter=1;iter<niter+1;iter++){
    double backup_p=p;
    
    double action_init=-515e0*log(p)-485e0*log(1-p)+100e0*pow(p-0.9e0,2e0);
    
    double dp = (double)rand()/RAND_MAX;
    dp=(dp-0.5e0)*step_size*2e0;
    p=p+dp;
    
    double action_fin=-515e0*log(p)-485e0*log(1-p)+100e0*pow(p-0.9e0,2e0);
    
    /* Metropolis test */
    double metropolis = (double)rand()/RAND_MAX;
    if((p >= 0e0)&&(p<=1e0)&&(exp(action_init-action_fin) > metropolis))
      /* accept */
      naccept=naccept+1;
    
    else 
      /* reject */
      p=backup_p;
    
    /* data output */	
    printf("%.10f  %f\n",p,(double)naccept/iter);}
  
}
