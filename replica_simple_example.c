#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int nbeta=2000;
const int niter=1000000;
const double step_size=0.1e0;
const double dbeta=0.5e0;

double calc_f(double x) {
  double fx = (x-1e0)*(x-1e0)*((x+1e0)*(x+1e0)+0.01e0);

  return fx;
}

int main(void){
  int naccept[nbeta];
  double x[nbeta]; 
  double beta[nbeta];

  srand((unsigned)time(NULL)); 

  /*********************************/
  /* Set the initial configuration */
  /*********************************/      
  for(int ibeta=0;ibeta<nbeta;ibeta++){
    x[ibeta]=0e0;
    beta[ibeta]=(double)(ibeta+1)*dbeta;
    naccept[ibeta]=0e0;
  }
  /*************/
  /* Main loop */
  /*************/
  for(int iter=1;iter<niter+1;iter++){
    for(int ibeta=0;ibeta<nbeta;ibeta++){
      double backup_x=x[ibeta];    
      double action_init=calc_f(x[ibeta])*beta[ibeta];
      
      double dx = (double)rand()/RAND_MAX;
      dx=(dx-0.5e0)*step_size*2e0;
      x[ibeta]=x[ibeta]+dx;
      
      double action_fin=calc_f(x[ibeta])*beta[ibeta];
      /*******************/
      /* Metropolis test */
      /*******************/
      double metropolis = (double)rand()/RAND_MAX;    
      if(exp(action_init-action_fin) > metropolis)
        /* accept */
        naccept[ibeta]=naccept[ibeta]+1;
      else 
        /* reject */
        x[ibeta]=backup_x;
    }
    for(int ibeta=0;ibeta<nbeta-1;ibeta++){
      double action_init=calc_f(x[ibeta])*beta[ibeta]+calc_f(x[ibeta+1])*beta[ibeta+1];
      double action_fin=calc_f(x[ibeta])*beta[ibeta+1]+calc_f(x[ibeta+1])*beta[ibeta];
      /*******************/
      /* Metropolis test */
      /*******************/
      double metropolis = (double)rand()/RAND_MAX;
      if(exp(action_init-action_fin) > metropolis){
        /* accept = exchange */
        double backup_x=x[ibeta];
        x[ibeta]=x[ibeta+1];
        x[ibeta+1]=backup_x;
      }
    }

      /***************/
      /* data output */
      /***************/
    printf("%f    %f    %f\n",x[19],x[199],x[1999]);}
}
