#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int nbeta=200;
const int niter=5000;
const double dbeta=0.5e0;
const int ncity=100;
const int ninit=1;// 0 -> read "100_cities.txt"; 1 -> random config.

double calc_distance(double x[2][ncity+1],int ordering[nbeta][ncity+1],int ibeta) {
  double distance = 0e0;
  double r1,r2;
  for(int icity=0;icity<ncity;icity++){
    int i=ordering[ibeta][icity];
    int j=ordering[ibeta][icity+1];
    r1=(x[0][i]-x[0][j]);
    r2=(x[1][i]-x[1][j]);
    distance=distance+sqrt(r1*r1+r2*r2);
  }
  //x[i][0]=x[i][ncity]
  //ordering[i][0]=0
  //ordering[i][ncity]=ncity 
  //(the first and the last are the same.)
  return distance;
}


int main(void){
  int ordering[nbeta][ncity+1];
  //ordering[i][0]=0
  //ordering[i][ncity]=ncity    
  double beta[nbeta];
  double x[2][ncity+1]; //location of the cities on 2d plane. x[i][0]=x[i][ncity].
  int naccept[nbeta];
  double minimum_distance=100e0;
  int minimum_ordering[ncity+1];
  srand((unsigned)time(NULL)); 
 
  /*********************************/
  /* Set the initial configuration */
  /*********************************/      
  for(int ibeta=0;ibeta<nbeta;ibeta++){
    for(int icity=0;icity<ncity+1;icity++){
      ordering[ibeta][icity]=icity;
    }
    beta[ibeta]=(double)(ibeta+1)*dbeta;
  }
  if(ninit==0){
    FILE *file;
    double read1,read2;
    file=fopen("100_cities.txt","r");
    for(int icity=0;icity<ncity;icity++){
      fscanf(file,"%lf",&read1);
      fscanf(file,"%lf",&read2);
      x[0][icity]=read1;
      x[1][icity]=read2;
    }
    fclose(file);
  }else if(ninit==1){
    FILE *input_config;
    double read1,read2;
    int read_ordering;
    input_config=fopen("input_config.txt","r");  
    for(int icity=0;icity<ncity;icity++){
      fscanf(input_config,"%lf",&read1);
      fscanf(input_config,"%lf",&read2);
	x[0][icity]=read1;
	x[1][icity]=read2;
    }
    for(int ibeta=0;ibeta<nbeta;ibeta++){
      for(int icity=0;icity<ncity+1;icity++){
	fscanf(input_config,"%i",&read_ordering);
	ordering[ibeta][icity]=read_ordering;
      }
    }
    
  }else if(ninit==2){
    for(int icity=0;icity<ncity;icity++){
      x[0][icity]=(double)rand()/RAND_MAX;
      x[1][icity]=(double)rand()/RAND_MAX;
    }
  }
  x[0][ncity]=x[0][0];
  x[1][ncity]=x[1][0];
  /*************/
  /* Main loop */
  /*************/
  FILE *outputfile = fopen("output.txt", "w");
  for(int iter=1;iter<niter+1;iter++){
    for(int ibeta=0;ibeta<nbeta;ibeta++){
      int info_kl=1;
      int k,l;
      while(info_kl==1){
        k=(int)((double)rand()/RAND_MAX*(ncity-1))+1;
        l=(int)((double)rand()/RAND_MAX*(ncity-1))+1;
        if(k!=l){info_kl=0;}
      }
      /*******************************/
      /* Metropolis for each replica */
      /*******************************/
      double action_init=calc_distance(x,ordering,ibeta)*beta[ibeta];
      int temp=ordering[ibeta][k];
      ordering[ibeta][k]=ordering[ibeta][l];
      ordering[ibeta][l]=temp;
      double action_fin=calc_distance(x,ordering,ibeta)*beta[ibeta];
      /*******************/
      /* Metropolis test */
      /*******************/
      double metropolis = (double)rand()/RAND_MAX;    
      if(exp(action_init-action_fin) > metropolis){
         /* accept */
         naccept[ibeta]=naccept[ibeta]+1;
      }else{ 
         /* reject */
         temp=ordering[ibeta][k];
       ordering[ibeta][k]=ordering[ibeta][l];
       ordering[ibeta][l]=temp;
       }}
    /*********************/
    /* Exchange replicas */
    /*********************/    
    for(int ibeta=0;ibeta<nbeta-1;ibeta++){
      double action_init=calc_distance(x,ordering,ibeta)*beta[ibeta]+calc_distance(x,ordering,ibeta+1)*beta[ibeta+1];
      double action_fin=calc_distance(x,ordering,ibeta)*beta[ibeta+1]+calc_distance(x,ordering,ibeta+1)*beta[ibeta];
      /*******************/
      /* Metropolis test */
      /*******************/
      double metropolis = (double)rand()/RAND_MAX;
      if(exp(action_init-action_fin) > metropolis){
        /* accept = exchange */
        for(int icity=0;icity<ncity;icity++){        
          int backup_ordering=ordering[ibeta][icity];
          ordering[ibeta][icity]=ordering[ibeta+1][icity];
          ordering[ibeta+1][icity]=backup_ordering;
        }
      }
    } 
    
    /***************/
    /* data output */
    /***************/    
    if(iter%500==0){
      double distance_200=calc_distance(x,ordering,199);
      double distance_100=calc_distance(x,ordering,99);
      if(distance_200<minimum_distance){
        minimum_distance=distance_200;
        for(int icity=0;icity<ncity+1;icity++){
          minimum_ordering[icity]=ordering[nbeta-1][icity];
        }
      }
      printf("%i   %lf     %lf\n",iter,distance_200,minimum_distance);    
      fprintf(outputfile,"%i   %lf     %lf\n",iter,distance_200,minimum_distance);    
    }
  }
  for(int icity=0;icity<ncity+1;icity++){
    printf("%lf     %lf\n",x[0][ordering[nbeta-1][icity]],x[1][ordering[nbeta-1][icity]]);
    fprintf(outputfile,"%lf     %lf\n",x[0][ordering[nbeta-1][icity]],x[1][ordering[nbeta-1][icity]]);
  }
  fclose(outputfile);

  FILE *output_config = fopen("output_config.txt", "w");
  for(int icity=0;icity<ncity;icity++){
    fprintf(output_config,"%lf   %lf\n",x[0][icity],x[1][icity]);
  }
  for(int ibeta=0;ibeta<nbeta;ibeta++){
    for(int icity=0;icity<ncity+1;icity++){
      fprintf(output_config,"%i  ",ordering[ibeta][icity]);
    }
  }
  fclose(output_config);
}
