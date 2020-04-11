/*******************************************/
/*** 2d Ising model with Wolff algorithm ***/
/*******************************************/
#include <iostream>
#include <cmath>
#include<fstream>
const long int niter=1000000;
const int nx=64;   //number of sites along x-direction
const int ny=64;   //number of sites along y-direction 
const double coupling_J=1.0e0;
const double coupling_h=0.1e0;
const double temperature=5e0;
const int nskip=100;    //Frequency of measurement
const int nconfig=10000000;    //Frequency of saving configs; 0->don't save
/*********************************/
/*** Calculation of the action ***/
/*********************************/
double calc_action(const int spin[nx][ny],const double coupling_J,const double coupling_h,const double temperature){
  
  double action=0e0;
  int sum1=0; 
  int sum2=0;
  for(int ix=0; ix!=nx; ix++){
      int ixp1=(ix+1)%nx; //ixp1=ix+1; be careful about the boundary condition.
      for(int iy=0; iy!=ny; iy++){
      int iyp1=(iy+1)%ny; //iyp1=iy+1; be careful about the boundary condition. condition.
      sum1=sum1+spin[ix][iy];
      sum2=sum2+spin[ix][iy]*spin[ixp1][iy]+spin[ix][iy]*spin[ix][iyp1];
    }
  }
  action=(sum2*coupling_J+sum1*coupling_h)/temperature*(-1e0);
  
  return action;
}
/*************************************/
/*** Calculation of the total spin ***/
/*************************************/
int calc_total_spin(const int spin[nx][ny]){

  int total_spin=0;

  for(int ix=0; ix!=nx; ix++){
     for(int iy=0; iy!=ny; iy++){
       total_spin = total_spin + spin[ix][iy];
    }
  }

  return total_spin;
}
/*****************************/
/*** Construct the cluster ***/
/*****************************/
int make_cluster(const int spin[nx][ny],const double coupling_J,
  const double temperature,int& n_cluster,int (&i_cluster)[nx*ny][2]){
  
  int in_or_out[nx][ny];
  for(int ix=0; ix!=nx; ix++){
    for(int iy=0; iy!=ny; iy++){
      in_or_out[ix][iy]=1;
    }
  }
  //in_or_out[ix][iy] = 1 -> not in the cluster; 0 -> in the cluster. 
  //choose a point randomly.                                                                         
  double rand_site = (double)rand()/RAND_MAX;
  rand_site=rand_site*nx*ny;
  int ix=(int)rand_site/ny;
  int iy=(int)rand_site%ny;
  in_or_out[ix][iy]=0;
  i_cluster[0][0]=ix;
  i_cluster[0][1]=iy;
  int spin_cluster=spin[ix][iy];
  n_cluster=1;
  double probability=1e0-exp(-2e0*coupling_J/temperature);
  int k=0;
  while(k < n_cluster){
    ix=i_cluster[k][0];
    iy=i_cluster[k][1];
    int ixp1=(ix+1)%nx; //ixp1=ix+1; be careful about the boundary condition.                                                                                      
    int iyp1=(iy+1)%ny; //iyp1=iy+1; be careful about the boundary condition.                                                                                      
    int ixm1=(ix-1+nx)%nx; //ixm1=ix-1; be careful about the boundary condition.                                                                                   
    int iym1=(iy-1+ny)%ny; //iym1=iy-1; be careful about the boundary condition.                                                                                   

    if(spin[ixp1][iy]==spin_cluster){
      if(in_or_out[ixp1][iy]==1){
        if((double)rand()/RAND_MAX < probability){
          i_cluster[n_cluster][0]=ixp1;
          i_cluster[n_cluster][1]=iy;
          n_cluster=n_cluster+1;
          in_or_out[ixp1][iy]=0;
        }
      }
    }
    if(spin[ix][iyp1]==spin_cluster){
      if(in_or_out[ix][iyp1]==1){
        if((double)rand()/RAND_MAX < probability){
          i_cluster[n_cluster][0]=ix;
          i_cluster[n_cluster][1]=iyp1;
          n_cluster=n_cluster+1;
          in_or_out[ix][iyp1]=0;
        }
      } 
    }
    if(spin[ixm1][iy]==spin_cluster){
      if(in_or_out[ixm1][iy]==1){
        if((double)rand()/RAND_MAX < probability){
          i_cluster[n_cluster][0]=ixm1;
          i_cluster[n_cluster][1]=iy;
          n_cluster=n_cluster+1;
          in_or_out[ixm1][iy]=0;
        }
      } 
    }
    if(spin[ix][iym1]==spin_cluster){
      if(in_or_out[ix][iym1]==1){
        if((double)rand()/RAND_MAX < probability){
          i_cluster[n_cluster][0]=ix;
          i_cluster[n_cluster][1]=iym1;
          n_cluster=n_cluster+1;
          in_or_out[ix][iym1]=0;
        }
      }
    }
    k=k+1;
  }
  return spin_cluster;
}
/************/
/*** Main ***/
/************/
int main()
{
  int spin[nx][ny];
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  for(int ix=0; ix!=nx; ix++){
    for(int iy=0; iy!=ny; iy++){
      spin[ix][iy]=1;
    }
  }
  /*****************/
  /*** Main part ***/
  /*****************/
  std::ofstream outputfile("output.txt");
  std::ofstream outputconfig("output_config.txt");

  for(long int iter=0; iter!=niter; iter++){
    int n_cluster;
    int i_cluster[nx*ny][2];
    int spin_cluster = make_cluster(spin,coupling_J,temperature,n_cluster,i_cluster);
    double metropolis = (double)rand()/RAND_MAX;
    if(exp(-2e0/temperature*coupling_h*spin_cluster*n_cluster) > metropolis){  
      for(int k=0; k!=n_cluster; k++){
        int ix=i_cluster[k][0];
        int iy=i_cluster[k][1];
        spin[ix][iy]=-spin[ix][iy];
      }
    }
    
    int total_spin = calc_total_spin(spin);
    double energy = calc_action(spin,coupling_J,coupling_h,temperature)*temperature;
    /*******************/
    /*** data output ***/
    /*******************/    
    if((iter+1)%nskip == 0){
      std::cout << std::fixed << std::setprecision(4) 
        << total_spin << "  " 
        << energy  << std::endl;
      outputfile << std::fixed << std::setprecision(4) 
        << total_spin << "  " 
        << energy  << std::endl;
    }
    /*********************/
    /*** config output ***/
    /*********************/
    if(nconfig > 0){
      if((iter+1)%nconfig == 0){
        for(int ix=0; ix!=nx; ix++){
          for(int iy=0; iy!=ny; iy++){
            outputconfig  << ix << ' ' << iy << ' ' << spin[ix][iy] << ' '<< std::endl;
          }
        }
      }
    }
  }
  outputfile.close();
  outputconfig.close();
  return 0;
}
  


