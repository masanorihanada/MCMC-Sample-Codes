#include <iostream>
#include <cmath>
#include <complex>
#include<fstream>
const int nmat=10;
const int niter=100;
const int ninit=0; //ninit=1 -> new config; ninit=0 -> old config
const int ntau=20;
const double dtau=0.05e0;
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
double calc_action(const std::complex<double> phi[nmat][nmat]){
  
  std::complex<double> phi2[nmat][nmat];
  
  //phi2=phi*phi
  for(int imat=0;imat!=nmat;imat++){
    for(int jmat=0;jmat!=nmat;jmat++){
      for(int kmat=0;kmat!=nmat;kmat++){
        phi2[imat][jmat]=phi2[imat][jmat]+phi[imat][kmat]*phi[kmat][jmat];
      }
    }
  }
  
  double action=0e0;
  //Tr phi^2 term
  for(int imat=0; imat!=nmat; ++imat){
    action=action+0.5e0*phi2[imat][imat].real();
  }
  //Tr phi^4 term
  for(int imat=0;imat!=nmat;imat++){
    for(int jmat=0;jmat!=nmat;jmat++){
      action=action+0.25e0*(phi2[imat][jmat]*phi2[jmat][imat]).real();
    }
  }
  //overall normalization 
  action=action*nmat;
  
  return action;
}
/**************************************/
/*** Calculation of the Hamiltonian ***/
/**************************************/
double calc_hamiltonian(const std::complex<double> phi[nmat][nmat],
                        const std::complex<double> P_phi[nmat][nmat]){
  double ham=calc_action(phi);
  
  for(int imat=0;imat!=nmat;imat++){
    for(int jmat=0;jmat!=nmat;jmat++){
      ham=ham+0.5e0*(P_phi[imat][jmat]*P_phi[jmat][imat]).real();
    }
  }
  
  return ham;
}
/********************************/
/*** Calculation of the force ***/
/********************************/
int calc_delh(std::complex<double> (&delh)[nmat][nmat],
               const std::complex<double> phi[nmat][nmat]){
  
  std::complex<double> phi2[nmat][nmat],phi3[nmat][nmat];
  
  //*** phi2=phi*phi, phi3=phi*phi*phi ***
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      phi2[imat][jmat]=std::complex<double>(0e0,0e0);
      phi3[imat][jmat]=std::complex<double>(0e0,0e0);
    }
  }
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      for(int kmat=0; kmat!=nmat; kmat++){
        phi2[imat][jmat]=phi2[imat][jmat]+phi[imat][kmat]*phi[kmat][jmat];
      }
    }
  }
  
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      for(int kmat=0; kmat!=nmat; kmat++){
        phi3[imat][jmat]=phi3[imat][jmat]+phi2[imat][kmat]*phi[kmat][jmat];
      }
    }
  }
  //!*** delh=dH/dphi ***
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      delh[imat][jmat]=(phi[imat][jmat]+phi3[imat][jmat])*(double)nmat;
    }
  }  
  return 0;
}
/***************************/
/*** Molecular evolution ***/
/***************************/
int Molecular_Dynamics(std::complex<double> (&phi)[nmat][nmat],
                       double& ham_init,double& ham_fin){
  
  std::complex<double> P_phi[nmat][nmat];
  std::complex<double> delh[nmat][nmat];
  
  for(int imat=0; imat!=nmat-1; imat++){
    for(int jmat=imat; jmat!=nmat; jmat++){
      double r1,r2;
      BoxMuller(r1,r2);
      P_phi[imat][jmat]=r1/sqrt(2e0)+r2/sqrt(2e0)*std::complex<double>(0e0,1e0);
      P_phi[jmat][imat]=r1/sqrt(2e0)-r2/sqrt(2e0)*std::complex<double>(0e0,1e0);
    }
  }
  for(int imat=0; imat!=nmat; imat++){
    double r1,r2;
    BoxMuller(r1,r2);
    P_phi[imat][imat]=r1;
  }
  
  //*** calculate Hamiltonian ***
  ham_init=calc_hamiltonian(phi,P_phi);
  //*** first step of leap frog ***
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      phi[imat][jmat]=phi[imat][jmat]+P_phi[imat][jmat]*0.5e0*dtau;
    }
  }
  //*** 2nd, ..., Ntau-th steps ***
  for(int step=1; step!=ntau; step++){
    
    calc_delh(delh,phi);
    for(int imat=0; imat!=nmat; imat++){
      for(int jmat=0; jmat!=nmat; jmat++){
        P_phi[imat][jmat]=P_phi[imat][jmat]-delh[imat][jmat]*dtau;
        phi[imat][jmat]=phi[imat][jmat]+P_phi[imat][jmat]*dtau;
      }
    }
  }
  //*** last step of leap frog ***
  calc_delh(delh,phi);
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      P_phi[imat][jmat]=P_phi[imat][jmat]-delh[imat][jmat]*dtau;
      phi[imat][jmat]=phi[imat][jmat]+P_phi[imat][jmat]*0.5e0*dtau;
    }
  }
  //*** calculate Hamiltonian ***
  ham_fin=calc_hamiltonian(phi,P_phi);
  
  return 0;
}

int main()
{
  srand((unsigned)time(NULL));
  /*********************************/
  /* Set the initial configuration */
  /*********************************/
  std::complex<double> phi[nmat][nmat]; 
  if(ninit==1){
    for(int imat=0; imat!=nmat; imat++){
      for(int jmat=0; jmat!=nmat; jmat++){
        phi[imat][jmat]=std::complex<double>(0e0,0e0);
      }
    }
  }
  if(ninit==0){
    //read old data
    std::ifstream configfile("configuration.dat");
    
    for(int imat=0; imat!=nmat; imat++){
      for(int jmat=0; jmat!=nmat; jmat++){
        configfile >> phi[imat][jmat];
      }
    }
    configfile.close();
  }
  
  double sum_action=0e0;
  /*****************/
  /*** Main part ***/
  /*****************/
  std::ofstream outputfile("output.txt");
  int naccept=0;//counter for the number of acceptance
  for(int iter=0; iter!=niter; iter++){
    std::complex<double> backup_phi[nmat][nmat];
    for(int imat=0; imat!=nmat; imat++){
      for(int jmat=0; jmat!=nmat; jmat++){
        backup_phi[imat][jmat]=phi[imat][jmat];
      }
    }
    double ham_init,ham_fin;
    Molecular_Dynamics(phi,ham_init,ham_fin);
    double metropolis = (double)rand()/RAND_MAX;
    if(exp(ham_init-ham_fin) > metropolis){
      //accept
      naccept=naccept+1;
    }else{
      //reject
      for(int imat=0; imat!=nmat; imat++){
        for(int jmat=0; jmat!=nmat; jmat++){
          phi[imat][jmat]=backup_phi[imat][jmat];
        }
      }
    }
    
    /*******************/
    /*** data output ***/
    /*******************/
    double action=calc_action(phi);
    sum_action=sum_action+action;
    
    std::cout << std::fixed << std::setprecision(6) 
      << sum_action/((iter+1)*nmat*nmat) << "   "
      << ((double)naccept)/(iter+1) 
      << std::endl;
     
    outputfile << std::fixed << std::setprecision(6) 
      << sum_action/((iter+1)*nmat*nmat) << "   "
      << ((double)naccept)/(iter+1) 
      << std::endl;
    
  }
  outputfile.close();
  
  //save final configuration
  std::ofstream configfile("configuration.dat");
  for(int imat=0; imat!=nmat; imat++){
    for(int jmat=0; jmat!=nmat; jmat++){
      configfile << phi[imat][jmat] << std::endl;
    }
  }
  configfile.close();
  
  return 0;
}
