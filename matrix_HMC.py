import numpy as np
import os
import sys
import copy
# seed of random numbers is the system time by default
# Numpy has Gaussian random number generator, so we don't use Box-Muller
NMAT=10
NITER=1000
NINIT=0 #ninit=1 -> new config; ninit=0 -> old config
NTAU=20
DTAU=0.05
##################################
### Calculation of the action  ###
##################################
#When you change the action, you should also change dH/dx, specified in "calc_delh".
def calc_action(phi):
    phi2 = np.dot(phi,phi)
    phi4 = np.dot(phi2,phi2)
    action = (np.trace(0.5*phi2+0.25*phi4)).real * NMAT
    return action
######################################
### Calculation of the Hamiltonian ###
######################################
def calc_hamiltonian(phi,P_phi):
    ham = calc_action(phi)
    ham += 0.5 * (np.trace(np.dot(P_phi,P_phi))).real
    return ham
################################
### Calculation of the force ###
################################
def calc_delh(phi):
    phi2 = np.dot(phi,phi)
    phi3 = np.dot(phi,phi2)
    delh = (phi + phi3) * NMAT
    return delh
############################
### Molecular evolution  ###
############################
def Molecular_Dynamics(phi):
    # Mementum P_phi is chosen to be Gaussian
    P_phi = np.zeros((NMAT,NMAT),dtype = 'complex_')
    for imat in range(NMAT-1):
        for jmat in range(imat+1,NMAT):
            r1 = np.random.randn()
            r2 = np.random.randn()
            P_phi[imat,jmat] = r1/np.sqrt(2.0) + r2/np.sqrt(2.0) * 1j
            P_phi[jmat,imat] = r1/np.sqrt(2.0) - r2/np.sqrt(2.0) * 1j
    for imat in range(NMAT):
            P_phi[imat,imat] = np.random.randn()
    # calculate Hamiltonian
    ham_init = calc_hamiltonian(phi,P_phi)
    ## first step of leap frog
    phi += P_phi * 0.5 * DTAU
    # 2nd, ..., Ntau-th steps
    for step in range(NTAU):
        delh=calc_delh(phi)
        P_phi += -delh * DTAU
        phi += P_phi * DTAU
    # last step of leap frog
    delh=calc_delh(phi)
    P_phi += -delh * DTAU
    phi += P_phi * 0.5 * DTAU
    # calculate Hamiltonian again
    ham_fin = calc_hamiltonian(phi,P_phi)
    return phi,ham_init,ham_fin
#####################################
### Set the initial configuration ###
#####################################
if(NINIT==1):
    phi = np.zeros((NMAT,NMAT),dtype = 'complex_')
if(NINIT==0):#read old data
    if(os.path.exists('config.dat')):
        phi = np.zeros((NMAT,NMAT),dtype = 'complex_')
        for read in open('config.dat').readlines():
            read = read[0:-3].split(' ')
            imat = int(read[0])
            jmat = int(read[1])

            if(imat == jmat):
                phi[imat,jmat] = float(read[2])
            else:
                phi[imat,jmat] = float(read[2]) + float(read[3])*1j
    else:
        print('no input configuration')
        sys.exit()
#################
### Main part ###
#################
path_output = 'output.txt'
with open(path_output, mode='w') as output:
    sum_action=0.0
    naccept = 0#counter for the number of acceptance
    for iter in range(NITER):
        backup_phi = copy.deepcopy(phi)
        phi,ham_init,ham_fin = Molecular_Dynamics(phi)
        metropolis = np.random.uniform(0,1)
        if(np.exp(ham_init-ham_fin) > metropolis):#accept
            naccept += 1
        else:#reject
            phi = copy.deepcopy(backup_phi)
        ###################
        ### data output ###
        ###################
        action = calc_action(phi)
        sum_action += action
        if((iter+1)%10 == 0):#output the results every ten steps.
            print(sum_action/((iter+1)*NMAT*NMAT),naccept/(iter+1))
            print(sum_action/((iter+1)*NMAT*NMAT),naccept/(iter+1), file=output)


#########################
### save final config ###
#########################
path_output_config = 'config.dat'
with open(path_output_config, mode='w') as output_config:
    for imat in range(NMAT):
        for jmat in range(NMAT):
            temp_re = (phi[imat,jmat]).real
            temp_im = (phi[imat,jmat]).imag
            print(imat,jmat,temp_re,temp_im, file=output_config)


