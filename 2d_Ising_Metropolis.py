import numpy as np
import os
import sys
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default

NITER = 4096000
NX = 64 #number of sites along x-direction
NY = 64 #number of sites along y-direction
COUPLING_J = 1.0
COUPLING_h = 0.1
TEMPERATURE = 5.0
NSKIP = 4096 #Frequency of measurement
NCONFIG = 0 #0 -> read 'input_config.txt'; 1 -> all up; -1 -> all down
#################################
### Calculation of the action ###
#################################
# original version 
# TOO SLOW!! (not vectorized)
#def calc_action(spin,COUPLING_J,COUPLING_h,TEMPERATURE):
#    action=0.0
#    sum1=np.sum(spin)
#    sum2=0
#    for ix in range(NX):
#        ixp1 = (ix+1)%NX #ixp1=ix+1; be careful about the boundary condition.
#        for iy in range(NY):
#            iyp1 = (iy+1)%NY #iyp1=iy+1; be careful about the boundary condition.
#            #sum1=sum1+spin[ix,iy]
#            sum2=sum2+spin[ix,iy]*spin[ixp1,iy]+spin[ix,iy]*spin[ix,iyp1]
#    action=(sum2*COUPLING_J+sum1*COUPLING_h)/TEMPERATURE*(-1.0)
#    return action
def calc_action(spin,COUPLING_J,COUPLING_h,TEMPERATURE):
    action=0.0
    spin2=np.array(spin)
    sum1=np.sum(spin2)

    spin_x=np.roll(spin2, 1, axis=1) # shift to x-direction with pbc
    spin_y=np.roll(spin2, 1, axis=0) # shift to y-direction with pbc
    sum2=np.sum(spin*spin_x + spin*spin_y)

    action=(sum2*COUPLING_J+sum1*COUPLING_h)/TEMPERATURE*(-1.0)
    return action
###############################################
### Calculation of the change of the action ###
###  when the spin at (ix,iy) is flipped    ###
###############################################
def calc_action_change(spin,COUPLING_J,COUPLING_h,TEMPERATURE,ix,iy):
    action_change=0.0

    ixp1=(ix+1)%NX #ixp1=ix+1; be careful about the boundary condition.
    iyp1=(iy+1)%NY #iyp1=iy+1; be careful about the boundary condition.
    ixm1=(ix-1+NX)%NX #ixm1=ix-1; be careful about the boundary condition.
    iym1=(iy-1+NY)%NY #iym1=iy-1; be careful about the boundary condition.

    sum1_change=2*spin[ix,iy]
    sum2_change=2*spin[ix,iy]*spin[ixp1,iy]+2*spin[ix,iy]*spin[ix,iyp1]+2*spin[ix,iy]*spin[ixm1,iy]+2*spin[ix,iy]*spin[ix,iym1]
      
    action_change=(sum2_change*COUPLING_J+sum1_change*COUPLING_h)/TEMPERATURE
    
    return action_change
#####################################
### Calculation of the total spin ###
#####################################
def calc_total_spin(spin):
    return np.sum(np.array(spin))
#####################################
### Set the initial configuration ###
#####################################
if(NCONFIG==1):
    spin = np.ones((NX,NY))
elif(NCONFIG==-1):
    spin = np.ones((NX,NY))
    spin=spin*(-1)
elif(NCONFIG==0):
    if(os.path.exists('input_config.txt')):
        spin = np.empty((NX,NY))
        for read in open('input_config.txt').readlines():
            read = read[:-2].split(' ')
            ix = int(read[0])
            iy = int(read[1])
            spin[ix,iy] = read[2]
    else:
        print('no input configuration')
        sys.exit()
#################
### Main part ###
#################
path_output = 'output.txt'
with open(path_output, mode='w') as output:

    naccept = 0 #counter for the number of acceptance
    for iter in range(NITER):
        #choose a point randomly.
        rand_site = np.random.uniform(0,NX*NY)
        ix :int = int(rand_site/NX)
        iy :int = int(rand_site%NY)
        action_change = calc_action_change(spin,COUPLING_J,COUPLING_h,TEMPERATURE,ix,iy)
        metropolis = np.random.uniform(0,1)
        if(np.exp(-action_change) > metropolis):
            #accept
            spin[ix,iy]=-spin[ix,iy]
            naccept=naccept+1;
        #else:
            #reject

        total_spin = calc_total_spin(spin)
        energy = calc_action(spin,COUPLING_J,COUPLING_h,TEMPERATURE)*TEMPERATURE
        ###################
        ### data output ###
        ###################
        if((iter+1)%NSKIP == 0):
            print(total_spin,energy,naccept/(iter+1),file=output)
            print(total_spin,energy,naccept/(iter+1))
            
######################################
### 2D plot of final configuration ###
######################################
plt.figure()
plt.imshow(spin,interpolation='nearest',vmin=0,vmax=1,cmap='jet')
plt.show()
            
#########################
### save final config ###
#########################
path_output_config = 'output_config.txt'
with open(path_output_config, mode='w') as output_config:
    for ix in range(NX):
        for iy in range(NY):
            print(ix,iy,spin[ix,iy], file=output_config)

