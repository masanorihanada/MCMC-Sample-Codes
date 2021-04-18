import numpy as np
import copy
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default
# Numpy has Gaussian random number generator, so we don't use Box-Muller

NITER=100000
av_x=-0.0930181
av_y=0.0475899
av_xx=1.06614
av_yy=1.28152
av_xy=-0.504944
NSAMPLE=100
NDIM=2 #don't chage this parameter.
STEP_A=0.1
NSKIP=10
#################################
### Calculation of the action ###
#################################
def calc_action(A,mu):
    action = \
    A[0][0]*((mu[0]-av_x)*(mu[0]-av_x)+av_xx-av_x*av_x)\
    +A[1][1]*((mu[1]-av_y)*(mu[1]-av_y)+av_yy-av_y*av_y)\
    +A[0][1]*((mu[0]-av_x)*(mu[1]-av_y)+av_xy-av_x*av_y)*2.0\
    -np.log(A[0][0]*A[1][1]-A[0][1]*A[1][0])

    action=0.5*action*NSAMPLE

    #add the contribution from the prior distribution
    action += \
    0.5*A[0][0]*A[0][0] + 0.5*A[1][1]*A[1][1] + A[0,1]*A[0,1]\
    + 0.5*mu[0]*mu[0] + 0.5*mu[1]*mu[1]
    return action
  
#####################################
### Set the initial configuration ###
#####################################
#NEVER choose a singular initial condition.
A = [[1.0,0.0],
     [0.0,1.0]]
mu =[0.0,0.0]
A = np.array(A)
mu = np.array(mu)
#################
### Main part ###
#################
naccept = 0 #counter for the number of acceptance

# for plot
A00=[]
A11=[]
A01=[]
mu0=[]
mu1=[]

for iter in range(NITER):
    #############################
    ### Gibbs sampling for mu ###
    #############################
    r = np.random.randn()
    mu[0] = r/np.sqrt(1.0+NSAMPLE*A[0][0])
    mu[0] += -NSAMPLE/(1.0+NSAMPLE*A[0][0])*(A[0][1]*(mu[1]-av_y)-A[0][0]*av_x)
    r = np.random.randn()
    mu[1] = r/np.sqrt(1.0+NSAMPLE*A[1][1])
    mu[1] += -NSAMPLE/(1.0+NSAMPLE*A[1][1])*(A[1][0]*(mu[0]-av_x)-A[1][1]*av_y)
    ########################
    ### Metropolis for A ###
    ########################
    backup_A = copy.deepcopy(A)
    action_init = calc_action(A,mu)
    dx = np.random.uniform(-STEP_A, STEP_A)
    A[0][0] += dx
    dx = np.random.uniform(-STEP_A, STEP_A)
    A[1][1] += dx
    dx = np.random.uniform(-STEP_A, STEP_A)
    A[1][0] += dx
    A[0][1] += dx
    action_fin = calc_action(A,mu)
    #######################
    ### Metropolis test ###
    #######################
    metropolis = np.random.uniform(0,1)
    if(np.exp(action_init-action_fin) > metropolis):#accept
        naccept += 1
    else:#reject
        A = copy.deepcopy(backup_A)
    ###################
    ### data output ###
    ###################
    if(iter%NSKIP == 0):
        print(A[0][0],A[1][1],A[0][1],mu[0],mu[1],naccept/(iter+1))
        if( iter > 1000 ):
          A00.append(A[0][0])
          A11.append(A[1][1])
          A01.append(A[0][1])
          mu0.append(mu[0])
          mu1.append(mu[1])

############
### plot ###
############
plt.figure(figsize=(12,6))
plt.axes().set_aspect('equal')
##
ax1 = plt.subplot(2, 3, 1)
ax1.hist(A00,bins=50,density=True)
ax1.set_xlabel("$A_{11}$")
##
ax2 = plt.subplot(2, 3, 2)
ax2.hist(A11,bins=50,density=True)
ax2.set_xlabel("$A_{22}$")
##
ax3 = plt.subplot(2, 3, 3)
ax3.hist(A01,bins=50,density=True)
ax3.set_xlabel("$A_{12}=A_{21}$")
##
ax4 = plt.subplot(2, 3, 4)
ax4.hist(mu0,bins=50,density=True)
ax4.set_xlabel("$\mu_1$")
##
ax5 = plt.subplot(2, 3, 5)
ax5.hist(mu1,bins=50,density=True)
ax5.set_xlabel("$\mu_2$")
##
plt.show()



