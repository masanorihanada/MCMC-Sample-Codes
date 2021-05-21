import numpy as np
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default
# Numpy has Gaussian random number generator, so we don't use Box-Muller

NITER = 500000
NTAU = 10
DTAU = 0.1
##################################
### Calculation of the action  ###
##################################
#When you change the action, you should also change dH/dx, specified in "calc_delh".
def calc_action(x):
    return 0.5*x*x
######################################
### Calculation of the Hamiltonian ###
######################################
def calc_hamiltonian(x,p):
    ham = calc_action(x)
    ham += 0.5*p*p
    return ham
#############################
### Calculation of dH/dx  ###
#############################
def calc_delh(x):
    delh = x
    return delh
#############################
### Molecular evolution  ###
#############################
def Molecular_Dynamics(x):
    # Momentum p is chosen to be Gaussian ##
    p = np.random.randn()
    # calculate Hamiltonian ##
    ham_init = calc_hamiltonian(x,p)
    # first step of leap frog ##
    x += p*0.5*DTAU
    # 2nd, ..., Ntau-th steps ##
    for step in range(1,NTAU):
        delh = calc_delh(x)
        p += -delh*DTAU
        x += p*DTAU
        
    # last step of leap frog ##
    delh = calc_delh(x)
    p += -delh*DTAU
    x += p*0.5e0*DTAU

    # calculate Hamiltonian again ##
    ham_fin=calc_hamiltonian(x,p)
    return x,ham_init,ham_fin
#####################################
### Set the initial configuration ###
#####################################
x = 0
naccept = 0
sum_xx = 0.0 #sum of x^2, useed for <x^2>

# data for plot
data_for_plot = []

#################
### Main Loop ###
#################
for iter in range(NITER):
    backup_x=x
    x,ham_init,ham_fin = Molecular_Dynamics(x)
    metropolis = np.random.uniform(0,1)
    if(np.exp(ham_init-ham_fin) > metropolis):#accept
        naccept += 1
    else:#reject
        x = backup_x
    ###################
    ### data output ###
    ###################
    sum_xx=sum_xx+x*x
    # output x, <x^2>, acceptance
    print(x,sum_xx/(iter+1),naccept/(iter+1))

    # save configuration for plot
    data_for_plot.append(x)

###########################################################
### plot the histogram and compare with analytic result ###
###########################################################
def gauss(x):
    return 1.0/np.sqrt(2.0*np.pi)*np.exp( -x*x/2.0 )

y=np.arange(-4.0,4.0,0.1)
plt.xlabel("$x$")
plt.ylabel("$P(x)$")
plt.plot(y,gauss(y), label="exact")
plt.hist(data_for_plot,bins=80,density=True, label="MCMC")
plt.legend()
plt.show() 
