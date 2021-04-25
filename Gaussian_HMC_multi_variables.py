import numpy as np
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default
# Numpy has Gaussian random number generator, so we don't use Box-Muller

NITER = 100000
NTAU = 20
DTAU = 0.5
NDIM = 3

A=[
[1.0,1.0,1.0],
[1.0,2.0,1.0],
[1.0,1.0,2.0]
]
A=np.array(A)
##################################
### Calculation of the action  ###
##################################
#When you change the action, you should also change dH/dx, specified in "calc_delh".
def calc_action(x,A):
    return np.dot(x,np.dot(A,x))*0.5
######################################
### Calculation of the Hamiltonian ###
######################################
def calc_hamiltonian(x,p,A):
    ham = calc_action(x,A)
    ham += 0.5*np.dot(p,p)
    return ham
#############################
### Calculation of dH/dx  ###
#############################
def calc_delh(x):
    delh = np.dot(A,x)
    return delh
############################
### Molecular evolution  ###
############################
def Molecular_Dynamics(x,A):
    # Mementum p is chosen to be Gaussian
    p = np.random.randn(NDIM)
    # calculate Hamiltonian
    ham_init = calc_hamiltonian(x,p,A)
    # first step of leap frog
    x += p*0.5*DTAU
    # 2nd, ..., Ntau-th steps
    for step in range(1,NTAU):
        delh = calc_delh(x)
        p += -delh*DTAU
        x += p*DTAU
    # last step of leap frog
    delh=calc_delh(x)
    p += -delh*DTAU
    x += p*0.5e0*DTAU
    # calculate Hamiltonian again
    ham_fin = calc_hamiltonian(x,p,A)
    return x,ham_init,ham_fin
#####################################
### Set the initial configuration ###
#####################################
x = np.zeros(NDIM)
naccept = 0#counter for the number of acceptance
#################
### Main Loop ###
#################

# for plot
data_for_plot = []

for iter in range(NITER):
    # We have to use "np.copy" here 
    backup_x = np.copy(x)
    x,ham_init,ham_fin = Molecular_Dynamics(x,A)
    metropolis = np.random.uniform(0,1)
    if(np.exp(ham_init-ham_fin) > metropolis):#accept
        naccept = naccept+1
    else:#reject
        # We have to use "np.copy" here 
        x = np.copy(backup_x)
    ###################
    ### data output ###
    ###################
    if((iter+1)%10 == 0):#output the results every ten steps.
        print(*x,naccept/(iter+1))
        data_for_plot.append(x)

data1=[data_for_plot[i][0] for i in range(len(data_for_plot))]
data2=[data_for_plot[i][1] for i in range(len(data_for_plot))]
data3=[data_for_plot[i][2] for i in range(len(data_for_plot))]


############
### plot ###
############
plt.figure(figsize=(15,5))
plt.axes().set_aspect('equal')
##
ax1 = plt.subplot(1, 3, 1)
ax1.scatter(data1,data2,c="black", s=20, marker="+")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_xlim([-8,8])
ax1.set_ylim([-8,8])
ax1.set_title("x-y")
ax1.grid(linestyle = "--")
##
ax2 = plt.subplot(1, 3, 2)
ax2.scatter(data1,data3,c="blue", s=20, marker="+")
ax2.set_xlim([-8,8])
ax2.set_ylim([-8,8])
ax2.set_xlabel("x")
ax2.set_ylabel("z")
ax2.set_title("x-z")
ax2.grid(linestyle = "--")
##
ax3 = plt.subplot(1, 3, 3)
ax3.scatter(data2,data3,c="red", s=20, marker="+")
ax3.set_xlim([-8,8])
ax3.set_ylim([-8,8])
ax3.set_xlabel("y")
ax3.set_ylabel("z")
ax3.set_title("y-z")
ax3.grid(linestyle = "--")
plt.show()

