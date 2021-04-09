import numpy as np
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default

NITER = 1000000
STEP_SIZE = 0.5
#####################################
### Set the initial configuration ###
#####################################
x = 0
naccept = 0

## data for plot
data_for_plot = []

#################
### Main Loop ###
#################
for iter in range(NITER):
    backup_x = x
    action_init = 0.5*x*x
    dx = np.random.uniform(-STEP_SIZE, STEP_SIZE)
    #x = x + dx
    x += dx
    action_fin = 0.5*x*x
    #######################
    ### Metropolis test ###
    #######################
    metropolis = np.random.uniform(0,1)
    if(np.exp(action_init-action_fin) > metropolis):#accept
      #naccept = naccept+1
      naccept += 1
    else:#reject
      x = backup_x
    print(x,naccept/(iter+1))

    ## save configuration for plot
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
