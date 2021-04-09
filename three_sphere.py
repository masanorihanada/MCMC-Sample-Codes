import numpy as np
# seed of random numbers is the system time by default

NITER = 100000

PI = np.pi #numpy knows PI.
sum_z=0.0
n_in = 0
#################
### Main Loop ###
#################
for iter in range(NITER):
    x = np.random.uniform(0, 1)
    y = np.random.uniform(0, 1)
    if(x*x+y*y < 1.0):
        #n_in=n_in+1
        n_in += 1
        z = np.sqrt(1.0-x*x-y*y)
        #sum_z = sum_z + z
        sum_z += z
    if(n_in > 0):#We get error if N_in=0
        print(iter+1,sum_z/n_in*2e0*PI)


