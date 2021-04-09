import numpy as np
# seed of random numbers is the system time by default

NITER = 100000
n_in = 0
#################
### Main Loop ###
#################
for iter in range(NITER):
    x = np.random.uniform(0, 1)
    y = np.random.uniform(0, 1)
    if(x*x+y*y < 1.0):
        #n_in = n_in + 1
        n_in += 1
    print(iter+1,n_in/(iter+1))


