import numpy as np
# seed of random numbers is the system time by default

NITER = 100000
sum_y=0.0
#################
### Main Loop ###
#################
for iter in range(NITER):
    x = np.random.uniform(0, 1)
    y = np.sqrt(1.0-x*x)
    #sum_y = sum_y + y
    sum_y += y
    print(iter+1,sum_y/(iter+1))


