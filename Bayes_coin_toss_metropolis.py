import numpy as np
import math
import matplotlib.pyplot as plt

NITER = 10000
step_size = 0.03
p = 0.5

#################
### Main part ###
#################
naccept=0
p_hist=[]
for ite in range(1,NITER):
  backup_p = p
  action_init = -515.0*math.log(p)-485.0*math.log(1-p)+100.0*(p-0.9)**2

  dp = np.random.rand()
  dp = (dp-0.5)*step_size*2.0
  p = p + dp

  if( p<0 or p>1 ): # automatically reject
    p=backup_p
  else:
    action_fin = -515.0*math.log(p)-485.0*math.log(1-p)+100.0*(p-0.9)**2
  
    # metropolis test 
    metropolis = np.random.rand()
    if(np.exp(action_init - action_fin) > metropolis):
      #accept
      naccept=naccept+1
    else:
      p=backup_p

  if( ite >= 100): # discard first 100 samples for thermalization
    p_hist.append(p)

  # data output
  print(p, naccept/ite)

##########################
### plot the posterior ###
##########################
plt.xlabel("$p$")
plt.xlim(0.4,0.65)
plt.hist(p_hist, bins=30, density=True)
plt.show()












