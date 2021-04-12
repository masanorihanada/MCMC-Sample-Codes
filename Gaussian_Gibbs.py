import numpy as np
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default
# Numpy has Gaussian random number generator, so we don't use Box-Muller

NITER = 10000
A = [[1,0,1.0,1.0],
     [1.0,2.0,1.0],
     [1.0,1.0,2.0]]
A = np.array(A)
#####################################
### Set the initial configuration ###
#####################################
x = 0
y = 0
z = 0
#################
### Main part ###
#################

# for plot
data1=[]
data2=[]
data3=[]

for iter in range(NITER):
    #update x
    sigma = 1.0/np.sqrt(A[0][0])
    mu = -A[0][1]/A[0][0]*y-A[0][2]/A[0][0]*z
    r = np.random.randn()
    x = sigma*r+mu
    #update y
    sigma = 1.0/np.sqrt(A[1][1])
    mu = -A[1][0]/A[1][1]*x-A[1][2]/A[1][1]*z
    r = np.random.randn()
    y = sigma*r+mu
    #update z
    sigma = 1.0/np.sqrt(A[2][2])
    mu = -A[2][0]/A[2][2]*x-A[2][1]/A[2][2]*y
    r = np.random.randn()
    z = sigma*r+mu
    #output x,y,z
    print(x,y,z)
    data1.append(x)
    data2.append(y)
    data3.append(z)

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
###
ax3 = plt.subplot(1, 3, 3)
ax3.scatter(data2,data3,c="red", s=20, marker="+")
ax3.set_xlim([-8,8])
ax3.set_ylim([-8,8])
ax3.set_xlabel("y")
ax3.set_ylabel("z")
ax3.set_title("y-z")
ax3.grid(linestyle = "--")
plt.show()

