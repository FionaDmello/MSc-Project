import numpy as np
import math
import matplotlib.pyplot as plt
import time

#Homogeneous line-like dark matter distribution, constant force (1D). No dummy particles, just the study system. Time scale in order of million years, distance scale in Mpc.

#Defining constants
pi = 3.1415926
n = 101																		
rbar =	30.6695 * (10**21)																							
M = 0.226899154 * (10**45)									 

#Defining time steps
T = 1000																									
dt = 0.01
t = np.linspace(0,T,num = 1001,retstep = True) 	
tlen = len(t[0])							

#Defining array:
d = np.zeros(n)
Ndisz = np.zeros((n-1))
x = np.zeros((n,tlen))
v = np.zeros((n,tlen))
X = np.zeros(n)
V = np.zeros(n)
Ndis = np.zeros((n-1))
I = np.zeros((n-1))
ind = np.zeros((n-1), dtype = 'int')
xlen = len(x)

#Initial conditions
x[(n-1)/2][0] = 0													
v[0][0] = 0	
d[(n-1)/2] = (1 / ((1.0)**49))	#Mpc
                            												
#Creating 2D particle indices for plotting 
for i in xrange(n-1):
	if i == 0:
		ind[i] = 0
	else:
		ind[i] = ind[i-1] + 1																
									
#Setting other intial conditions:
for p in xrange((((n-1)/2)-1),-1,-1):
	if (p > 54):
		d[p] = d[p+1] + (0 * d[p+1]) 
		x[p][0] = x[p+1][0] - (d[p])	
	else:
		d[p] = d[p+1] 
		x[p][0] = x[p+1][0] - (d[p])

for p in xrange((((n-1)/2)+1),n):
	if (p < 156 ):
		d[p] = d[p-1] + (0 * d[p-1])
		x[p][0] = x[p-1][0] + (d[p])
	else:
		d[p] = d[p-1] 
		x[p][0] = x[p-1][0] + (d[p])

for p in xrange(0,n):
	v[p][0] = v[0][0]

#Defining and calculating the acceleration
def accel(M,x,i,K):
	accel_S = 0 																	
	for j in xrange(n):							 
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.0001)**0.5)													
		accel_S = accel_S + accel_I
	return (2 * pi * ((0.647686264*(10**-6))/ (rbar**2)) * accel_S)         

#Defining function for velocity calculations
def vel(M,dt,v,x,i,K):
	if K == 0:
		return v[i][0]
	else:
		return  v[i][K-1] + accel(M,x,i,K-1) * dt

#Defining function for position calculations
def pos(dt,x,i,K):
	return x[i][K] + vel(M,dt,v,x,i,K) * dt

#Calculating positions and velocities
if __name__ == "__main__":
	start = time.time()  
	for k in xrange(1,tlen):
		for i in xrange(0,n):
			v[i][k] = vel(M,dt,v,x,i,k)
			x[i][k] = pos(dt,x,i,k-1)

	X0 = x.T[tlen-tlen]
	V0 = v.T[tlen-tlen]

	X = x.T[tlen-1]	
	V = v.T[tlen-1]

	end = time.time()
	print end - start
	
	plt.figure()
	plt.plot(X0,V0,'g.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()
	
	plt.figure()	
	plt.plot(X,V,'g.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

