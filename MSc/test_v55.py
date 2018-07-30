import numpy as np
import math
import matplotlib.pyplot as plt
import time

#Second run of Dark matter line-like distributin with an inhomogeneity, relevant conditions at the boundary applied, constant force (1D). Time scale in order of million years, distance scale in Mpc.

#Defining constants
pi = 3.1415926
n = 1501
rbar = 34.61 * (10**24)																								
M = 397.2 * (10**45)

#Defining time steps
T = 1000																						
dt = 0.01 		
t = np.linspace(0,T,num = 1001,retstep = True) 
tlen = len(t[0])

#Defining arrays:
x = np.zeros((n,len(t[0])))
v = np.zeros((n,len(t[0])))
NdisI = np.zeros((n-1))
I = np.zeros((n-1))
ind = np.zeros((n-1), dtype = 'int')
xlen = len(x) 
																									
#Setting intial conditions:
Lx = np.loadtxt("1501X46.txt", unpack = True)	#Loading information from the end of the previous run as initial conditions for the current run.
Lv = np.loadtxt("1501V46.txt", unpack = True)
Ndis0 = np.loadtxt("1501Ndis0.txt", unpack = True)
d = np.loadtxt("1501d.txt", unpack = True)

for p in xrange(0,n):
	x[p][0] = Lx[p]
	v[p][0] = Lv[p]

ind[0] = 0
#Creating 2D particle indices for plotting 
for i in xrange(n-1):
	if i == 0:
		ind[i] = 0
	else:
		ind[i] = ind[i-1] + 1
	
#Defining and calculating the acceleration
def accel(M,x,i,K):
	accel_S = 0 																			
	for j in xrange(n):							 
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.00001)**0.5)													
		accel_S = accel_S + accel_I
	return (2 * pi * ((2.1502 * (10**-6))/ (rbar**2)) * accel_S)

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
		for i in xrange(500,1001):
			v[i][k] = vel(M,dt,v,x,i,k)
			x[i][k] = pos(dt,x,i,k-1)
		for j in xrange(0,500):									
			x[j][k] = x[500][k] - ((500 - j) * d[j])		
			v[j][k] = v[500][k]
		for j in xrange(1001,n):									
			x[j][k] = x[1000][k] + ((j - 1000) * d[j])			
			v[j][k] = v[1000][k]

	Xn = x.T[tlen - 1]
	Vn = v.T[tlen - 1]
	
	#Xf = np.vstack((X,Xn))
	#Vf = np.vstack((V,Vn))

	for i in xrange(500,1001):
		NdisI[i] = Xn[i+1] - Xn[i]
		I[i] = (((((NdisI[i])**2)**0.5) - Ndis0[i])/Ndis0[i])
	
	np.savetxt("1501X47.txt",Xn)
	np.savetxt("1501V47.txt",Vn)

	end = time.time()
	print end - start
	
	plt.figure()
	for i in xrange(500,1001):
		plt.plot(Xn[i],Vn[i],'b.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(ind[i],I[i],'r.')
		plt.xlabel('Index i(from 500-1001)')
		plt.ylabel('Stretch Parameter, D')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(Xn[i],ind[i],'g.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Index i(from 500-1001)')
	plt.show()
	"""
	plt.figure()
	plt.hist(X,200,facecolor = 'blue')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Number of particles')
	plt.show()
	"""
