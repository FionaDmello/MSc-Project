import numpy as np
import math
import matplotlib.pyplot as plt
import time

#Second run of Dark matter line-like distributin with an inhomogeneity, relevant conditions at the boundary applied, constant force (1D). Time scale in order of million years, distance scale in Mpc.

#Defining constants
pi = 3.1415926
n = 211
rbar = 30.6695 * (10**21)																								
M = 0.226899154 * (10**45)

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
Lx = np.loadtxt("X'''9.txt", unpack = True)	#Loading information from the end of the previous run as initial conditions for the current run.
Lv = np.loadtxt("V'''9.txt", unpack = True)
X = np.loadtxt("X'''9.txt", unpack = True)
V = np.loadtxt("V'''9.txt", unpack = True)
Ndis0 = np.loadtxt("Ndis00.txt", unpack = True)
d = np.loadtxt("d2.txt", unpack = True)

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
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.0001)**0.5)													
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
		for i in xrange(55,156):
			v[i][k] = vel(M,dt,v,x,i,k)
			x[i][k] = pos(dt,x,i,k-1)
		for j in xrange(0,55):									
			x[j][k] = x[55][k] - ((55 - j) * d[j])		
			v[j][k] = v[55][k]
		for j in xrange(156,n):									
			x[j][k] = x[155][k] + ((j - 155) * d[j])			
			v[j][k] = v[155][k]

	Xn = x.T[tlen - 1]
	Vn = v.T[tlen - 1]
	
	Xf = np.vstack((X,Xn))
	Vf = np.vstack((V,Vn))

	for i in xrange(55,155):
		NdisI[i] = Xn[i+1] - Xn[i]
		I[i] = (((((NdisI[i])**2)**0.5) - Ndis0[i])/Ndis0[i])
	
	np.savetxt("X'''10.txt",Xn)
	np.savetxt("V'''10.txt",Vn)
	np.savetxt("X'''.txt",Xf)
	np.savetxt("V'''.txt",Vf)

	end = time.time()
	print end - start
	
	plt.figure()
	for i in xrange(55,156):
		plt.plot(Xn[i],Vn[i],'b.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(55,155):	
		plt.plot(ind[i],I[i],'k.')
		plt.xlabel('Index i(from 55-156)')
		plt.ylabel('Stretch Parameter, D')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(ind[i],Xn[i],'g.')
		plt.xlabel('Index i(from 55-156)')
		plt.ylabel('Position (Mpc)')
	plt.show()

	plt.figure()
	plt.hist(X,150,facecolor = 'blue')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Number of particles')
	plt.show()
