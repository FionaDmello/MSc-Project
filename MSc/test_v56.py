import numpy as np
import math
import matplotlib.pyplot as plt
import time

#First Run of Dark matter line-like distribution with an inhomogeneity, relevant conditions at the boundary applied, constant force (1D). Time Scale - order of million years(MY), Distance scale - Mega parsec.

#Defining constants
pi = 3.1415926
n=1501						#Total number of particles.
rbar = 34.61 * (10**24)				# in meters 																	
M = 397.2 * (10**45)
l = 0.03				#Mass carried by each mass-element.					 

#Defining time steps
T = 1000								#Total no. of time steps																	
dt = 0.01								#Each time step is 0.01 MY.
t = np.linspace(0,T,num = 1001,retstep = True) 	
tlen = len(t[0])							

#Defining array:
d = np.zeros(n)
Ndis0 = np.zeros((n-1))
x = np.zeros((n,tlen))
v = np.zeros((n,tlen))
X = np.zeros(n)
V = np.zeros(n)
NdisI = np.zeros((n-1))
I = np.zeros((n-1))
ind = np.zeros((n-1), dtype = 'int')
xlen = len(x)

#Initial conditions
x[(n-1)/2][0] = 0													
v[0][0] = 0	
d[(n-1)/2] = ((1.332) / ((1 + l)**250))			#in Mpc, distance between the central particle and its nearest neighbors.
                            												
#Creating 2D particle indices for plotting 
for i in xrange(n-1):
	if i == 0:
		ind[i] = 0
	else:
		ind[i] = ind[i-1] + 1																
									
#Setting other intial conditions:
for p in xrange((((n-1)/2)-1),-1,-1):
	if (p > 499):							#Setting up the intial positions of all particles with a mass density inhomogeneity at the center.
		d[p] = d[p+1] + (l * d[p+1])
		x[p][0] = x[p+1][0] - (d[p])	
	else:
		d[p] = d[p+1]
		x[p][0] = x[p+1][0] - (d[p])

for p in xrange((((n-1)/2)+1),n):
	if (p < 1001 ):
		d[p] = d[p-1] + (l * d[p-1])
		x[p][0] = x[p-1][0] + (d[p])
	else:
		d[p] = d[p-1]
		x[p][0] = x[p-1][0] + (d[p])

for p in xrange(0,n): 
	v[p][0] = v[0][0]	#Initial velocities being set to zero.
	if (p < (n-1)):		#Saving the initial value of distance between neighboring particles of the system in Ndis0.
		Ndis0[p] = x[p+1][0] - x[p][0]

np.savetxt("1Ndis0.txt",Ndis0)
np.savetxt("1d.txt",d)

#Defining and calculating the acceleration
def accel(M,x,i,K):
	accel_S = 0 																	
	for j in xrange(n):							 
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.00001)**0.5)												
		accel_S = accel_S + accel_I										  	#Sum of effective acceleration components
	return (2 * pi * ((2.1502*(10**-6))/ (rbar**2)) * accel_S)         # in Mpc/(MY)^2.

#Defining function for velocity calculations
def vel(M,dt,v,x,i,K):
	if K == 0:
		return v[i][0]
	else:
		return  v[i][K-1] + accel(M,x,i,K-1) * dt						# in Mpc/(MY).

#Defining function for position calculations
def pos(dt,x,i,K):
	return x[i][K] + vel(M,dt,v,x,i,K) * dt								# in Mpc.

#Calculating positions and velocities
if __name__ == "__main__":
	start = time.time()  
	for k in xrange(1,tlen):
		for i in xrange(500,1001):
			v[i][k] = vel(M,dt,v,x,i,k)
			x[i][k] = pos(dt,x,i,k-1)
		for j in xrange(0,500): 							
			x[j][k] = x[500][k] - ((500 - j)*d[j])		#Setting position of dummy particles wrt the nearest real particle at every time step.
			v[j][k] = v[500][k]							#Making the dummy particle follow the nearest real particle in terms of velocity.
		for j in xrange(1001,n):									
			x[j][k] = x[1000][k] + ((j - 1000)*d[j])			
			v[j][k] = v[1000][k]

	X0 = x.T[tlen-tlen]							#Information from the initial state, zeroth time.
	V0 = v.T[tlen-tlen]

	X = x.T[tlen-1]								#Information at the end of 10 MY.
	V = v.T[tlen-1]

	for i in xrange(500,1001):						#Calculating the distance of separation between neighbors at after 10 MY.
		NdisI[i] = X[i+1] - X[i]
		I[i] = (((((NdisI[i])**2)**0.5) - Ndis0[i])/Ndis0[i])
	
	np.savetxt("1X1.txt",X)					#Saving information about the system at 10 MY, to be used as initial conditions for the next run.
	np.savetxt("1V1.txt",V)
	end = time.time()
	print end - start
	
	#Plotting phase space plots at 0 and 10 MY. Also plotting the D vs index graph. 
	plt.figure()
	plt.plot(X0,V0,'r.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(x[i][0],v[i][0],'k.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()	
	plt.plot(X,V,'m.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(X[i],V[i],'r.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(ind[i],I[i],'b.')
		plt.xlabel('Index i(from 500-1001)')
		plt.ylabel('Stretch Parameter, D')
	plt.show()

	plt.figure()
	for i in xrange(500,1001):	
		plt.plot(X[i],ind[i],'g.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Index i(from 500-1001)')
	plt.show()

