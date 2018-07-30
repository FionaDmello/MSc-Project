import numpy as np
import math
import matplotlib.pyplot as plt
import time

#First Run of Dark matter line-like distribution with an inhomogeneity, relevant conditions at the boundary applied, constant force (1D). Time Scale - order of million years(MY), Distance scale - Mega parsec.

#Defining constants
pi = 3.1415926
n = 211								#Total no.of particles.																	
rbar = 30.6695 * (10**21)				# in meters = 1Mpc																		
M = 0.226899154 * (10**45)				#Mass carried by each mass-element.					 

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
d[(n-1)/2] = (1 / ((1.12)**24))			#in Mpc, distance between the central particle and its nearest neighbors.
                            												
#Creating 2D particle indices for plotting 
for i in xrange(n-1):
	if i == 0:
		ind[i] = 0
	else:
		ind[i] = ind[i-1] + 1																
									
#Setting other intial conditions:
for p in xrange((((n-1)/2)-1),-1,-1):
	if (p > 79):							#Setting up the intial positions of all particles with a mass density inhomogeneity at the center.
		d[p] = d[p+1] + (0.12 * d[p+1])
		x[p][0] = x[p+1][0] - (d[p])	
	else:
		d[p] = d[p+1]
		x[p][0] = x[p+1][0] - (d[p])

for p in xrange((((n-1)/2)+1),n):
	if (p < 131 ):
		d[p] = d[p-1] + (0.12 * d[p-1])
		x[p][0] = x[p-1][0] + (d[p])
	else:
		d[p] = d[p-1]
		x[p][0] = x[p-1][0] + (d[p])

for p in xrange(0,n): 
	v[p][0] = v[0][0]	#Initial velocities being set to zero.
	if (p < (n-1)):		#Saving the initial value of distance between neighboring particles of the system in Ndis0.
		Ndis0[p] = x[p+1][0] - x[p][0]

np.savetxt("Ndis00.txt",Ndis0)
np.savetxt("d2.txt",d)

#Defining and calculating the acceleration
def accel(M,x,i,K):
	accel_S = 0 																	
	for j in xrange(n):							 
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.0001)**0.5)												
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
		for i in xrange(55,156):
			v[i][k] = vel(M,dt,v,x,i,k)
			x[i][k] = pos(dt,x,i,k-1)
		for j in xrange(0,55): 							
			x[j][k] = x[55][k] - ((55 - j)*d[j])		#Setting position of dummy particles wrt the nearest real particle at every time step.
			v[j][k] = v[55][k]							#Making the dummy particle follow the nearest real particle in terms of velocity.
		for j in xrange(156,n):									
			x[j][k] = x[155][k] + ((j - 155)*d[j])			
			v[j][k] = v[155][k]

	X0 = x.T[tlen-tlen]							#Information from the initial state, zeroth time.
	V0 = v.T[tlen-tlen]

	X = x.T[tlen-1]								#Information at the end of 10 MY.
	V = v.T[tlen-1]

	for i in xrange(55,156):						#Calculating the distance of separation between neighbors at after 10 MY.
		NdisI[i] = X[i+1] - X[i]
		I[i] = (((((NdisI[i])**2)**0.5) - Ndis0[i])/Ndis0[i])
	
	np.savetxt("X'''1.txt",X)					#Saving information about the system at 10 MY, to be used as initial conditions for the next run.
	np.savetxt("X'''.txt",X)
	np.savetxt("V'''1.txt",V)
	np.savetxt("V'''.txt",V)
	end = time.time()
	print end - start
	
	#Plotting phase space plots at 0 and 10 MY. Also plotting the D vs index graph. 
	plt.figure()
	plt.plot(X0,V0,'b.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(x[i][0],v[i][0],'g.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()	
	plt.plot(X,V,'r.')
	plt.xlabel('Position (Mpc)')
	plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(X[i],V[i],'r.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(ind[i],I[i],'k.')
		plt.xlabel('Index i(from 55-156)')
		plt.ylabel('Stretch Parameter, D')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(ind[i],X[i],'m.')
		plt.xlabel('Index i(from 55-156)')
		plt.ylabel('Position (Mpc)')
	plt.show()
