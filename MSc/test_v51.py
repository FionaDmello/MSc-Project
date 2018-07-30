import numpy as np
import math
import matplotlib.pyplot as plt
import time

#Second and multiple run of Dark matter line with uniform distance discretization, boundary conditions applied, 1/r form of force. The change in mass density distribution (form of inhomogeneity) is gradual. Time Scale MY. Canis Major Dwarf 25000LY.

#Defining constants
pi = 3.1415926
n = 211
a = 2.351679 * (10**20)																								
M = 1.0740 * (10**42)

#Defining time steps
T = 1000																										#Run-time.
dt = 0.001 		#1BY
t = np.linspace(0,T,num = 1001,retstep = True) 		#creates a 2D array. t[0] element of which is the array containing 1000 points, the next element showing what the step size between elements t[0].
tlen = len(t[0])

#Defining arrays
x = np.zeros((n,len(t[0])))
v = np.zeros((n,len(t[0])))
NdisI = np.zeros((n-1))
I = np.zeros((n-1))
ind = np.zeros((n-1), dtype = 'int')
xlen = len(x) 
																					#No. of time intervals created.					
#Setting intial conditions
Lx = np.loadtxt("X'''10.txt", unpack = True)		#Relativistic speeds at 0.08. No flip.
Lv = np.loadtxt("V'''10.txt", unpack = True)
X = np.loadtxt("X'''10.txt", unpack = True)
V = np.loadtxt("V'''10.txt", unpack = True)
Ndisz = np.loadtxt("Ndis02.txt", unpack = True)
b = np.loadtxt("b3.txt", unpack = True)

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
	accel_S = 0 																			#Sum of the effective acceleration components.
	for j in xrange(n):							 
		accel_I = M * ((x[j][K] - x[i][K])/(((x[j][K] - x[i][K])**2)+0.00001)**0.5)													
		accel_S = accel_S + accel_I
	return (2 * pi * ((0.647686264*(10**-6))/ (a**2)) * accel_S)

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
		for j in xrange(54,-1,-1):									
			x[j][k] = x[55][k] - ((55 - j)*b[j])		
			v[j][k] = v[55][k]
		for j in xrange(156,n):									
			x[j][k] = x[155][k] + ((j - 155)*b[j])			
			v[j][k] = v[155][k]


	Xt = x.T[tlen-501]
	Vt = v.T[tlen-501]

	Xn = x.T[tlen - 1]
	Vn = v.T[tlen - 1]
	
	Xf = np.vstack((X,Xn))
	Vf = np.vstack((V,Vn))

	for i in xrange(55,156):
		NdisI[i] = Xn[i+1] - Xn[i]
		I[i] =((((NdisI[i] - Ndisz[i])**2)**0.5)/Ndisz[i])
	
	np.savetxt("X'''11.txt",Xn)
	np.savetxt("V'''11.txt",Vn)
	np.savetxt("X'''.txt",Xf)
	np.savetxt("V'''.txt",Vf)

	end = time.time()
	print end - start
	
	plt.figure()
	for i in xrange(55,156):
		plt.plot(Xt[i],Vt[i],'g.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()
	
	plt.figure()
	for i in xrange(55,156):
		plt.plot(Xn[i],Vn[i],'b.')
		plt.xlabel('Position (Mpc)')
		plt.ylabel('Velocity (Mpc/MY)')
	plt.show()

	plt.figure()
	for i in xrange(55,156):	
		plt.plot(ind[i],I[i],'k.')
		plt.xlabel('Index i(from 55-156)')
		plt.ylabel('Stretch Parameter, D')
	plt.show()
