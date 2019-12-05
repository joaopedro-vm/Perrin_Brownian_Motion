#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# *************************************************************************
#
# Physics Institute, University of Brasília, Brasil
# João Pedro Valeriano Miranda
# 
# A simple script to calculate de mean square deviations in the
# brownian motion trajectories measured by Jean Perrin, plot the MSD by the
# interval time, and estimate the Avogadro's number by means of the
# Einstein's relation.
#
# *************************************************************************

# Necessary additional packages.
import numpy as np # Array operations.
import matplotlib.pyplot as plt # Plotting tools.
from scipy.optimize import curve_fit # Fitting the MSD vs interval time plot.

R = 8.314 # Universal gas constant.
T = 290 # System temperature.
eta = 0.0011 # Water viscosity.
r = 0.52e-6 # Particless' radius.

# List to store the trajectories.
d = []

# Import trajectories.
for i in range(1,4):
    d.append(np.genfromtxt(r"traj%i.txt"%i, delimiter=","))

# Turn trajectories list into an array.
d = np.array([d[i] for i in range(3)], dtype=object)

# Function to calculate the 2D MSD over a 'dat' trajectory using the minimal
# time interval.
def msd2D(dat):
    
    x = np.array(dat[:,0])
    y = np.array(dat[:,1])
    
    x = np.sum([(x[i+1] - x[i])**2 for i in range(len(x)-1)])/(len(x)-1)
    y = np.sum([(y[i+1] - y[i])**2 for i in range(len(y)-1)])/(len(y)-1)
    
    return x+y

###########################################################################

# Array to store MSD data.
msd = np.zeros((0,2))

N_dt = 5

# Fill the 'msd' array with data for the first 'N_dt' time intervals.
for dat in d:
    for i in range(1, N_dt):
        for j in range(0, i):
        
            k = dat[j::i]
        
            msd = np.append(msd, [[msd2D(k[:,1:]), k[1,0]-k[0,0]]], axis=0)

# Estimate and print the Avogadro's number usign the Eistein's relation.
N = 2*msd[:,1]/msd[:,0]*R*T/(3*np.pi*eta*r)
print("N_Avogadro = %.1e +/- %.0e" % (np.mean(N), np.std(N)/np.sqrt(len(N))))

###########################################################################

# Reinitialize 'msd' array.
msd = np.zeros((0,2))

# Calculates the MSD for all possible time intervals.
for dat in d:
    for i in range(len(dat)-1):
        for j in range(i+1, len(dat)):
        
            msd = np.append(msd, [[(dat[j,1]-dat[i,1])**2+(dat[j,2]-dat[i,2])**2, dat[j,0]-dat[i,0]]], axis=0)

# Function take the average and standard deviation of the MSD for each time
# interval.            
def t_mean(m):
    
    m = m[m[:,1].argsort()]
    
    p = 0
    l = []
    for i in m[:,1]:
        
        if (i != p):
            p = i
            l.append(p)
    
    mean = np.zeros((0,3))
    
    while (len(mean) < len(l)):
        
        e = []
        
        j = np.shape(mean)[0]
        
        for i in range(len(m)):
            
            if (m[i,1] == l[j]):
                
                e.append(m[i,0])
                
            else:
                break
        
        mean = np.append(mean, [[l[j],np.mean(e),np.std(e)/np.sqrt(len(e))]], axis=0)
        
        m = m[i:]
        
    return mean

m = t_mean(msd)

def f(x, a):
    return a*x

###########################################################################

#Plots

popt, pcov = curve_fit(f, m[:18,0], m[:18,1])
plt.plot(m[:,0], popt[0]*m[:,0], color="red", zorder=1)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.errorbar(m[:,0], m[:,1], yerr=m[:,2], color="black", ecolor="grey", zorder=0)
plt.xlabel("Tempo (s)")
plt.ylabel(r"Desvio Quadrático Médio (m$^{2}$)")

plt.savefig(r"msd-t-perrin.pdf")
