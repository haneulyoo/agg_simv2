#!/bin/python
"""
Created on Sun Oct 18 12:03:28 2015

@author: cat
"""

# ODE Method for solving the system


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import gridspec

# Parameters from von der Haar 2008
total_HSP104 = 28591 # unshocked conditions
HSP104_deg = 0.011363 # min^-1
total_Pab1 = 100115 # unshocked conditions
#Pab1_deg = 0.000891 # min^-1; basically 0 for our purposes

def deriv(z, t):
    Ea = 40
    k1 = 10.*np.exp(Ea*(1-(298./T))) # Deactivation
    k2 = .001 # Reactivation
    k3 = 1. # Protein synthesis
    k4 = HSP104_deg # Protein degradation
    
    dPab = -k1*z[0] + k2*z[1]*z[2]
    diPab = -k2*z[1]*z[2] + k1*z[0]
    dC = k3*(z[1]/20.) - k4*z[2]    
    
    return np.array([dPab, diPab, dC])

T = 298                 
time1 = np.arange(0, 20, .1)
zinit = np.array([total_Pab1, 0, total_HSP104])
z1 = odeint(deriv, zinit, time1)
Tlist = [T, T]

T = 298
time2 = np.arange(20, 30, .1)
z2 = odeint(deriv, z1[-1], time2)
Tlist.append(T)

T = 298
time3 = np.arange(30, 50, .1)
z3 = odeint(deriv, z2[-1], time3)
Tlist.append(T)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))






# Plots
names = ['Pab1', 'iPab1', 'C']
colors = ['royalblue', 'firebrick', 'gold']
#f, ax = plt.subplots(figsize=(12,5))
#for i in xrange(3):
#    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
#ax.set_xlabel('time')
#ax.set_ylabel('Species Concentration')
#plt.legend(loc=0)
#f2, ax2 = plt.subplots(figsize=(12,1))
#ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
#ax2.set_xlabel('time')
#ax2.set_ylabel('$\Delta$ T (K)')
#ax2.set_ylim(-1, 30)


f = plt.figure(figsize=(10, 5))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
ax2 = plt.subplot(gs[0])
ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
ax2.set_ylabel('$\Delta$ T (K)')
ax2.set_ylim(290, 320)
plt.setp(ax2.get_xticklabels(), visible=False)
ax = plt.subplot(gs[1], sharex=ax2)
for i in xrange(3):
    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
ax.set_xlabel('time (min?)')
ax.set_ylabel('Species Concentration')
ax.legend(loc=0)
plt.tight_layout()
#plt.savefig('nohs.pdf')