# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:51:47 2015

@author: cat
"""

# ODE Method for solving the system


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Constants

V = 37 #units: cubic microns; do cells change size during heat stress? probably
# protein abundance -> using amounts from Kulak et al 2014 becuase its a methods paper
total_Pab = 95624. # molecules/cell
total30_HSP104 = 52655. # molecules/cell, should be roughly double that 15-60
                       # minutes of mild heat shock according to Kawai 1999

def deriv(z, t):
#    k1 = 10. # base Deactivation rate
#    k2 = 1. # Reactivation rate
#    k3 = .1 # Protein synthesis
#    k4 = .01 # Protein degradation
    
#    # Small number parameters
    k1 = .01 # Deactivation
    k2 = .002 # Reactivation
    k3 = .003 # Protein synthesis
    k4 = .001 # Protein degradation
    
    dPab = -1*k1*z[0]*deltaT + k2*z[1]*z[2]
    diPab = k1*z[0]*deltaT - k2*z[1]*z[2]
    dC = k3*z[1] - k4*z[2]    
    
    return np.array([dPab, diPab, dC])

deltaT = 1                  
time1 = np.arange(0, 100, .1)
#zinit = np.array([.85*total_Pab, .15*total_Pab, total30_HSP104])
zinit = np.array([85, 15, 20])
z1 = odeint(deriv, zinit, time1)

deltaT = 2
time2 = np.arange(100, 400, .1)
z2 = odeint(deriv, z1[-1], time2)

deltaT = 1
time3 = np.arange(400, 600, .1)
z3 = odeint(deriv, z2[-1], time3)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))






# Plots
names = ['Pab1', 'iPab1', 'C']
colors = ['royalblue', 'firebrick', 'gold']
for i in xrange(3):
    plt.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
#plt.plot(times, np.ones(len(times))*total_Pab, c='k')
plt.plot(times, np.ones(len(times))*100)
plt.xlabel('t')
plt.ylabel('Species Count')
plt.legend(loc=0)
#plt.xlim(99, 101)
plt.show()


    
#def deriv2(z, t):
#    k1 = .2
#    k2 = .1
#
#    dA = -1*k1*z[0] - k2*z[0]
#    dB = k1*z[0]
#    dC = k2*z[0]
#    return np.array([dA, dB, dC])
    
#time = np.arange(0, 50, 1)
#zinit = np.array([100., 0., 0.])
#z = odeint(deriv2, zinit, time)
#
#names = ['A', 'B', 'C']
#colors = ['royalblue', 'firebrick', 'gold']
#for i in xrange(3):
#    plt.plot(time, z[:, i], label=names[i], c=colors[i], linewidth=2)
#plt.xlabel('t')
#plt.ylabel('Species Count')
#plt.legend(loc=0)
#plt.ylim(-2, 100)
#plt.show()