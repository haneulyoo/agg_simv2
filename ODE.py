# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:51:47 2015

@author: cat
"""

# ODE Method for solving the system


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def deriv(z, t):
    k1 = .1 # Deactivation
    k2 = .1 # Reactivation
    k3 = 1 # Protein synthesis
    k4 = .1 # Protein degradation
    T = 1
    
    dPab = -k1*z[0]*T + k2*z[1]*z[2]
    diPab = -k2*z[1]*z[2] + k1*z[0]*T
    dC = k3*z[1] - k4*z[2]    
    
    return np.array([dPab, diPab, dC])
                    
time = np.arange(0, 100, 1)
yinit = np.array([100., 0., 50])
z = odeint(deriv, yinit, time)

names = ['Pab1', 'iPab1', 'C']
colors = ['royalblue', 'firebrick', 'gold']
for i in xrange(3):
    plt.plot(time, z[:, i], label=names[i], c=colors[i], linewidth=2)
plt.xlabel('t')
plt.ylabel('Species Count')
plt.legend(loc=0)
plt.show()

    
def deriv2(z, t):
    k1 = .2
    k2 = .1

    dA = -1*k1*z[0] - k2*z[0]
    dB = k1*z[0]
    dC = k2*z[0]
    return np.array([dA, dB, dC])
    
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