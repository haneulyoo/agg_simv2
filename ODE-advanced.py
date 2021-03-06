#!/bin/python
"""
Motivation: previous model (ODE-basic) recapitulates relative levels of Pab1,
insoluble Pab1, and chaperone on long time scales, but fails to capture the
expected delay between heat shock induced protein aggregation and chaperone
translation. This iteration of the model includes explicit mRNA-protein inter-
actions in an attempt to address this shortcoming.


Species: Pab, iPab, C, mRNAC, Pab_mRNAC
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import gridspec
#import time


# Parameters from von der Haar 2008
total_HSP104 = 28591 # unshocked conditdions
HSP104_deg = 0.011363 # min^-1
total_Pab1 = 100115 # unshocked conditions
#Pab1_deg = 0.000891 # min^-1; basically 0 for our purposes
base_HSP104_mRNA = 4.7

def deriv(z, t):
    Ea = 40. # assembly reaction activation energy (arb. units for now)
    Ea2 = 100. # mRNA production activation energy (arb. units for now)
    k1 = 1.*np.exp(Ea*(1-(303./T))) # Deactivation (M^-1*min^-1)
    k2 = .00005 # Reactivation 
    k3 = 105. # Protein synthesis M^-1*min^-1
    k4 = HSP104_deg # Protein degradation
    k5 = 5.*np.exp(Ea2*(1-(303./T))) # mRNA production rate (min^-1)
    k6 = .001 # Pab-mRNA binding rate M^-2*min^-1
    km6 = .01 # Pab-mRNA unbinding rate M^-1*min^-1
    
    dPab = -k1*z[0] + k2*z[1]*z[2] - k6*z[0]*z[3] + km6*z[4]
    diPab = k1*z[0] - k2*z[1]*z[2]
    dC = k3*z[3] - k4*z[2]
    dmRNAC = k5 - k6*z[0]*z[3] + km6*z[4]
    dPab_mRNAC = k6*z[0]*z[3] - km6*z[4]
    
    return np.array([dPab, diPab, dC, dmRNAC, dPab_mRNAC])


T = 317                 
time1 = np.arange(0, 5, .01)
zinit = np.array([total_Pab1, 0, total_HSP104, 5, 0])
z1 = odeint(deriv, zinit, time1)
Tlist = [T, T]

T = 317
time2 = np.arange(5, 30, .01)
z2 = odeint(deriv, z1[-1], time2)
Tlist.append(T)

T = 303
time3 = np.arange(30, 45, .01)
z3 = odeint(deriv, z2[-1], time3)
Tlist.append(T)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))


# Plots
names = ['$Pab1$', '$iPab1$', '$C$', '$mRNA_C$', '$Pab1:mRNA_C$']
colors = ['royalblue', 'firebrick', 'gold', 'darkgreen', 'k']

f = plt.figure(figsize=(10, 5))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
ax2 = plt.subplot(gs[0])
ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
ax2.set_ylabel('$\Delta$ T (K)')
ax2.set_ylim(290, 320)
plt.setp(ax2.get_xticklabels(), visible=False)
ax = plt.subplot(gs[1], sharex=ax2)
for i in xrange(5):
    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
ax.set_xlabel('time (min?)')
ax.set_ylabel('Species Concentration')
ax.legend(loc=0)
#ax.set_yscale('log')
plt.tight_layout()