#!/bin/python
"""
Motivation: previous model (ODE-basic) recapitulates relative levels of Pab1,
insoluble Pab1, and chaperone on long time scales, but fails to capture the
expected delay between heat shock induced protein aggregation and chaperone
translation. This iteration of the model includes explicit mRNA-protein inter-
actions in an attempt to address this shortcoming.


Species: Pab, iPab, C, mRNAC, Pab_mRNAC, mRNAB, Pab_mRNAB
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import gridspec


# Parameters from von der Haar 2008
total_HSP104 = 28591 # unshocked conditdions
HSP104_deg = 0.011363 # min^-1
total_Pab1 = 100115 # unshocked conditions
#Pab1_deg = 0.000891 # min^-1; basically 0 for our purposes
base_HSP104_mRNA = 4.7
total_cellular_mRNA = 36000 #*
estimated_HSP104_hs = 250 #molecules/cell


def deriv(z, t):
    Ea = 80. # assembly reaction activation energy (arb. units for now)
    Ea2 = 125. # mRNA production activation energy (arb. units)
    k1 = .01*np.exp(Ea*(1-(303./T))) # Deactivation (M^-1*min^-1)
    k2 = .000001 # Reactivation 
    k3 = 1000 # Protein synthesis M^-1*min^-1
    #k4 = HSP104_deg # Protein degradation
    k4 = 0.005 #Protein degradation
    k5 = .1*np.exp(Ea2*(1-(303./T))) # mRNA production rate (min^-1)
    k6 = .003 # Pab-mRNA binding rate M^-2*min^-1
    km6 = .03 # Pab-mRNA unbinding rate M^-1*min^-1
    k7 = .1 # mRNA decay rate    
    
    Pab = z[0]
    iPab = z[1]
    C = z[2]
    mRNAC = z[3]
    Pab_mRNAC = z[4]
    mRNAB = z[5]
    Pab_mRNAB = z[6]
    
    dPab = -k1*Pab + k2*iPab*C - k6*Pab*mRNAC + km6*Pab_mRNAC - k6*Pab*mRNAB + km6*Pab_mRNAB
    diPab = k1*Pab - k2*iPab*C
    dC = k3*mRNAC - k4*C
    dmRNAC = k5 - k6*Pab*mRNAC + km6*Pab_mRNAC - k7*mRNAC
    dPab_mRNAC = k6*Pab*mRNAC - km6*Pab_mRNAC
    dmRNAB = km6*Pab_mRNAB - k6*Pab*mRNAB
    dPab_mRNAB = k6*Pab*mRNAB - km6*Pab_mRNAB
    
    return np.array([dPab, diPab, dC, dmRNAC, dPab_mRNAC, dmRNAB, dPab_mRNAB])


T = 303 
time1 = np.arange(0, 10.0, .01)
zinit = np.array([total_Pab1, 0, total_HSP104, 5, 0, total_cellular_mRNA, 0])
z1 = odeint(deriv, zinit, time1)
Tlist = [T, T]

T = 317
time2 = np.arange(10.0, 20.0, .01)
z2 = odeint(deriv, z1[-1], time2)
Tlist.append(T)

T = 303
time3 = np.arange(20.0, 80.0, .01)
z3 = odeint(deriv, z2[-1], time3)
Tlist.append(T)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))

print 'Total C mRNA after heat shock: ' + str(z2[-1, 3] + z2[-1, 4])

# Plots
names = ['$Pab1$', '$iPab1$', '$C$', '$free mRNA_C$', '$Pab1:mRNA_C$', '$free mRNA_B$', '$Pab:mRNA_B$']
colors = ['royalblue', 'firebrick', 'gold', 'darkgreen', 'k', 'indigo', 'darkorange']

f = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax2 = plt.subplot(gs[0])
ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
ax2.set_ylabel('$\Delta$ T (K)')
ax2.set_ylim(290, 320)
plt.setp(ax2.get_xticklabels(), visible=False)
ax = plt.subplot(gs[1], sharex=ax2)
for i in xrange(len(names)):
    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
#ax.plot(times, final[:, 3] + final[:, 4], label='total $mRNA_C$', color='papayawhip', linewidth=3)
ax.set_xlabel('time (min)')
ax.set_ylabel('Species Count')
ax.set_yscale('log')
ax.legend(loc='upper right', fontsize=10)
plt.tight_layout()

print final[-1, :]