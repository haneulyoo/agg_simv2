#!/bin/python
"""
Motivation: What are some possible explanations for why we see less protein
product after pH-altered heat shock?

Possible reasons (not mutually exclusive nor exhaustive):
1. Less transcript produced
2. Less translational efficiency (translated mRNA/total mRNA, relative measure)
3. Less ribosomes (translated mRNA/total mRNA, absolute measure)
4. More mRNA degradation (same amount of transcript produced but steady-state
   value is reduced.)
5. More protein degradation (same amount of protein produced but steady-state
   value is reduced.)


Species: Pab, iPab, C, mRNAC, Pab_mRNAC, B, mRNAB, Pab_mRNAB
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
    Ea = 50. # assembly reaction activation energy (arb. units for now)
    Ea2 = 150. # mRNA production activation energy (arb. units for now)
    k1 = 10.*np.exp(Ea*(1-(303./T))) # Deactivation (M^-1*min^-1)
    k2 = .0001 # Reactivation 
    k3 = .1 # Protein synthesis
    k4 = HSP104_deg # Protein degradation
    k5 = 50.*np.exp(Ea2*(1-(303./T))) # mRNA production rate (min^-1)
    k6 = .001 # Pab-mRNA binding rate
    km6 = 1. # Pab-mRNA unbinding rate
    k7 = .01 # mRNA decay rate    
    
    Pab = z[0]
    iPab = z[1]
    C = z[2]
    mRNAC = z[3]
    Pab_mRNAC = z[4]
    mRNAB = z[5]
    Pab_mRNAB = z[6]
    B = z[7]
    
    dPab = -k1*Pab + k2*iPab*C - k6*Pab*mRNAC + km6*Pab_mRNAC - k6*Pab*mRNAB + km6*Pab_mRNAB
    diPab = k1*Pab - k2*iPab*C
    dC = k3*mRNAC - k4*C
    dmRNAC = k5 - k6*Pab*mRNAC + km6*Pab_mRNAC - k7*mRNAC
    dPab_mRNAC = k6*Pab*mRNAC - km6*Pab_mRNAC
    dmRNAB = km6*Pab_mRNAB - k6*Pab*mRNAB
    dPab_mRNAB = k6*Pab*mRNAB - km6*Pab_mRNAB
    dB = k3*Pab_mRNAB - k4*B
    
    return np.array([dPab, diPab, dC, dmRNAC, dPab_mRNAC, dmRNAB, dPab_mRNAB, dB])


T = 317                 
time1 = np.arange(0, 1.0, .001)
zinit = np.array([total_Pab1, 0, total_HSP104, 5, 0, 100000, 0, 40000])
z1 = odeint(deriv, zinit, time1)
Tlist = [T, T]

T = 317
time2 = np.arange(1.0, 4.0, .001)
z2 = odeint(deriv, z1[-1], time2)
Tlist.append(T)

T = 303
time3 = np.arange(4.0, 7.0, .001)
z3 = odeint(deriv, z2[-1], time3)
Tlist.append(T)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))

print(final[4000])

# Plots
names = ['$Pab1$', '$iPab1$', '$C$', '$mRNA_C$', '$Pab1:mRNA_C$', '$mRNA_B$', '$Pab:mRNA_B$', '$B$']
colors = ['royalblue', 'firebrick', 'gold', 'darkgreen', 'k', 'peru', 'darkorange', 'indigo']

f = plt.figure(figsize=(9, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax2 = plt.subplot(gs[0])
ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
ax2.set_ylabel('$\Delta$ T (K)')
ax2.set_ylim(290, 320)
plt.setp(ax2.get_xticklabels(), visible=False)
ax = plt.subplot(gs[1], sharex=ax2)
for i in range(8):
    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
ax.set_xlabel('time (min?)')
ax.set_ylabel('Species Count')
ax.legend(loc=0)
plt.tight_layout()