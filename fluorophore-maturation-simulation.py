"""
Motivation: As in previous version, but using an 'effective' Pab1 concentration rather than an actual 

Species: Pab, iPab, C, mRNAC, Pab_mRNAC
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import gridspec


# Parameters from von der Haar 2008
total_HSP104 = 28591 # unshocked conditions
HSP104_deg = 0.011363 # min^-1
#total_Pab1 = 100115 # unshocked conditions
total_Pab1 = 13000
base_HSP104_mRNA = 4.7
total_cellular_mRNA = 36000 #&
estimated_HSP104_hs = 250 #molecules/cell



def deriv(z, t):
    Ea = 80.                             # activation energy (arb. units)
    k1 = .1*np.exp(Ea*(1-(303./T)))    # mRNA production (fxn of temperature)
    k2 = 1                           # Protein synthesis min^-1 
    k3 = . 5                           # Fluorophore maturation
    k4 = .1                           # Protein degradation/dilution
    k5 = 0.005                          # mRNA degradation
    
    # Initial conditions
    mRNAC = z[0]
    C = z[1]
    Cx = z[2]
    

    # ODEs
    dmRNAC = k1 - k5*mRNAC
    dC = k2*mRNAC - k4*C - k3*C
    dCx = k3*C - k4*Cx
        
    return np.array([dmRNAC, dC, dCx])


T = 303 
time1 = np.arange(0, 10.0, .01)
zinit = np.array([1, 0, 0])
z1 = odeint(deriv, zinit, time1)
Tlist = [T, T]

T = 317
time2 = np.arange(10.0, 30.0, .01)
z2 = odeint(deriv, z1[-1], time2)
Tlist.append(T)

T = 303
time3 = np.arange(30.0, 250.0, .01)
z3 = odeint(deriv, z2[-1], time3)
Tlist.append(T)

times = np.concatenate((time1, time2, time3))
final = np.vstack((z1, z2, z3))

#print('Total C mRNA after heat shock: ' + str(z2[-1, 3] + z2[-1, 4]))

# Plots
names = ['$Chaperone mRNA$', '$Chaperone$', '$Visible Chaperone$']
colors = ['k', 'grey', 'firebrick']

f = plt.figure(figsize=(8, 5))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax2 = plt.subplot(gs[0])
ax2.step([time1[0], time1[-1], time2[-1], time3[-1]], Tlist, c='k')
ax2.set_ylabel('$\Delta$ T (K)')
ax2.set_ylim(290, 320)
plt.setp(ax2.get_xticklabels(), visible=False)
ax = plt.subplot(gs[1], sharex=ax2)
for i in range(len(names)):
    ax.plot(times, final[:, i], label=names[i], c=colors[i], linewidth=2)
#ax.plot(times, final[:, 0] + final[:, -1], c='b', linewidth=2, label='Total active Pab1')
#ax.plot(times, final[:, 3] + final[:, 4], label='total $mRNA_C$', color='papayawhip', linewidth=3)
ax.set_xlabel('time (min)')
ax.set_ylabel('Species Count')
#ax.set_yscale('log')
ax.legend(loc='upper right', fontsize=10)
plt.tight_layout()

print(final[-1, :])