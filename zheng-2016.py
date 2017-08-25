"""
Motivation: As in previous version, but using an 'effective' Pab1 concentration rather than an actual 

Species: Pab, iPab, C, mRNAC, Pab_mRNAC
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Parameters
k1 = 166.8 # per min per au, client:hsp70 off rate
k3 = k1
k2 = 2.783 # per min, HSP70:HSF1 off rate
k4 = 0.0464 # per min HSP70:UP off rate
k5 = 4.64e-7 # degradation rate of HSP70:UP
beta = 1.778 # Transcription activation rate
Kd = 0.0022 # Dissociation constant of HSF1-DNA interaction
kdil = 2.78e-3 # Dilution rate of YFP; assuming 50% growth rate
n = 3 # Hill coefficient


# Initial values
HSP_0 = 1 # Free HSP70
HSF1_0 = 0 # Free HSF1
HSP_HSF1_0 = 1/500.
HSP_UP_0 = 0
YFP_0 = 3 # 



def deriv(z, t):

    
    HSP = z[0]
    HSF1 = z[1]
    HSP_HSF1 = z[2]
    HSP_UP = z[3]
    UP = z[4]
    YFP = z[5]
        
    dHSPdt = k2*HSP_HSF1 - k1*HSP*HSF1 + (k4 + k5)*HSP_UP - k3*HSP*UP + \
    beta*(HSF1**n / (Kd**n + HSF1**n))
    dHSF1dt = k2*HSP_HSF1 - k1*HSP*HSF1
    dHSP_HSF1dt = -k2*HSP_HSF1 + k1*HSP*HSF1
    dHSP_UPdt = -(k4 + k5)*HSP_UP + k3*HSP*UP
    dUPdt = k4*HSP_UP - k3*HSP*UP
    dYFPdt = beta*(HSF1**n / (Kd**n + HSF1**n)) - kdil*YFP
    
    
    return(np.array([dHSPdt, dHSF1dt, dHSP_HSF1dt, dHSP_UPdt, dUPdt, dYFPdt]))




# Inital test: T = 39C
T = [37, 39, 40, 42, 45]
z = []

for t in T:
    UP_0 = 0.0024*np.exp(0.215*t) # generates initial value of UP; empirical
    time = np.arange(0, 200, 0.1)
    zinit = np.array([HSP_0, HSF1_0, HSP_HSF1_0, HSP_UP_0, UP_0, YFP_0])
    z.append(odeint(deriv, zinit, time))


## Plots
names = ['37', '39', '40', '42', '45']
colors = ['royalblue', 'firebrick', 'gold', 'darkgreen', 'k']

f = plt.figure(figsize=(10,5))
ax = f.add_subplot(111)
for i in range(len(z)):
    ax.plot(time, z[i][:,5], color = colors[i], label = names[i])
ax.legend()
plt.title("Reporter level vs. time")