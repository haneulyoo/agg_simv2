#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 09:19:22 2017
@author: Triandafillou

Simulation of a simple positive autoregulatory loop with the species X
activating its own transcription
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

plt.style.use('ggplot')

# parameters
beta_max = 3
k_deg = 1
K = 1
n = 2

initial_X = 2

def deriv(z, t):
    X = z[0]
    
    dX = (beta_max * X ** n) / (K ** n + X ** n) - k_deg*X

    return(np.array([dX]))
    

t = np.arange(0, 10, 0.01)
zinit = np.array([initial_X])

z_f = odeint(deriv, zinit, t)
z_f2 = odeint(deriv, np.array([0.5]), t)
z_f3 = odeint(deriv, np.array([0.7739]), t)

#plt.figure(figsize=(4,2.5))
#plt.plot(t, z_f, 'r--', t, z_f2, 'b--', t, z_f3, 'g--')
#plt.xlabel("Time")
#plt.ylabel("[X]")
#plt.tight_layout
#plt.show



# plot initial vs. final
Xi = np.arange(0.37, 0.4, 0.0001)
Xf = []

for i in Xi:
    z_f = odeint(deriv, np.array([i]), t)
    Xf.append(z_f[-1])

plt.figure(figsize=(4,2.5))
plt.plot(Xi, Xf, '.')
plt.xlabel("Inital [X]")
plt.ylabel("Final [X]")
plt.tight_layout
plt.show

