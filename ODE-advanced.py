# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 10:03:45 2015

@author: cat
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import gridspec

# Parameters from von der Haar 2008
total_HSP104 = 28591 # unshocked conditions
HSP104_deg = 0.011363 # min^-1
total_Pab1 = 100115 # unshocked conditions
#Pab1_deg = 0.000891 # min^-1; basically 0 for our purposes
base_HSP104_mRNA = 4.7