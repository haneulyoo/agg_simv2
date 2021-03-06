#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from gsim_A import Species, TConstInduction, ConstInduction, UniDeg, Network

def main():
    """Generate and simulate a model of lemmings approaching and jumping off
    a cliff.
    
    Now with adjustable temperature!
    """
    # Input species and reactions
    L = Species("Lemming", 0)
    arrival = TConstInduction("Induction", None, [L], -1e-22, 0.1)
    jump = UniDeg("Degredation", [L], None, -1e-22, 0.1)
    temp_fxn = [(100,350), (150,298)]
    #temp_fxn = [(1,298)]
    cliff = Network([L], [arrival, jump])
    # Simulate and parse data    
    x = cliff.simulate(0, 200, temp_fxn, "dummy_file2.dat")
    y = x[:,2]    
    y1 = [x[i,2] for i in range(len(x)) if x[i,0] > 50 and x[i,0] < 100]
    y2 = [x[i,2] for i in range(len(x)) if x[i,0] > 100 and x[i,0] < 150]
    y3 = [x[i,2] for i in range(len(x)) if x[i,0] > 150]
    
    #print len(y1), len(y2), len(y3), len(x)
    #print "Mean1: %f, Mean1/Var1: %f" % (np.mean(y1), np.mean(y1)/np.var(y1))
    #print "Mean2: %f, Mean2/Var2: %f" % (np.mean(y2), np.mean(y2)/np.var(y2))
    #print "Mean3: %f, Mean3/Var3: %f" % (np.mean(y3), np.mean(y3)/np.var(y3))
    #plotting function -> plots temperature over timecourse w/ shared x axis    
    
    # Generate plots
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=(1,3))
    axis1 = plt.subplot(gs[1])
    axis1.step(x[:,0], y)
    plt.xlabel("Time (s)")
    axis1.set_ylabel("Lemmings")
    axis0 = plt.subplot(gs[0], sharex=axis1)
    axis0.step(x[:,0],x[:,1], c='r')
    axis0.set_ylabel("Temperature")
    plt.setp(axis0.get_xticklabels(), visible=False)
    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":
    main()