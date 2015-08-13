#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from gsim_A import Species, TConstInduction, ConstInduction, UniDeg, Network

def main():
    """Generate and simulate a model of lemmings approaching and jumping off
    a cliff.
    """
    # Input species and reactions
    L = Species("Lemming", 0)
    arrival = TConstInduction("Induction", None, [L], -1.0e-23, 1.0)
    jump = UniDeg("Degredation", [L], None, -1.0e-23, 0.1)
    temp_fxn = [(50,350)]
    cliff = Network([L], [arrival, jump])
    x = cliff.simulate(0, 100, temp_fxn, "dummy_file2.dat")
    y = x[:,2]    
    y1, y2 = x[:50,2], x[50:,2]
    print "Mean1: %f, Mean1/Var1: %f" % (np.mean(y1), np.mean(y1)/np.var(y1))
    print "Mean2: %f, Mean2/Var2: %f" % (np.mean(y2), np.mean(y2)/np.var(y2))
    
    #plotting function -> plots temperature over timecourse w/ shared x axis    
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=(1,3))
    axis1 = plt.subplot(gs[1])
    axis1.step(x[:,0], y)
    plt.xlabel("Time (s)")
    axis1.set_ylabel("Lemmings")
    axis0 = plt.subplot(gs[0], sharex=axis1)
    axis0.step(x[:,0],x[:,1])
    axis0.set_ylabel("Temperature")
    plt.setp(axis0.get_xticklabels(), visible=False)
    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":
    main()