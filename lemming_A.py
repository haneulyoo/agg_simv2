#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from gsim_A import Species, TConstInduction, ConstInduction, UniDeg, Network

def main():
    """Generate and simulate a model of lemmings approaching and jumping off
    a cliff.
    """
    # Input species and reactions
    L = Species("Lemming", 0)
    arrival = TConstInduction("Induction", None, [L], -1.0e-23, 1.0)
    jump = UniDeg("Degredation", [L], None, -1.0e-23, 0.1)
    temp_fxn = [(500,350)]
    cliff = Network([L], [arrival, jump])
    x = cliff.simulate(0, 1000, temp_fxn, "dummy_file2.dat")
    y = x[:,2]    
    y1, y2 = x[:500,2], x[500:,2]
    print "Mean1: %f, Mean1/Var1: %f" % (np.mean(y1), np.mean(y1)/np.var(y1))
    print "Mean2: %f, Mean2/Var2: %f" % (np.mean(y2), np.mean(y2)/np.var(y2))
    plt.subplot(212)
    plt.step(x[:,0], y)
    plt.title("Lemmings vs. time")
    plt.xlabel("Time (s)")
    plt.ylabel("Number of Lemmings")
    plt.legend(loc=0)
    plt.subplot(211)
    plt.step(x[:,0],x[:,1], label='Time')
    plt.ylabel("Temperature")
    plt.show()
    
if __name__ == "__main__":
    main()