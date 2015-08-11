#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from gsim import ConstInduction, Network, Species, UniDeg

def main():
    """Generate and simulate a model of lemmings approaching and jumping off
    a cliff.
    """
    L = Species("Lemming", 0)
    arrival = ConstInduction("Induction", None, [L], 1.0)
    jump = UniDeg("Degredation", [L], None, 0.1)
    cliff = Network([L], [arrival, jump])
    x = cliff.simulate(0, 200, "dummy_file2.dat")
    y = x[:,1]
    print "Mean: %f, Mean/Variance: %f \n (Poisson mean/var = 1)" % (np.mean(y), np.mean(y)/np.var(y))
    plt.subplot(121)
    plt.step(x[:,0], y)
    plt.title("Lemmings vs. time")
    plt.xlabel("Time (s)")
    plt.ylabel("Number of Lemmings")
    plt.subplot(122)
    plt.title("Count distribution")
    plt.hist(y, normed=True)
    plt.xlabel("Lemming population")
    plt.show()
    
if __name__ == "__main__":
    main()