# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:34:40 2015

@author: Triandafillou
"""

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
    x, y = cliff.simulate(200, "dummy_file")

    print np.mean(y), np.mean(y)/np.var(y)
    plt.subplot(121)
    plt.step(x, y)
    plt.title("Lemmings vs. time")
    plt.xlabel("Time (s)")
    plt.ylabel("Number of Lemmings")
    plt.subplot(122)
    plt.title("Lemming count distribution")
    plt.hist([i[0] for i in y], normed=True)
    plt.xlabel("Lemming population")
    plt.show()
    
if __name__ == "__main__":
    main()