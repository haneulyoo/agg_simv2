# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:48:31 2015

@author: cat
"""

import numpy as np
import matplotlib.pyplot as plt
from gsim_A import Reaction, Species, Network, UniDeg, TDReaction

class Pab1Deactivation(TDReaction):
    """Reaction representing Pab aggregation. Temperature dependent.
    
    Pab1 -> iPab1
    """
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()

class Pab1Reactivation(Reaction):
    """A reaction which converts a single inactivated species to two active.
    
    iPab1 + C -> Pab1 + C
    In this particular instance disaggregation and reactivation are coupled.
    !!!For now -> the order of Species in the input/output list should be
    [inactivated, preserved] ([iAA/AA, C] in the above reaaction.
    """
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()
        
class CProduction(Reaction):
    """Chaperone translation reaction.
    
    0 -> C
    
    Chaperone translation is inversely proportionl to Pab1 concentration 
    (currently implemented as 1 / [Pab1]. Pab1 should be the input for the
    reaction.
    """
    def prop(self, T):
        try:            
            x = self.baserate*(1. / self.ip[0].count)
        except ZeroDivisionError:
            x = self.baserate
        return x
        
    def perform(self):
        self.op[0].produce()
        
class CDegradation(Reaction):
    """Chaperone degradation reaction.
    
    C -> 0
    
    Rate is a general protein degradation rate.
    """
    def perform(self):
        self.ip[0].destroy()
        
        
        
def main():
    # Rates
    k1 = .01 # Deactivation
    k2 = .01 # Reactivation
    k3 = .01 # Protein synthesis
    k4 = .01 # Protein degradation    
    temp_fxn = [(0,298)]
    # Species
    Pab1 = Species("Pab1", 100)
    C = Species("Chaperone", 0)
    iPab1 = Species("iPab1", 0)
    # Reactions
    aggregation = Pab1Deactivation("agg", [Pab1], [iPab1], None, k1)
    disaggregation = Pab1Reactivation("disagg", [iPab1, C], [Pab1, C], None, k2)
    Ctranslation = CProduction("C Translation", [Pab1], [C], None, k3)
    Cdeg = CDegradation("C Degradation", [C], [], None, k4)
    # System
    species_list = [Pab1, C, iPab1]
    rxn_list = [aggregation, disaggregation, Ctranslation, Cdeg]
    system = Network(species_list, rxn_list)
    # Simulation
    x = system.simulate(0, 1000, temp_fxn, None)
    # Plots
    colors = ['royalblue', 'gold', 'firebrick']
    fig, (axis2, axis) = plt.subplots(2, 1, figsize=(7,7))
    for i in range(2,len(species_list)+2):
        axis.plot(x[:100,0], x[:100,i], label=species_list[i-2].name, c=colors[i-2], linewidth=3)
    axis.legend(loc=0)
    axis.set_xlim(0,1)#x[-1,0])
    axis.set_xlabel("Time")
    axis.set_ylabel("Molecular species count")
    axis2.step(x[:, 0], x[:, 1], linewidth=3)
    axis2.set_xlabel("Time (s)")
    axis2.set_ylabel("Temperature")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()