#!/usr/bin/python
'''
Simulates reactivation of deactivated species in a temperature dependent fashion.
'''
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from gsim import Reaction, Species, Network
from inactivation_sim import Inactivation, Temp_Dimerization

class Reactivation(Reaction):
    """A reaction which converts a single inactivated species to two active.
    
    iAA + C -> A + A + C
    In this particular instance disaggregation and reactivation are coupled.
    !!!For now -> the order of Species in the input/output list should be
    [inactivated, preserved] ([iAA/A, C] in the above reaaction.
    """
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()
        self.op[0].produce()
        
class HeatInducedProduction(Reaction):
    """A reaction where the rate of production is temperature dependent.
    
    0 -> A
    """
    def __init__(self, name, ip, op, rate, T):
        self.name = name
        self.rate = rate
        self.ip = ip
        self.op = op
        self.temp = T
   
    def prop(self):
        return self.rate*self.temp
        
    def perform(self):
        self.op[0].produce()
        
def main():
    # Species
    A = Species("A", 100)
    AA = Species("AA", 0)
    iAA = Species("iAA", 0)
    C = Species("C", 0)
    # Temperature independent reactions
    inactivation = Inactivation("Inactivation", [AA], [iAA], 0.005)
    disagg = Reactivation("Disagg", [iAA, C], [A, C], 0.005)
    # Simulation and Temperature-Dependent Reaction
    # Temperature ramp = 300, 320, 300
    nuc = Temp_Dimerization("Nucleation", [A], [AA], 0.00001, 300)
    hip = HeatInducedProduction("Heat Induced Production", [None], [C], 0.001, 300)
    sp_list = [A, AA, iAA, C]
    system = Network(sp_list, [nuc, inactivation, disagg, hip])
    x = system.simulate(0, 50, "None")
    
    fig, axis = plt.subplots()    
    for i in range(1,len(sp_list)+1):
        axis.step(x[:,0], x[:,i], label=sp_list[i-1].name)
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main()