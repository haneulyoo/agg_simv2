#!/usr/bin/python
'''
Simulates reactivation of deactivated species in a temperature dependent fashion.
'''
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from gsim import Reaction, Species, Network, UniDeg
from inactivation_sim import Temp_Dimerization

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
    Also dependent on the concentration of inactivated molecules in the
    simulation. Include this parameter in the ip list.
    """
    def __init__(self, name, ip, op, rate, T):
        self.name = name
        self.rate = rate
        self.ip = ip
        self.op = op
        self.temp = T
   
    def prop(self):
        return self.rate*self.temp*self.ip[0].count
        
    def perform(self):
        self.op[0].produce()
        
class HeatInducedInactivation(Reaction):
    """A reaction in which a species is inactivated at a rate proportional to
    the concentration of previously inactivated species and the temperature.
    
    A -> iA
    Generic rate is in units s^-1; the propensity of the reaction is dependent
    linearly on the number of activated molecules, and also on the square root
    of the number of inactivated molecules, which roughly approximates the
    surface area of the growing aggregated (inactivated) bulk species.
    """    
    def __init__(self, name, ip, op, rate, T):
        self.name = name
        self.rate = rate
        self.ip = ip
        self.op = op
        self.temp = T
    
    def prop(self):
        if self.op[0].count == 0:
            return self.rate*self.ip[0].count
        else:
            return self.rate*self.ip[0].count*np.sqrt(self.op[0].count)*self.temp
    
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()
        
        
def main():
    """Simulate temperature dependent aggregation and disaggregation.
    
    Clearly need to explicitly write temperature into the simulate method in
    the future.
    """
    # Rates
    k1 = 0.01 #disaggregation rate
    k2 = 0.1 #chaperone degradation rate
    k3 = 0.00001 #dimerization rate
    k4 = 0.01 #heat induced chaperone production rate
    k5 = 0.01 #heat inactivation rate of assembler
    T1 = 200
    T2 = 220
    # Species
    A = Species("A", 100)
    AA = Species("AA", 0)
    iAA = Species("iAA", 0)
    C = Species("C", 10)
    # Temperature independent reactions
    disagg = Reactivation("Disagg", [iAA, C], [A, C], k1)
    deg = UniDeg("C Degredation", [C], [None], k2)
    # Simulation and Temperature-Dependent Reaction
    # Temperature ramp = 300, 320, 300
    nuc = Temp_Dimerization("Nucleation", [A], [AA], k3, T1)
    hip = HeatInducedProduction("Heat Induced Production", [iAA], [C], k4, T1)
    inactivation = HeatInducedInactivation("Inactivation", [AA], [iAA], k5, T1)
    sp_list = [A, AA, iAA, C]
    rxn_list= [nuc, inactivation, disagg, hip, deg]
    system = Network(sp_list, rxn_list)
    x = system.simulate(0, 20, "None")
    #nuc = Temp_Dimerization("Nucleation", [A], [AA], k3, 30)
    #hip = HeatInducedProduction("Heat Induced Production", [iAA], [C], k4, T2)
    #inactivation = HeatInducedInactivation("Inactivation", [AA], [iAA], k5, T2)
    #y = system.simulate(x[-1,0], 50, "None")
    #nuc = Temp_Dimerization("Nucleation", [A], [AA], k3, 10)
    #hip = HeatInducedProduction("Heat Induced Production", [iAA], [C], k4, T1)
    #inactivation = HeatInducedInactivation("Inactivation", [AA], [iAA], k5, T1)
    #z = system.simulate(y[-1,0], 70, "None")
    
    #final = np.vstack((x, y, z))
    #x2 = [i[-1,0] for i in [x, y, z]]
    #y2 = [T1, T2, T1]
   
    fig, axis = plt.subplots(figsize=(10,10)) 
    for i in range(1,len(sp_list)+1):
        axis.step(x[:,0], x[:,i], label=sp_list[i-1].name)
    plt.legend(loc=0)
    plt.xlim(0,x[-1,0])
    plt.xlabel("Time")
    plt.ylabel("Molecular species count")
    #axis2 = axis.twinx()
    #axis2.step(x2,y2, c='darkgray')
    #plt.ylim(T1,T2)
    #plt.fill_between([x2[0], x2[1]], [220,220], alpha=0.5, color='darkgray')
    #plt.ylabel("Temperature (T)")
    #plt.savefig("reactivation_plot.pdf")
    plt.show()
    
if __name__ == "__main__":
    main()