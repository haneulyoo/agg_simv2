# -*- coding: utf-8 -*-
"""
Modeling only inactivation and reactivation - no dimerization step
"""

import numpy as np
import matplotlib.pyplot as plt
from gsim_A import Reaction, TDReaction, TConstInduction, Species, Network, UniDeg

class MonomerReactivation(Reaction):
    """A reaction which converts a single inactivated species to two active.
    
    iA + C -> A + C
    In this particular instance disaggregation and reactivation are coupled.
    !!!For now -> the order of Species in the input/output list should be
    [inactivated, preserved] ([iA/A, C] in the above reaaction.
    """
    def perform(self):
        self.ip[0].destroy()
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
            return self.rate*self.ip[0].count*self.temp*np.sqrt(self.op[0].count)
    
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
    k2 = 0.01 #chaperone degradation rate
    k4 = 0.01 #heat induced chaperone production rate
    k5 = 0.01 #heat inactivation rate of assembler
    T1 = 1
    T2 = 0
    T3 = 20
    # Species
    A = Species("A", 100)
    iA = Species("iA", 0)
    C = Species("C", 10)
    # Temperature independent reactions
    disagg = MonomerReactivation("Disagg", [iA, C], [A, C], k1)
    deg = UniDeg("C Degredation", [C], [None], k2)
    # Simulation and Temperature-Dependent Reaction
    hip = HeatInducedProduction("Heat Induced Production", [iA], [C], k4, T3)
    inactivation = HeatInducedInactivation("Inactivation", [A], [iA], k5, T3)
    sp_list = [A, iA, C]
    rxn_list= [inactivation, disagg, hip, deg]
    system = Network(sp_list, rxn_list)
    x = system.simulate(0, 70, "None")
    #hip = HeatInducedProduction("Heat Induced Production", [iA], [C], k4, T3)
    #inactivation = HeatInducedInactivation("Inactivation", [A], [iA], k5, T3)
    #y = system.simulate(x[-1,0], 400, "None")
    #hip = HeatInducedProduction("Heat Induced Production", [iA], [C], k4, T1)
    #inactivation = HeatInducedInactivation("Inactivation", [A], [iA], k5, T1)
    #z = system.simulate(y[-1,0], 70, "None")
    
    #final = np.vstack((x,y))
    #final = np.vstack((x, y, z))
    #x2 = [i[-1,0] for i in [x, y, z]]
    #y2 = [T1, T2, T1]
    colors = ['b','r','c']
    fig, axis = plt.subplots(figsize=(7,7)) 
    for i in range(1,len(sp_list)+1):
        axis.step(x[:,0], x[:,i], label=sp_list[i-1].name, c=colors[i-1])
    plt.legend(loc=0)
    plt.xlim(0,x[-1,0])
    plt.xlabel("Time")
    plt.ylabel("Molecular species count")
    #axis2 = axis.twinx()
    #axis2.step(x2,y2, c='darkgray')
    #plt.ylim(T1,T2)
    #plt.fill_between([x2[0], x2[1]], [220,220], alpha=0.5, color='darkgray')
    #plt.ylabel("Temperature (T)")
    plt.savefig("reactivation_p1(wsqrt).pdf")
    plt.show()
    
if __name__ == "__main__":
    main()