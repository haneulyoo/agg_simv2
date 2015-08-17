# -*- coding: utf-8 -*-
"""
Modeling only inactivation and reactivation - no dimerization step
"""

import numpy as np
import matplotlib.pyplot as plt
from gsim_A import Reaction, Species, Network, UniDeg

class MonomerReactivation(Reaction):
    """A reaction which converts a single inactivated species to two active.
    
    iAA + C -> A + A + C
    In this particular instance disaggregation and reactivation are coupled.
    !!!For now -> the order of Species in the input/output list should be
    [inactivated, preserved] ([iA/A, C] in the above reaaction.
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
    def prop(self, T):
        return self.baserate*T*self.ip[0].count
        
    def perform(self):
        self.op[0].produce()
        
class HeatInducedInactivation(Reaction):
    """A reaction in which a species is inactivated at a rate proportional to
    the concentration of previously inactivated species and the temperature.
    
    AA -> iAA
    Generic rate is in units s^-1; the propensity of the reaction is dependent
    linearly on the number of activated molecules, and also on the square root
    of the number of inactivated molecules, which roughly approximates the
    surface area of the growing aggregated (inactivated) bulk species.
    Temperature scaling currently at 1/12.
    """
    def prop(self, T):
        if self.op[0].count == 0:
            return self.baserate*self.ip[0].count
        else:
            return self.baserate*self.ip[0].count*T*0.09*np.sqrt(self.op[0].count)
    
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()
        
class Temp_Dimerization(Reaction):
    """A reaction in which a species dimerizes at a rate proportional to
    temperature.

    A + A -> AA
    Input should be a list of one species (in the above example A) and the
    output should be a list of a single dimerized species).
    """        
    def prop(self, T):
        return self.ip[0].count*(self.ip[0].count-1)*0.09*self.baserate*T
    
    def perform(self):
        self.ip[0].destroy()
        self.ip[0].destroy()
        self.op[0].produce()        
        
        
        
def main():
    """Simulate temperature dependent aggregation and disaggregation."""
    # Rates and Temperatures
    k1 = 0.01 #disaggregation rate
    k2 = 0.5 #chaperone degradation rate
    k3 = 0.00001 #dimerization rate
    k4 = 0.01 #heat induced chaperone production rate
    k5 = 0.001 #heat inactivation rate of assembler
    temp_fxn = [(1,298)]
    # Species
    A = Species("A", 100)
    AA = Species("AA", 0)
    iAA = Species("iAA", 0)
    C = Species("C", 10)
    # Temperature independent reactions
    disagg = MonomerReactivation("Disagg", [iAA, C], [A, C], None, k1)
    deg = UniDeg("C Degredation", [C], [None], None, k2)
    # Simulation and Temperature-Dependent Reaction
    nuc = Temp_Dimerization("Dimerization", [A], [AA], None, k3)
    hip = HeatInducedProduction("Heat Induced Production", [iAA], [C], None, k4)
    inactivation = HeatInducedInactivation("Inactivation", [AA], [iAA], None, k5)
    sp_list = [A, AA, iAA, C]
    rxn_list= [inactivation, disagg, hip, deg, nuc]
    system = Network(sp_list, rxn_list)
    x = system.simulate(0, 100, temp_fxn, "None")
    
    #x2 = [i[-1,0] for i in [x, y, z]]
    #y2 = [T1, T2, T1]
    colors = ['b','g','r','c']
    fig, axis = plt.subplots(figsize=(7,7)) 
    for i in range(2,len(sp_list)+2):
        axis.step(x[:,0], x[:,i], label=sp_list[i-2].name, c=colors[i-2])
    plt.legend(loc=0)
    plt.xlim(0,x[-1,0])
    plt.xlabel("Time")
    plt.ylabel("Molecular species count")
    #axis2 = axis.twinx()
    #axis2.step(x2,y2, c='darkgray')
    #plt.ylim(T1,T2)
    #plt.fill_between([x2[0], x2[1]], [220,220], alpha=0.5, color='darkgray')
    #plt.ylabel("Temperature (T)")
    #plt.savefig("reactivation_p1(wsqrt).pdf")
    plt.show()
    
if __name__ == "__main__":
    main()