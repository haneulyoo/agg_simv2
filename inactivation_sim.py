#!/usr/bin/python
'''
Simulates temperature-dependent aggregation where a nucleation step
(approximated by dimerization of the species in a temperature dependent fashion)
and dimers are subsequently inactivated at a rate proportional to the square
root of the number of inactivated species present. 
'''
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import gsim
#from gsim import ConstInduction, Network, Reaction, Species, UniDeg

class Inactivation(gsim.Reaction):
    """A reaction in which a species is inactivated at a rate proportional to
    the concentration of previously inactivated species.
    
    A -> iA
    Generic rate is in units s^-1; the propensity of the reaction is dependent
    linearly on the number of activated molecules, and also on the square root
    of the number of inactivated molecules, which roughly approximates the
    surface area of the growing aggregated (inactivated) bulk species.
    """
    #this could be recoded to change an internal state rather than introducing
    #a new species?
    def prop(self):
        if self.op[0].count == 0:
            return self.rate*self.ip[0].count
        else:
            return self.rate * self.ip[0].count * np.sqrt(self.op[0].count)
    
    def perform(self):
        self.ip[0].destroy()
        self.op[0].produce()

        
class Temp_Dimerization(gsim.Reaction):
    """A reaction in which a species dimerizes at a rate proportional to
    temperature.

    A + A -> AA
    Input should be a list of one species (in the above example A) and the
    output should be a list of a single dimerized species).
    """
    
    def __init__(self, name, ip, op, rate, T):
        self.name = name
        self.rate = rate
        self.ip = ip
        self.op = op
        self.temp = T
        
    def prop(self):
        return self.ip[0].count*(self.ip[0].count-1)*self.rate*self.temp
    
    def perform(self):
        self.ip[0].destroy()
        self.ip[0].destroy()
        self.op[0].produce()
        
        
def main():
    # Species in the system
    A = gsim.Species("A", 100)
    AA = gsim.Species("AA", 0)
    iAA = gsim.Species("iAA", 0)
    # Possible reactions
    inactivation = Inactivation("Inactivation", [AA], [iAA], 0.005)
    nuc = Temp_Dimerization("Nucleation", [A], [AA], 0.00001, 300)
    system = gsim.Network([A, AA, iAA], [nuc, inactivation])
    x = system.simulate(0, 50, "None")
    
    nuc = Temp_Dimerization("Nucleation", [A], [AA], 0.00001, 320)
    y = system.simulate(x[-1,0], 100, "None")
    
    final = np.vstack((x, y))
    print final, final.shape
    #temp_course = [320, 300]
    #t_add = 50
    #for T in range(len(temp_course)):
    #    t_total = 50
    #    nuc = Temp_Dimerization("Nucleation", [A], [AA], 0.0001, temp_course[T])
    #    system = gsim.Network([A, AA, iAA], [nuc, inactivation])
    #    x2 = system.simulate(t_total, t_total + t_add, "None")
    #    np.hstack(x1, x2)
    #    t_total += t_add
    
    #print x
    
    fig, axis = plt.subplots()    
    for i in range(1,4):
        axis.step(final[:,0], final[:,i])
    plt.show()

if __name__ == "__main__":
    main()