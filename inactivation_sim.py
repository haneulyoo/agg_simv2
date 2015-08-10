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
        return self.rate * self.ip[0].count * np.sqrt(self.op[0].count)
    
    def perform(self, ip, op):
        self.ip[0].destroy()
        self.op[0].create()

        
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
        return self.count*(self.count-1)*self.rate*T
    
    def perform(self, ip, op, rate):
        self.ip[0].destroy()
        self.ip[0].destroy()
        self.op[0].produce()
        
        
def main():
    # Species in the system
    A = gsim.Species("A", 100)
    AA = gsim.Species("AA", 0)
    iAA = gsim.Species("iAA", 0)
    # Possible reactions
    nuc = Temp_Dimerization("Nucleation", [A], [AA], 0.1)
    inactivation = Inactivation("Inactivation", [AA], [iAA], 1.0)
    # System Loop
    system = gsim.Network([A, AA, iAA], [nuc, inactivation])
    x, y = system.simulate()