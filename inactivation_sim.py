# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:26:02 2015

@author: Triandafillou
"""

import numpy as np
import matplotlib.pyplot as plt
import gsim
#from gsim import ConstInduction, Network, Reaction, Species, UniDeg

class Inactivation(gsim.Reaction):
    """A reaction in which a species is inactivated at a rate proportional to
    the concentration of previously inactivated species.
    
    A -> iA
    """
    def perform(self, ip, op, rate):
        self.ip[0].destroy()
        self.op[0].destroy()

        
class Temp_Dimerization(gsim.Reaction):
    """A reaction in which a species dimerizes at a rate proportional to
    temperature.

    A + A -> AA
    Input should be a list of one species (in the above example A) and the
    output should be a list of a single dimerized species)
    """
    def prop(self):
        return self.count*(self.count-1)*self.rate
    
    def perform(self, ip, op, rate):
        self.ip[0].destroy()
        self.ip[0].destroy()
        self.op[0].produce()