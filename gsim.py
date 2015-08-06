#!/usr/bin/python
"""
Basic implementation of the Gillespie algorithm.
__author__
"""

import numpy as np

class Species(object):
    """A chemical species. Has an inherent count."""
    def __init__(self, count):
        self.count = count
    
    def produce(self):
        self.count += 1
        
    def destroy(self):
        self.count -= 1

class Reaction(object):
    """A Reaction used in Gillespie simulation.
    
    A general class which defines the required outputs for a reaction
    (propensity); inputs include name of the reaction (str), input (list of
    species objects), output (list of species objects), rate (float), and
    potentially a description of the reaction (str).
    """
    def __init__(self, name, ip, op, rate):
        self.name = name
        self.rate = rate
        self.ip = ip
        self.op = op
        
    def prop(self, ip, rate):
        alpha = 1
        for i in ip:
            alpha *= i.count
        return alpha*rate
        
class ConstInduction(Reaction):
    """Specific Reaction instance with a rate independent of species concentration"""
    def prop(self, rate):
        return rate
        
        
class Network(object):
    """Holds a reaction network which can be acted on during the Gillespie
    loop. Initializing a network requires a list of Species (with
    associated inital counts) and a list of possible Reactions (with associated
    rates, inputs, outputs, etc.)
    Make sure to import numpy as np
    """
    def __init__(self, sp_list, rxn_list):
        self.species = sp_list
        self.reactions = rxn_list
        
    def simulate(self, t_max, file_path):
        """Still very exploratory/could be crashy. Make sure to double check
        rate/times so that the simulation doesn't get out of control.
        """
        f.open(file_path, "w")
        f.write(str(self.species) + "\n")
        t = 0
        n_steps = 0
        while t < t_max:
            r1, r2 = np.random.uniform(), np.random.uniform()
            alpha = 0
            for i in self.reactions:
                alpha += i.prop()
            try:
                tau = 1./alpha
            except ZeroDivisionError:
                print "Propensity equal to zero at step: %d, time: %d" % (n_steps, t)
            
       print "Simulation finished after %d steps" % n_steps     
       f.close() 
        
#def main(sp_dict, )