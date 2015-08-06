#!/usr/bin/python
"""
Basic implementation of the Gillespie algorithm.

Designed for unimolecular creation/destruction reactions
"""

import numpy as np
import matplotlib.pyplot as plt

class Species(object):
    """A chemical species. Has an inherent property count."""
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
        
    def prop(self):
        alpha = 1
        for i in self.ip:
            alpha *= i.count
        return alpha*self.rate
        
    def perform(self, ip, op, rate):
        #insert specific reaction details here
        pass
        
class ConstInduction(Reaction):
    """Specific Reaction instance with a rate independent of species concentration"""
    def prop(self):
        return self.rate
        
    def perform(self):
        #quick fix now, better fix later:
        self.op[0].produce()
        
class UniDeg(Reaction):
    """Specific Reaction instance of degredation with a rate dependent on the
    concentration of the species being degraded.
    """
    def perform(self):
        #quick fix now, better fix later
        self.ip[0].destroy()
        
        
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
        
    def update_species(self, sp_list, rxn_list, alpha_vec, rn):
        pass
        
    def simulate(self, t_max, file_path):
        """Still very exploratory/could be crashy. Make sure to double check
        rate/times so that the simulation doesn't get out of control.
        Current implementation is specifically for two unimolecular creation/
        destruction reactions.
        """
        #f = open(file_path, "w")
        #f.write(str(self.species) + "\n")
        t = 0
        n_steps = 0
        #temporary data storage:
        data = []
        time = []
        print self.species
        while t < t_max:
            #again, generalization
            data.append(self.species[0].count)
            time.append(t)
            r1, r2 = np.random.uniform(), np.random.uniform()
            alpha_vec = [i.prop() for i in self.reactions]
            alpha = np.sum(alpha_vec)
            try:
                tau = 1./alpha*np.log(1./r1)
            except ZeroDivisionError:
                print "Propensity equal to zero at step = %d, time = %d" % (n_steps, t)
            t = t + tau
            n_steps += 1
            #update species here - could generalize much much more:
            if r2 < 1./alpha:
                self.reactions[0].perform()
            else:
                self.reactions[1].perform()
            if n_steps % 10 == 0:
                print "Number of steps = %d" % n_steps
            
        print "Simulation finished after %d steps" % n_steps     
        #f.close()
        return time, data
        
def main():
    L = Species(0)
    arrival = ConstInduction("Induction", None, [L], 1)
    jump = UniDeg("Degredation", [L], None, 0.1)
    cliff = Network([L], [arrival, jump])
    x, y = cliff.simulate(100, "dummy_file")
    #plots
    plt.step(x, y)
    plt.xlabel("Time (s)")
    plt.ylabel("Number of Lemmings")
    plt.show()  
    
if __name__ == "__main__":
    main()