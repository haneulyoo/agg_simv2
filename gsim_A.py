#!/usr/bin/python
"""
Basic implementation of the Gillespie algorithm.
"""

import numpy as np

class Species(object):
    """A chemical species. Has the inherent properties name and count."""
    def __init__(self, name, count):
        self.count = count
        self.name = name
    
    def produce(self):
        self.count += 1
        
    def destroy(self):
        self.count -= 1

class Reaction(object):
    """Defines the base class for a temperature-independent reaction.
    
    A general class which defines the required outputs for a reaction
    (propensity); inputs include name of the reaction (str), input (list of
    species objects), output (list of species objects), base rate (float), and
    potentially a description of the reaction (str).
    The base class has no perform function! You must explicitly define the
    reaction for each new reaction (for now?)
    Reaction.prop now takes temperature as an input (but does not use)
    """
    def __init__(self, name, ip, op, Ea, baserate):
        self.name = name
        self.ip = ip
        self.op = op
        self.Ea = Ea
        self.baserate = baserate
    
    def prop(self, T):
        """Return the propensity. Doesn't use T."""
        alpha = self.baserate
        for i in self.ip:
            alpha *= i.count
            return alpha
        
    def perform(self):
        #raise error here to ensure definition of the reaction?
        pass
    
class TDReaction(Reaction):
    """Defines the base class for a temperature dependent reaction.
    
    Currently implementation uses the Arrehenius equation to determine the rate
    """
    def prop(self, T):
        alpha = self.baserate*T
        for i in self.ip:
            alpha *= i.count
            return alpha
            
#    def prop(self, T):
#        alpha = self.baserate*np.exp(-self.Ea/1.381e-23*(1-(298/T))) #Kb in units of J/K; base temperature assumed to be 298 K
#        for i in self.ip:
#            alpha *= i.count
#            return alpha
        
class ConstInduction(Reaction):
    """Production from nothing.

    0 -> A
    """
    def prop(self, T):
        return self.baserate
        
    def perform(self):
        self.op[0].produce()
        
class TConstInduction(TDReaction):
    """Production from nothing at a temperature dependent rate."""
    #def prop(self, T):
    #    return self.baserate*np.exp(-self.Ea/1.381e-23*(1-(298/T)))
    
    def prop(self, T):
        return self.baserate*T
    
    def perform(self):
        self.op[0].produce()
        
class UniDeg(Reaction):
    """Specific Reaction instance of degredation.
    
    A -> 0
    The propensity is dependent on the concentration of the species being degraded.
    """
    def perform(self):
        self.ip[0].destroy()
        
        
class Network(object):
    """Holds a reaction network which can be acted on during the Gillespie
    loop.
    
    Initializing a network requires a list of Species (with
    associated inital counts) and a list of possible Reactions (with associated
    rates, inputs, outputs, etc.)
    Make sure to import numpy as np
    """
    def __init__(self, sp_list, rxn_list):
        self.species = sp_list
        self.reactions = rxn_list
        
    def determine_temperature(self, temp_fxn, t):
        """Given a step function of time and temperature and the current time
        return the current temperature.
        
        temp_fxn should be formatted as a list of tuples of the form
        (time, temperature) for all desired changes in temperature.
        Default T is 298 K.
        """
        if t < temp_fxn[0][0]:
            return 298
        elif t > temp_fxn[len(temp_fxn)-1][0]:
            return temp_fxn[len(temp_fxn)-1][1]
        else:
            for i in range(1,len(temp_fxn)):
                if t < temp_fxn[i][0]:
                    return temp_fxn[i-1][1]
        
        
    def simulate(self, t_start, t_end, temp_fxn, filename):
        """Implementation of the Gillespie SSA.
        
        No check on runtime. Make sure to double check rate/times so that
        simulation doesn't get out of control.
        temp_fxn should be formatted as a list of tuples of the form (time, temperature)
        for all desired changes in temperature; t_start < time < t_end.
        Simulation is assumed to start at 298 K. Start with (0,alternate_temp)
        to alter this behavior.
        """
        assert t_start < t_end
        t = t_start
        n_steps = 0
        time = [[t]]
        data = [[i.count for i in self.species]]
        temp = [[298]]
        while t <= t_end:
            T = self.determine_temperature(temp_fxn, t)
            r1, r2 = np.random.uniform(), np.random.uniform()
            alpha_vec = [i.prop(T) for i in self.reactions]
            alpha = np.sum(alpha_vec)
            try:
                tau = 1./alpha*np.log(1./r1)
            except ZeroDivisionError:
                print "Propensity equal to zero at step = %d, time = %d;\
                Simulation terminated." % (n_steps, t)
                break

            t += tau
            n_steps += 1
            #update species; there is definitely a better/faster way to do this
            #also find a way to ensure that a reaction happens at every time step
            #lead -> no two rows in the final array should be the same except time
            z = 0
            for i in range(len(alpha_vec)):
                z += alpha_vec[i]
                if float(z)/alpha == 1.0:
                    self.reactions[i].perform()
                    #print self.reactions[i].name
                    break
                elif r2 < float(z)/alpha:
                    self.reactions[i].perform()
                    #print self.reactions[i].name
                    break
            time.append([t]) 
            data.append([i.count for i in self.species])
            temp.append([T])
            if n_steps % 500 == 0:
                print "Step: %d" % n_steps
        
        print "Simulation finished after %d steps" % n_steps
        data = np.asarray(data)
        time = np.asarray(time)
        temp = np.asarray(temp)
        total = np.hstack((time, temp, data))
        #np.savetxt(filename, total, header = "Time " + str([i.name for i in self.species]))
        return total
