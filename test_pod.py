#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 20:47:53 2018

@author: marjamaa

This is an invented toy example on how this POD works. 

Circuit has an input current and output voltage, four resistor values 
out of seven are adjustable. So the system formed from the  circuit by
nodal analysis is 1-in 1-out and parametrized by four 
resistor values (conductances).

Circuit has six nodes -> five DoFs.

After POD analysis we could reduce this to less than five by discarding the 
least important 'modes'. 

See circuit.png for the topology of the circuit.
"""

import numpy as np
import matplotlib.pyplot as plt
from pod import make_reducer

G3,G4,G5 = 1,1,1


def Zmu(params):
    """ Parametrized part of the nodal analysis admittance matrix """
    G1, G2, G6, G7 = params

    Z = np.zeros(shape=(5,5))
                  
    Z[0,0] = G1+G3
    Z[1,1] = G1+G2+G4 
    Z[2,2] = G2+G5
    Z[3,3] = G5+G7
    Z[4,4] = G4+G6+G7
    
    Z[0,1] = -G1
    Z[1,2] = -G2
    Z[1,0] = -G1
    Z[2,1] = -G2
    Z[3,4] = -G7
    Z[4,3] = -G7
    return Z

# Constant part of nodal analysis matrix
Zconst = np.zeros(shape=(5,5))
Zconst[1,4] = -G4
Zconst[2,3] = -G5
Zconst[3,2] = -G5
Zconst[4,1] = -G4
                  
Zconst[0,0] = G3
Zconst[1,1] = G4 
Zconst[2,2] = G5
Zconst[3,3] = G5
Zconst[4,4] = G4

# Input current vector
J1 = np.zeros(shape=(5,1))
J1[0] = 1

# Limits for the input current and the conductances of the resistors
parameter_limits = [(1,10)] + 4*[(0.1, 10.0)]

def full_model(params):
    """ Calculate a solution for given parameter vector """

    i1 = params[0]
    conds = params[1:]
    
    Z = Zconst + Zmu(conds)
    J = J1*i1
    
    u = np.linalg.solve(Z,J)
    return u

def rom_builder(psi):
    """ Build a reduced order model. A function which 
    given a parameter vector returns a reduced order solution """

    rZconst = psi.T@Zconst@psi
    rJ1 = psi.T@J1
    
    def reduced_model(params):
        i1 = params[0]
        conds = params[1:]
        
        rZ = rZconst + psi.T@Zmu(conds)@psi
        rJ = rJ1*i1
        
        return np.linalg.solve(rZ, rJ)
    
    return reduced_model

def error_estimator(psi, reduced_model, params):
    """ Estimates the error of the reduced model by evaluating 
    norm(Z * psi * u_red - J) which should be zero if reduced order 
    model is perfectly accurate. """ 
    i1 = params[0]
    conds = params[1:]
    
    rsol = reduced_model(params)
    
    Z = Zconst + Zmu(conds)
    J = J1*i1
    
    return np.linalg.norm(Z@psi@rsol - J)

# Make a reducer instance using The Awesome Constructor Function(TM)
red = make_reducer(full_model=full_model, 
                   parameter_limits=parameter_limits,
                   error_estimator=error_estimator,
                   # The "cutoff" dimension, ie. take this many
                   # modes into the reduced order model.
                   poddim=3,
                   # or
                   # The error limit for the cutoff,
                   # podtol = 1e-9,
                   rom_builder=rom_builder,
                   # print stuff
                   verbose=True,
                   # If exceptions are thrown from argument functions, 
                   # rethrow them
                   rethrow_exceptions=True)

# Do an initial reduced order model using 2 samples
# this is very little
res = red.reduce(4)

# Improve the model iteratively
# by finding the maximum of the error function
for i in range(0,100):
    red.improve(res)
    if res.error < 1e-9:
        break
else:
    print("The error didn't converge to required tolerance.")
    print("This happens if poddim is too small or podtol "
          "is too big or too small amount of samples")
    print("Remainder error is {}".format(res.error))
    
# Plot the singular values of the system
plt.plot(np.arange(1,res.ss.shape[0]+1), res.ss)
plt.title('Magnitude of singular values')
plt.xlabel("Index of $\sigma$")
plt.ylabel("Magnitude of $\sigma$")
plt.xticks(np.arange(1,res.ss.shape[0]+1))
plt.show()

print("Original dimension {}, reduced dimension {}".format(res.psi.shape[0],
                                                           res.psi.shape[1]))

