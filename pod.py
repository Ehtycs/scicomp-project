#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 20:47:53 2018

@author: marjamaa

References:
1. Marjamäki, A., Coupling reduced order dynamic electromagnetic models with 
   circuit simulators, Available http://urn.fi/URN:NBN:fi:tty-201708241744

"""

import numpy as np
#import scipy.sparse as sps
import scipy.optimize as spopt
#from itertools import product, repeat
#from collections import namedtuple
#import functools as ft #, partial
#from operator import mul

def make_reducer(**kwargs):
    """ Common factory method to construct reducers
    Maybe in the future there are other types available"""
        
    reducer = PodReducer()

    # Take the kwargs and check if there are unknown kwargs provided
    # I was planning to extend this to check if some required args are missing
    # but that's probably not the way to go.
    kwkeys = set(kwargs.keys())
    # Reduced allowed_attrs is a frozen set of allowed attributes
    unkn = kwkeys.difference(reducer.allowed_attrs)
    if(len(unkn) > 0):
        raise ValueError("Unknown parameters supplied:"+str(unkn))

    # Set the attributes, eliminates the need for typing a constructor
    # method which just takes these and assigns them
    for key in kwkeys.intersection(reducer.allowed_attrs):
        setattr(reducer, key, kwargs[key])

    return reducer


def np_extend_to(arr, dim, const):
    """ Extend numpy array with zeros"""
    adim = arr.shape[0]
    return np.concatenate((arr, np.full(dim-adim, const)), axis=0)

def random_sampling(limits, n):
    """ Create n random parameter vectors from parameter range plimits"""
    ub,lb = np.array(limits).T
    return np.random.rand(n, ub.shape[0])*(ub-lb) + lb

def boundary_basis_sampling(limits):
    """ Create a sampling so that each one of the limits is "one" 
    at a time. The rest are "zero" 
    """
    eve = np.eye(len(limits))
    ub,lb = np.array(limits).T
    return (1-eve)*ub + eve*lb


#def calc_pset_dimensions(parameters):
#    """ Calculate the dimension and the amount of elements in the 
#    discrete parameter space
#    """
#    dim = len(parameters)
#    nelems = reduce(mul, (len(p) for p in parameters))
#    return (dim, nelems)

def check_null(val, msg):
    if val is None:
       raise ValueError(msg)
       
#class ParameterSpace(object):
#    
#    def __init__(self, psubs):
#        self.psubspaces = psubs
#       
#    def _buildSampleGenerator(self, psubs):
#        pass

       
class ReductionResult(object):
    """ Stores the reduction result so that it can be refined. """
    def __init__(self, psi, ss=None, snaps=None, err=None, pconf=None):
        self.psi = psi
        self.ss = ss
        self.snapshots = snaps
        self.error = err
        self.last_pconf = pconf

class PodReducer(object):
    
    allowed_attrs = frozenset([
            'full_model', 'verbose', 'podtol', 'poddim',
            'initial_parameters', 'parameter_limits',
            'rom_builder', 'error_estimator',
            'rethrow_exceptions'])
    
    def __init__(self):
        """ Populate default values """
        self.verbose = False
        self.endtol = 1e-6
        self.podtol = 1e-6
        self.poddim = None
        
        self.initial_parameters = None
        self.parameter_limits = None
        
        self.rethrow_exceptions = True
        
        self.rom_builder = None
        self.full_model = None
        self.error_estimator = None
        
    def check_req(self):
        """ Check requirements for standard reduce operation """
        if(self.initial_parameters is None):
            check_null(self.parameter_limits, 
                       "Parameter_limits and initial_parameters both"
                       "missing, provide one")
        
    def check_req_improve(self):
        """ Check requirements for improve operation """
        check_null(self.parameter_limits, "parameter_limits missing")
        check_null(self.rom_builder, "rom_builder missing")
        check_null(self.error_estimator, "error_estimator missing")
        
        
    def calculate_projector(self, snaps):
        """ Calculate the projector based on snapshots 
            Returns the projector and the singular values 
        """
        
        u, s, _ = np.linalg.svd(np.stack(snaps, axis=1).squeeze(), 
                                full_matrices=False)
        
        # cut off criteria:
        # [1] page 24. equation (3.7)
        if(self.poddim is None):
            mask = np.flip(np.cumsum(np.flip(s**2,0)), 0) > self.podtol 
            psi = (u @ np.diag(s))[:, mask]
        else:
            # fixed dimension cut-off
            psi = (u @ np.diag(s))[:,0:self.poddim]
        
        return psi, s

    def evaluate_fullmodel(self, param):
        """ Evaluate the full model and handle possible errors """ 
        try:
            return self.full_model(param)
        except Exception as e:
            self.msg("User provided full model"
                      " threw an exception: {}".format(e))
            if self.rethrow_exceptions:
                raise e
        return None

    def msg(self, msg):
        if(self.verbose):
            print(msg)
        
    def calculate_snapshots(self, parameters):
        """ calculate a bunch of snapshots """
        snaps = (self.evaluate_fullmodel(ps) 
                for ps in parameters)
        return [s for s in snaps if s is not None]
    
    def reduce(self, ninitp=10):
        """ Do the standard reduce :
            1. Take the parameter space and sample it, 
            2. calculate snapshots for each sample
            3. SVD and return the reduction result 
        """
        
        self.check_req()
        
        plims = self.parameter_limits

        if(self.initial_parameters is not None):
            params = self.initial_parameters
        else:
            params = random_sampling(plims, ninitp)

        self.msg("Calculating snapshots...")        

        snaps = self.calculate_snapshots(params)

        self.msg("Calculating projection mapping...")

        psi, ss = self.calculate_projector(snaps)

        return ReductionResult(psi, ss, snaps)
        
    def targetfun(self, pspacev, psi, reduced_model):
        """ The target function to minimize using scipy opt. 
        Maximise function == minimize -function """
        return -self.error_estimator(psi, reduced_model, pspacev)
    
    def initial_guess(self):
        """ Randomize an initial guess from the parameter space """
        lims = self.parameter_limits
        ub,lb = np.array(lims).T
        rans = np.random.rand(1, ub.shape[0])*(ub-lb) + lb
        return rans


    def improve(self, red_result_mut):
        """ Improve an existing reduced order model.
        
        Replaces the old model in red_redsult_mut.
        1. Get the old model, run a maximization algorithm on error 
           function.
        2. Evaluate the full model on the maximum, add result to snapshots
        3. Calculate new reduced order model.
        
        Return nothing to kind of signal that this function does naughty 
        mutations for its input parameter.
        """ 
        
        psi = red_result_mut.psi
        snaps = red_result_mut.snapshots    
        
        # If there exists a previous maximum, use that. If not randomize new
        if red_result_mut.last_pconf is None:
            pconf = self.initial_guess()
        else:
            pconf = red_result_mut.last_pconf
        
        redmodel = self.rom_builder(psi)

        for i in range(0,10):
            """ Sometimes this spopt.minimize hits a parameter value 
            which is not solvable (if stiffness matrix becomes singular etc...)
            Then this step retries with different initial guess until 
            things work out or fails after 10 tries. """
            try:
                           
                res = spopt.minimize(self.targetfun, pconf, 
                                     args=(psi, redmodel), 
                                     bounds=self.parameter_limits)
                break
            
            except Exception as e:
                self.msg("Minimization step threw an exception: {}".format(e))
                if self.rethrow_exceptions:
                    raise e
                # Randomize a new point
                pconf = self.initial_guess()
        else:
            raise RuntimeError("Minimization failed after 10 tries.")
        
        #print(res.nfev)
        #print(res.x)
        snap = self.evaluate_fullmodel(res.x)
    
        if(snap is not None):
            snaps.append(snap)
        
        psi,ss = self.calculate_projector(snaps)
        
        # Shamefully mutate the input value
        red_result_mut.psi = psi
        red_result_mut.ss = ss
        #red_result_mut.snapshots = snaps
        red_result_mut.error = -res.fun 
        red_result_mut.last_pconf = res.x
    