#I want to do some kind of intelligent grid search to find a solution within the basin of attraction of the global minimum.
#This solution can then be fed to a local method such as levenberg-marquardt in inv_lm.py
#Here I chose basinhopping as there was already a routine in scipy for this
import numpy as np
import functools
from private_routines import l2_norm_wrap
from scipy.optimize import basinhopping

def lsqrs_inv_global(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, niter = 1000, c_min = 100.0, c_def_step = 1.0, verbose=False):
    #freqs: vector of frequencies (Hz), length nfreqs
    #obs  : vector of observations at those frequencies (i.e. the measured dispersion curve), length N
    #vp_arr: Vector of P-wave velocities with length N
    #vs_arr: Vector of S-wave velocities with length N 
    #rho_arr: Vector of densities with length N
    #thk_arr: Vector of layer thicknesses of either length N-1 or length N with the last value np.inf to indicate that the bottom layer is the half-space. Both methods are equivalent.

    #Bounds on where to search, http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.basinhopping.html
    #Will only implement lower bound
    class MyBounds(object): 
        def __init__(self, x_min):
            self.xmin = x_min
        def __call__(self, **kwargs):
            #Will return false if any entry in proposed x is smaller than x_min.
            #This way this step will not be considered.
            x = kwargs["x_new"]
            is_correct = np.any(x < self.xmin)
            
            if not is_correct:
                print "rejecting unphysical step"
            
            return is_correct

    #do not alter input, copy it and use as initial
    vs_init_arr = np.copy(vs_arr)
    
    #Define bounds (prevent searching negative velocities which would crash the dispersion curve generating code)
    #The dispersion curve code starts looking at c_min. Rayleigh speed could be slightly lower than lowest Vs. 
    #Make sure that lowest possible Rayleigh wave speed is above c_min searched in dispersion curve code.
    #Empirical value, not based on research
    empirical_val = 1.2
    min_val       = empirical_val*c_min  
    mybounds      = MyBounds(min_val)
    
    #resid_f will just take a vector vs of length N and return a residual vector of size nfreqs
    resid_f  = functools.partial(l2_norm_wrap, freqs, obs, vp_arr, rho_arr, thk_arr, c_min = c_min, c_def_step=c_def_step, verbose=verbose)

    ret = basinhopping(resid_f, vs_init_arr, niter=niter, accept_test=mybounds)    
    
    return ret.x
    

