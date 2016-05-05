#Inversion using the levenberg-marquardt routine
from scipy.optimize import leastsq #nonlinear least squares package implementing levenberg-marquardt for instance
from private_routines import resid_wrap, jacob_vs_wrap
import functools
import numpy as np

#Should include a sanity check on the combination of vp, vs and rho to make sure it is physical
def lsqrs_inv_lm(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, c_min = 100.0, c_def_step = 1.0, perturb_pct = 0.1, verbose=False):
    #freqs: vector of frequencies (Hz), length nfreqs
    #obs  : vector of observations at those frequencies (i.e. the measured dispersion curve), length N
    #vp_arr: Vector of P-wave velocities with length N
    #vs_arr: Vector of S-wave velocities with length N 
    #rho_arr: Vector of densities with length N
    #thk_arr: Vector of layer thicknesses of either length N-1 or length N with the last value np.inf to indicate that the bottom layer is the half-space. Both methods are equivalent.

    #do not alter input, copy it and use as initial
    vs_init_arr = np.copy(vs_arr)
    
    #resid_f will just take a vector vs of length N and return a residual vector of size nfreqs
    resid_f  = functools.partial(resid_wrap, freqs, obs, vp_arr, rho_arr, thk_arr, c_min = c_min, c_def_step=c_def_step, verbose=verbose)
    
    jacob_vs = functools.partial(jacob_vs_wrap, freqs, obs, vp_arr, rho_arr, thk_arr, c_min = c_min, c_def_step=c_def_step, perturb_pct = perturb_pct, verbose=verbose)
    
    vs_fin, success = leastsq(resid_f,vs_init_arr, Dfun = jacob_vs)
    
    if not success:
        raise Exception("Not succesful, code=%i", success)
    
    return vs_fin