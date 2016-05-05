from surfwave_inv.inversion.inv_global import lsqrs_inv_global
from inv_lm import lsqrs_inv_lm
import numpy as np

def lsqrs_inversion(freqs, obs, vp_arr, rho_arr, thk_arr, npart_global = 10, niter_global = 100, c_min = 100.0, c_def_step = 1.0, perturb_pct = 0.1, verbose=False):
    #freqs: vector of frequencies (Hz), length nfreqs
    #obs  : vector of observations at those frequencies (i.e. the measured dispersion curve), length N
    #vp_arr: Vector of P-wave velocities with length N
    #vs_arr: Vector of S-wave velocities with length N 
    #rho_arr: Vector of densities with length N
    #thk_arr: Vector of layer thicknesses of either length N-1 or length N with the last value np.inf to indicate that the bottom layer is the half-space. Both methods are equivalent.
    
    #First do a global search to get close to global minimum (basin of attraction)
    vs_pso_arr = lsqrs_inv_global(freqs, obs, vp_arr, rho_arr, thk_arr, niter = niter_global, c_min = c_min, c_def_step = c_def_step, verbose=verbose)

    #Use result as a start for local search through levenberg-marquardt
    try: 
        vs_ret_arr = lsqrs_inv_lm(freqs, obs, vp_arr, vs_pso_arr, rho_arr, thk_arr, c_min = c_min, c_def_step = c_def_step, perturb_pct = perturb_pct, verbose=verbose)
        lm_success = True
    except: #If gauss newton fails
        vs_ret_arr = vs_pso_arr
        lm_success = False
        
    retdict = {'vs_pso': vs_pso_arr,
               'vs_lm': vs_ret_arr,
               'lm_success': lm_success}
    
    return retdict