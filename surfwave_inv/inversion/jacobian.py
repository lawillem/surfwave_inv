#In this file we will estimate the Jacobian
#Right now I only implement the Jacobian for Vs, since Vs has the largest influence on the Rayleigh wave dispersion
from disp_crv_wrap import get_disp_curve
import numpy as np

#Should add functionality that if error code 3 is thrown from C code, smaller step sizes are automatically going to be tried from that point on
#Such a procedure is currently not implemented.
def jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, perturb_pct = 0.1, verbose=False, **kwargs):
    #perturb_pct is the percentage change in amplitude of the vs perturbation we use when approximating J 
    
    if perturb_pct <= 0.0:
        raise Exception("Unacceptable perturbation percentage provided")
    
    nl = vp_arr.size 
    nf = freqs.size
    
    J = np.zeros((nf, nl))
    disp_base = get_disp_curve(freqs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs) 
    for ilayer in xrange(nl): #Fill one col at a time
        if verbose: print "Starting with column %i of Jacobian"%ilayer
        
        vs_pert_arr  = np.copy(vs_arr)
        old_vs_layer = vs_arr[ilayer]
        vs_perturb   = (perturb_pct/100.0)*old_vs_layer
        vs_pert_arr[ilayer] += vs_perturb
        
        disp_pert = get_disp_curve(freqs, vp_arr, vs_pert_arr, rho_arr, thk_arr, verbose=verbose, **kwargs)
        derivs = (disp_pert-disp_base)/vs_perturb
        J[:, ilayer] = derivs
        
    return J        
