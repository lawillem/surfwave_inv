from surfwave_inv.inversion.disp_crv_wrap import get_disp_curve

def residual(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs):
    #freqs: vector of frequencies (Hz), length nfreqs
    #obs  : vector of observations at those frequencies (i.e. the measured dispersion curve), length N
    #vp_arr: Vector of P-wave velocities with length N
    #vs_arr: Vector of S-wave velocities with length N 
    #rho_arr: Vector of densities with length N
    #thk_arr: Vector of layer thicknesses of either length N-1 or length N with the last value np.inf to indicate that the bottom layer is the half-space. Both methods are equivalent.
    
    #The misfit function
    sim   = get_disp_curve(freqs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs)
    resid = sim-obs
    return resid
            