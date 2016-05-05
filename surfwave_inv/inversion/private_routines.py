from jacobian import jacobian_vs
from residual import residual
from l2_norm import l2_norm

#Some private convenience routines for the inversion functions
def l2_norm_wrap(freqs, obs, vp_arr, rho_arr, thk_arr, vs_arr, **kwargs):
    #vs is the model we want to invert for. 
    #when using partial from functools, can only prepend on left
    #do not want to modify original residual definition to keep order of arguments consistent
    
    import matplotlib.pyplot as plt
    plt.plot(vs_arr)
    plt.show()
        
    return l2_norm(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs)

def resid_wrap(freqs, obs, vp_arr, rho_arr, thk_arr, vs_arr, **kwargs):
    #vs is the model we want to invert for. 
    #when using partial from functools, can only prepend on left
    #do not want to modify original residual definition to keep order of arguments consistent
        
    return residual(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs)

def jacob_vs_wrap(freqs, obs, vp_arr, rho_arr, thk_arr, vs_arr, **kwargs):
    #same reason for wrap as in resid_wrap
        
    return jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs) 