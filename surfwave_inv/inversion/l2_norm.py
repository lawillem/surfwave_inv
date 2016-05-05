from residual import residual

def l2_norm(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs):
    
    res = residual(freqs, obs, vp_arr, vs_arr, rho_arr, thk_arr, **kwargs)
    return 0.5*res.dot(res)