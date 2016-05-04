import numpy as np

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999
def model():
    vp_arr  = np.array([650.0,
                        750.0,
                        1400.0,
                        1800.0,
                        2150.0,
                        2800.0]) #units m/s
      
    vs_arr  = np.array([194.0,
                        270.0,
                        367.0,
                        485.0,
                        603.0,
                        740.0]) #units m/s
      
    rho_arr = np.array([1820.0,
                        1860.0,
                        1910.0,
                        1960.0,
                        2020.0,
                        2090.0]) #units: kg/m**3
     
     
    thk_arr = np.array([2.0,
                        2.3,
                        2.5,
                        2.8,
                        3.2,
                        np.inf]) #units: m    
    
    return vp_arr, vs_arr, rho_arr, thk_arr