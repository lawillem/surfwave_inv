import numpy as np
import matplotlib.pyplot as plt
from surfwave_inv.inversion.disp_crv_wrap import get_disp_curve

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999
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

verbose    = False 
c_min      = 100.0
c_def_step = 10.0
minfreq    = 5.0
maxfreq    = 50.0
freqs      = np.linspace(minfreq, maxfreq, 19)

# #SINCE WE ARE DEALING WITH UNITLESS QUANTITIES (EPSILONS, ZETAS) WE SHOULD BE ABLE TO SCALE FROM m -> km AND GET SAME RESULT
# scaling = 1e-0 #All terms involving m are scaled. Should cancel out in the constants as they do not have units m
# vp_arr *= scaling
# vs_arr *= scaling
# rho_arr *= scaling
# thk_arr *= scaling
# c_min *= scaling
# c_def_step *= scaling

phase_vels  = get_disp_curve(freqs, vp_arr, vs_arr, rho_arr, thk_arr, verbose = verbose, c_min=c_min, c_def_step=c_def_step)

plt.plot(freqs, phase_vels,  'r', marker='*')
plt.title("Dispersion curve")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim([minfreq, maxfreq])
plt.ylim([0.0, 1000.0])
plt.show()