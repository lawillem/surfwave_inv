import numpy as np
import matplotlib.pyplot as plt
from surfwave_inv.inversion.jacobian import jacobian_vs
from get_model import model

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999
vp_arr, vs_arr, rho_arr, thk_arr = model()

verbose    = False 
c_min      = 100.0
c_def_step = 10.0
minfreq    = 5.0
maxfreq    = 50.0
freqs      = np.linspace(minfreq, maxfreq, 19)

perturb_pct_1 = 50.0
perturb_pct_2 = 25.0 #in paper this seems to be suggested
perturb_pct_3 = 10.0

J1  = jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, verbose = verbose, perturb_pct = perturb_pct_1, c_min=c_min, c_def_step=c_def_step)
J2  = jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, verbose = verbose, perturb_pct = perturb_pct_2, c_min=c_min, c_def_step=c_def_step)
J3  = jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, verbose = verbose, perturb_pct = perturb_pct_3, c_min=c_min, c_def_step=c_def_step)

c_def_step_fact = 0.1
c_def_step *= c_def_step_fact
J2_fine = jacobian_vs(freqs, vp_arr, vs_arr, rho_arr, thk_arr, verbose = verbose, perturb_pct = perturb_pct_2, c_min=c_min, c_def_step=c_def_step)

plt.figure(1)
plt.title("%e pct change"%perturb_pct_1)
plt.imshow(J1,interpolation='nearest')
plt.colorbar()
plt.figure(2)
plt.title("%e pct change"%perturb_pct_2)
plt.imshow(J2,interpolation='nearest')
plt.colorbar()
plt.figure(3)
plt.title("%e pct change"%perturb_pct_3)
plt.imshow(J3,interpolation='nearest')
plt.colorbar()
plt.figure(4)
plt.title("%e pct change FINER"%perturb_pct_2)
plt.imshow(J2_fine,interpolation='nearest')
plt.colorbar()
plt.show()