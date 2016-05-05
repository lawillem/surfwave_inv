import numpy as np
import matplotlib.pyplot as plt
from surfwave_inv.inversion.inversion import lsqrs_inversion
from surfwave_inv.inversion.disp_crv_wrap import get_disp_curve
from surfwave_inv.plot_helper.profiles import gen_depth_profiles
from get_model import model

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999, use as true
vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true = model()

#Some settings we use. Actually default settings, but leaving them here for easy changing.
verbose    = True
minfreq    = 5.0
maxfreq    = 50.0
freqs      = np.linspace(minfreq, maxfreq, 19)
perturb_pct= 0.1
c_def_step = 1.0

#'true' dispersion curve. This will be the curve we try to fit
obs = get_disp_curve(freqs, vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true, verbose = verbose)

#get initial model 1 by perturbing only the vs and keeping vp, rho and thk_arr at true values.
#somewhat arbitrary new vs
vs_arr_init_1      = np.copy(vs_arr_true)
vs_arr_init_1     *= 0.8
vs_arr_init_1[-1]  = vs_arr_true[-1]*1.2

sim_init_1   = get_disp_curve(freqs, vp_arr_true, vs_arr_init_1, rho_arr_true, thk_arr_true, verbose = verbose)
vs_arr_inv_1 = lsqrs_inversion(freqs, obs, vp_arr_true, vs_arr_init_1, rho_arr_true, thk_arr_true, c_def_step = c_def_step, perturb_pct = perturb_pct, verbose=False)
sim_inv_1    = get_disp_curve(freqs, vp_arr_true, vs_arr_inv_1, rho_arr_true, thk_arr_true, verbose = verbose)

#Make some depth plots
N                 = vp_arr_true.size
nlayer            = N-1
bdry_depths       = np.cumsum(thk_arr_true[0:nlayer])
deepest_bdry      = bdry_depths[-1]
z_max             = 1.2*deepest_bdry 
z_arr             = np.linspace(0.0, z_max, 1001)
vs_true_profile   = gen_depth_profiles(z_arr, vs_arr_true  , thk_arr_true)
vs_init_1_profile = gen_depth_profiles(z_arr, vs_arr_init_1, thk_arr_true)
vs_inv_1_profile  = gen_depth_profiles(z_arr, vs_arr_inv_1 , thk_arr_true)

#DO PLOTS FOR INVERSION CASE 1
plt.figure(1) #Dispersion
plt.plot(freqs, obs, 'r', label='True')
plt.plot(freqs, sim_init_1, 'b', label ='Init 1')
plt.plot(freqs, sim_inv_1 , 'g', label ='Inv  1')
plt.legend()
plt.title("Dispersion curve inv case 1")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim([minfreq, maxfreq])
plt.ylim([0.0, 1000.0])

plt.figure(2) #depth profiles
plt.plot(z_arr, vs_true_profile, 'r', label='True')
plt.plot(z_arr, vs_init_1_profile, 'b', label ='Init 1')
plt.plot(z_arr, vs_inv_1_profile , 'g', label ='Inv  1')
plt.legend()
plt.title("Depth profile inv case 1")
plt.xlabel("Depth (m)")
plt.ylabel("Vs (m/s)")
plt.xlim([0.0, z_max])
plt.ylim([0.0, 1000.0])

print "MAYBE DO GRID SEARCH FIRST TO GET CLOSE TO TRUE, THEN WHEN WE ARE IN BASIN OF ATTRACTION, USE GAUSS NEWTON"

plt.show()
