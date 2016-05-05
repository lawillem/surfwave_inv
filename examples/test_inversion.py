import numpy as np
import matplotlib.pyplot as plt
from surfwave_inv.inversion.inversion import lsqrs_inversion
from surfwave_inv.inversion.disp_crv_wrap import get_disp_curve
from surfwave_inv.plot_helper.profiles import gen_depth_profiles
from get_model import model

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999, use as true
vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true = model()

#Some settings we use. Actually default settings, but leaving them here for easy changing.
verbose      = False
minfreq      = 5.0
maxfreq      = 50.0
freqs        = np.linspace(minfreq, maxfreq, 19)
perturb_pct  = 0.1
c_def_step   = 1.0
niter_global = 125
#'true' dispersion curve. This will be the curve we try to fit
obs = get_disp_curve(freqs, vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true, verbose = verbose)

#Start inversion. Does not need initial guess for Vs, as we start with a global optimization
inv_ret_dict = lsqrs_inversion(freqs, obs, vp_arr_true, rho_arr_true, thk_arr_true, niter_global=niter_global, c_def_step = c_def_step, perturb_pct = perturb_pct, verbose=False)
vs_pso_1     = inv_ret_dict['vs_pso']
vs_lm_1      = inv_ret_dict['vs_lm']

sim_pso_1    = get_disp_curve(freqs, vp_arr_true, vs_pso_1, rho_arr_true, thk_arr_true, verbose = verbose)
sim_lm_1     = get_disp_curve(freqs, vp_arr_true, vs_lm_1, rho_arr_true, thk_arr_true, verbose = verbose)

#Make some depth plots
N                 = vp_arr_true.size
nlayer            = N-1
bdry_depths       = np.cumsum(thk_arr_true[0:nlayer])
deepest_bdry      = bdry_depths[-1]
z_max             = 1.2*deepest_bdry 
z_arr             = np.linspace(0.0, z_max, 1001)
vs_true_profile   = gen_depth_profiles(z_arr, vs_arr_true , thk_arr_true)
vs_pso_1_profile  = gen_depth_profiles(z_arr, vs_pso_1    , thk_arr_true)
vs_lm_1_profile   = gen_depth_profiles(z_arr, vs_lm_1     , thk_arr_true)

#DO PLOTS FOR INVERSION CASE 1
plt.figure(1) #Dispersion
plt.plot(freqs, obs, 'k', label='True')
plt.plot(freqs, sim_pso_1 , 'b', label ='PSO result 1')
plt.plot(freqs, sim_lm_1  , 'r', label ='LM result 1')
plt.legend()
plt.title("Dispersion curve inv case 1. PSO, then LM.")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim([minfreq, maxfreq])
plt.ylim([0.0, 1000.0])

plt.figure(2) #depth profiles
plt.plot(z_arr, vs_true_profile , 'k', label='True')
plt.plot(z_arr, vs_pso_1_profile, 'b', label ='PSO result 1')
plt.plot(z_arr, vs_lm_1_profile , 'r', label ='LM result 1')
plt.legend()
plt.title("Depth profile inv case 1. PSO, then LM.")
plt.xlabel("Depth (m)")
plt.ylabel("Vs (m/s)")
plt.xlim([0.0, z_max])
plt.ylim([0.0, 1000.0])

plt.show()

print "Now using slightly incorrect Vp and Rho"

#arbitrary changes (still true thickness)
vp_arr_est  = 0.87*vp_arr_true + 100.0
rho_arr_est = 1.2*rho_arr_true - 100.0

#Start inversion. Does not need initial guess for Vs, as we start with a global optimization
inv_ret_dict = lsqrs_inversion(freqs, obs, vp_arr_est, rho_arr_est, thk_arr_true, niter_global=niter_global, c_def_step = c_def_step, perturb_pct = perturb_pct, verbose=False)
vs_pso_2     = inv_ret_dict['vs_pso']
vs_lm_2      = inv_ret_dict['vs_lm']

sim_pso_2    = get_disp_curve(freqs, vp_arr_true, vs_pso_2, rho_arr_true, thk_arr_true, verbose = verbose)
sim_lm_2     = get_disp_curve(freqs, vp_arr_true, vs_lm_2, rho_arr_true, thk_arr_true, verbose = verbose)

vs_pso_2_profile  = gen_depth_profiles(z_arr, vs_pso_2    , thk_arr_true)
vs_lm_2_profile   = gen_depth_profiles(z_arr, vs_lm_2     , thk_arr_true)

#DO PLOTS FOR INVERSION CASE 2
plt.figure(3) #Dispersion
plt.plot(freqs, obs, 'k', label='True')
plt.plot(freqs, sim_pso_2 , 'b', label ='PSO result 2')
plt.plot(freqs, sim_lm_2  , 'r', label ='LM result 2')
plt.legend()
plt.title("Dispersion curve inv case 1. PSO, then LM.")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim([minfreq, maxfreq])
plt.ylim([0.0, 1000.0])

plt.figure(4) #depth profiles
plt.plot(z_arr, vs_true_profile , 'k', label='True')
plt.plot(z_arr, vs_pso_2_profile, 'b', label ='PSO result 2')
plt.plot(z_arr, vs_lm_2_profile , 'r', label ='LM result 2')
plt.legend()
plt.title("Depth profile inv case 1. PSO, then LM.")
plt.xlabel("Depth (m)")
plt.ylabel("Vs (m/s)")
plt.xlim([0.0, z_max])
plt.ylim([0.0, 1000.0])

plt.show()