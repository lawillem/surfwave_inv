import numpy as np
import matplotlib.pyplot as plt
from surfwave_inv.inversion.inv_lm import lsqrs_inv_lm
from surfwave_inv.inversion.disp_crv_wrap import get_disp_curve
from get_model import model

#table 1 in paper "Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves, Xia et al 1999, use as true
vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true = model()

verbose    = False
minfreq    = 5.0
maxfreq    = 50.0
freqs      = np.linspace(minfreq, maxfreq, 19)
perturb_pct= 0.1
c_def_step = 1.0
#'true' dispersion curve
obs = get_disp_curve(freqs, vp_arr_true, vs_arr_true, rho_arr_true, thk_arr_true, verbose = verbose)

#get initial model 1 by perturbing only the vs and keeping vp, rho and thk_arr at true values.
#somewhat arbitrary new vs
vs_arr_init_1      = np.copy(vs_arr_true)
vs_arr_init_1     *= 0.9
vs_arr_init_1[-1] *= 1.3

sim_init_1   = get_disp_curve(freqs, vp_arr_true, vs_arr_init_1, rho_arr_true, thk_arr_true, verbose = verbose)
vs_arr_inv_1 = lsqrs_inv_lm(freqs, obs, vp_arr_true, vs_arr_init_1, rho_arr_true, thk_arr_true, c_def_step = c_def_step, perturb_pct = perturb_pct, verbose=False)
sim_inv_1    = get_disp_curve(freqs, vp_arr_true, vs_arr_inv_1, rho_arr_true, thk_arr_true, verbose = verbose)

plt.figure(1) #Dispersion
plt.plot(freqs, obs, 'r', label='True')
plt.plot(freqs, sim_init_1, 'b', label ='Init 1')
plt.plot(freqs, sim_inv_1 , 'g', label ='Inv  1')
plt.legend()
plt.title("Dispersion curve")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase velocity (m/s)")
plt.xlim([minfreq, maxfreq])
plt.ylim([0.0, 1000.0])
plt.show()
