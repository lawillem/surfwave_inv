from ctypes import *
import numpy as np
import os

def get_disp_curve(freqs, vp_arr, vs_arr, rho_arr, thick_arr, verbose = 0, c_min = None, c_def_step=10.0, NQUAD=5):
    
    #get_disp_crv(int N, double* alphas, double* betas, double* rhos, double* ds, double* phase_vels, double* freqs, int nfreqs, int verbose)
    
    def load_func(): #Will be called only once. Load the function from the so
        lib = cdll.LoadLibrary(os.path.dirname(__file__) +'/../sk_disp_crv/_sk_disp_crv.so') #os.path.dirname(__file__) gets path of this wrapper function
        func = lib.get_disp_crv
        func.restype  = c_int
        func.argtypes = [c_int,                                                  #N 
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #alphas
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #betas
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #rhos
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #ds
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #phase_vels
                         np.ctypeslib.ndpointer(c_double, flags="C_CONTIGUOUS"), #freqs
                         c_int                                                 , #nfreqs
                         c_double                                              , #C_min
                         c_double                                              , #C_def_step
                         c_int                                                 , #NQUAD
                         c_int                                                 ] #verbose
        return func
    
    try: #See if we already loaded the shared lib and got func setup
        func = get_disp_curve.func
    except: #If not, load and setup
        func = load_func()
        get_disp_curve.func = func

    if c_min == None:
        c_min = 0.5*np.min(vs_arr) #rather conservative

    nfreqs = freqs.size
    N      = vp_arr.size
        
    #sanity check
    if vs_arr.size != N or rho_arr.size != N:
        raise Exception("Inconsistent material input array sizes")

    if thick_arr.size !=N-1:
        if thick_arr.size != N or thick_arr[-1] != np.inf:
            raise Exception("Thickness array not correct. Should either have length N with last value np.inf or size N-1")

    phase_vels = np.zeros(nfreqs, dtype='float64')
    
    try: 
        retval = func(c_int(N), 
                      np.ascontiguousarray(vp_arr),
                      np.ascontiguousarray(vs_arr),
                      np.ascontiguousarray(rho_arr),
                      np.ascontiguousarray(thick_arr),
                      np.ascontiguousarray(phase_vels),
                      np.ascontiguousarray(freqs),
                      c_int(nfreqs),
                      c_double(c_min),
                      c_double(c_def_step),
                      c_int(NQUAD),
                      c_int(verbose)
                      ) 
    except:
        raise Exception('Unexpected error during execution')
    
    if retval == 0:    
        return phase_vels
    else:
        raise Exception("The dispersion curve C-code exited with error code %i"%retval)