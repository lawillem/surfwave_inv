#Compute the dispersion spectrum on an input shotgather
#Code adapted from matlab code provided by Elita (Yunyue) Li  

import numpy as np

def get_disp_crv(data, ts, xs, fmin=4.0, fmax=50.0, cmin=10.0, cmax=400.0, dc = 5.0):
    #data: 2d shotgather with nt rows and xs cols
    #ts  : 1d array of times. 
    #xs  : 1d array of receiver locations
    
    def float_eq(a, b, epsilon=0.00000001): #convenience function to test for float equality. 
        return abs(a - b) < epsilon
    

    #First check x and t arrays are increasing to the right
    #for now enforce, could also correct
    if xs[-1] < xs[0] or ts[-1] < ts[0]:
        raise Exception("Assuming that ts and xs increase to the right")
    
    
    dt = ts[1] - ts[0]
    dx = xs[1] - xs[0]
    
    nt = ts.size
    nx = xs.size

    if fmax >= 0.5*1./dt:
        raise Exception("Requested freq above nyquist")
    
    #currently assuming ts and xs uniformly spaced. 
    #verify

    if not np.all(float_eq(ts, dt*np.arange(nt))):
        raise Exception("Time vector not uniformly spaced. Require for now.")
    
    if not np.all(float_eq(xs, dx*np.arange(nx))):
        raise Exception("Position vector not uniformly spaced. Require for now.")    
    
    #Done with checking
    xs = xs - xs[np.floor(nx/2)]

    #Do FFT in time
    data_fft  = np.fft.fft(data, axis=0)
    freqs_all = np.fft.fftfreq(nt, dt)
    
    freq_inds = np.logical_and(freqs_all >= fmin, freqs_all <= fmax)
    freqs     = freqs_all[freq_inds]
    nf        = freqs.size
    
    nc        = int(np.round((cmax-cmin)/dc))
    cs        = cmin + dc*np.arange(nc)
    
    
    subs_data_fft = data_fft[freq_inds,:]
    
    #slant stack to get the phase velocity spectrum
    phsspec   = np.zeros((nc, nf), dtype='complex128')
    for ic in xrange(nc):
        for ifr in xrange(nf):
            freq = freqs[ifr]
            c    = cs[ic]
            phsspec[-(ic+1), ifr] = np.sum(np.exp(1j*2*np.pi*freq*xs/c)*subs_data_fft[ifr,:]/np.abs(subs_data_fft[ifr,:]) ) 
    
    
    
    return phsspec, cs, freqs

    