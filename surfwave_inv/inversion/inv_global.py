#I want to do some kind of intelligent grid search to find a solution within the basin of attraction of the global minimum.
#This solution can then be fed to a local method such as levenberg-marquardt in inv_lm.py
#Here I chose particle swarm optimization (PSO) as there was a routine available for this in DEAP
import functools
from private_routines import l2_norm_wrap

#PSO imports
#http://deap.readthedocs.io/en/master/examples/pso_basic.html
#numpy version at: https://github.com/DEAP/deap/blob/master/examples/pso/basic_numpy.py
import operator
import random

import numpy as np

from deap import base
from deap import benchmarks
from deap import creator
from deap import tools

def lsqrs_inv_global(freqs, obs, vp_arr, rho_arr, thk_arr, npart = 0, niter = 1000, c_min = 100.0, c_def_step = 1.0, verbose=False):
    #freqs: vector of frequencies (Hz), length nfreqs
    #obs  : vector of observations at those frequencies (i.e. the measured dispersion curve), length N
    #vp_arr: Vector of P-wave velocities with length N
    #rho_arr: Vector of densities with length N
    #thk_arr: Vector of layer thicknesses of either length N-1 or length N with the last value np.inf to indicate that the bottom layer is the half-space. Both methods are equivalent.
    #npart: number of particles
    #niter: number of evolution steps for each particle
    
    #total number of dispersion curves will therefore by (npart * niter)

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Particle", np.ndarray, fitness=creator.FitnessMax, speed=list, 
        pmin=list, pmax=list, smin=None, smax=None, best=None)    

    #Generate particles
    def generate(size, pmin_arr, pmax_arr, smin, smax):
        part = creator.Particle(np.zeros(size)) #initialize 
        for ip in xrange(size): #generate initial position in each dimension based on bounds
            part[ip] = np.random.uniform(pmin_arr[ip], pmax_arr[ip])

        part.pmin  = pmin_arr
        part.pmax  = pmax_arr         
        part.speed = np.random.uniform(smin, smax, size)
        part.smin  = smin
        part.smax  = smax
        return part

    def updateParticle(part, best, phi1, phi2):
        #generate two random variables
        #right now the velocity of particles for any layer is same
        #should probably bias so that it is larger in the bottom and smaller on the top 
        #as there is naturally more variation in bottom

        if part.best is None or best is None: 
            #the forward model is sometimes problematic for the randomly generated models. 
            #In those cases we need to generate a new model.
            #If this happens during the first iteration, part.best or best may not have been set yet
            #randomly get new value
            for ip in xrange(len(part)): #generate initial position in each dimension based on bounds
                part[ip] = np.random.uniform(part.pmin[ip], part.pmax[ip])
                 
        else: #normal particle update
            u1 = np.random.uniform(0, phi1, len(part))
            u2 = np.random.uniform(0, phi2, len(part))
            
            v_u1 = u1 * (part.best - part)
            v_u2 = u2 * (best - part)
            part.speed += v_u1 + v_u2
            
            #check speed
            for i, speed in enumerate(part.speed):
                if speed < part.smin:
                    part.speed[i] = part.smin
                elif speed > part.smax:
                    part.speed[i] = part.smax
                    
            part += part.speed
            
            #check if position is within bounds. If not, correct.
            for i, pos in enumerate(part):
                if pos < part.pmin[i]:
                    part[i] = part.pmin[i]
                elif pos > part.pmax[i]:
                    part[i] = part.pmax[i]
                
    def particle_eval_wrap(individual, func):
        vs_arr    = np.zeros(len(individual))
        vs_arr[:] = individual[:] #copy data
        
        l2_norm = func(vs_arr)
        
        #the worse the estimate is, the higher the l2 norm is
        #we want to assign a low fitness to this, so use inverse
        retval  = 1.0/l2_norm
        return retval,
    
    
    N = vp_arr.size
    if npart == 0: #if it has not been initialized, set to number of layers 
        npart = N

    #Define bounds (prevent searching negative velocities which would crash the dispersion curve generating code)
    #The dispersion curve code starts looking at c_min. Rayleigh speed could be slightly lower than lowest Vs. 
    #Make sure that lowest possible Rayleigh wave speed is above c_min searched in dispersion curve code.
    #Empirical value, not based on research
    empirical_val = 1.4
    min_val_arr = empirical_val*c_min*np.ones_like(vp_arr)  
    max_val_arr = (1./np.sqrt(2))*vp_arr #any larger than this and elastic constant lambda becomes negative. Vs at 1/sqrt(2)*Vp will already be unrealistics though

    max_part_speed = 25.0
    phi_amp        = 2.0
    #l2_norm_wrap will just take a vector vs of length N and return a scalar
    l2_norm_func = functools.partial(l2_norm_wrap, freqs, obs, vp_arr, rho_arr, thk_arr, c_min = c_min, c_def_step=c_def_step, verbose=verbose)
    
    #init toolbox
    toolbox = base.Toolbox()
    toolbox.register("particle", generate, size=N, pmin_arr=min_val_arr, pmax_arr=max_val_arr, smin=-max_part_speed, smax=max_part_speed)
    toolbox.register("population", tools.initRepeat, list, toolbox.particle)
    toolbox.register("update", updateParticle, phi1=phi_amp, phi2=phi_amp)
    #toolbox.register("evaluate", l2_norm)    
    
    #main part of code
    pop = toolbox.population(n=npart)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    logbook = tools.Logbook()
    logbook.header = ["gen", "evals"] + stats.fields

    GEN = niter
    best = None

    for g in range(GEN):
        print "Particle swarm iteration %i of %i"%(g+1, GEN)
        for part in pop:
            #init
            calc_success = False
            while not calc_success:
                try:
                    part.fitness.values = particle_eval_wrap(part, l2_norm_func)
                    calc_success = True
                except: #very weird profile, makes C code crash (still poorly behaved sometimes)
                    toolbox.update(part, best) #Update and try again
                
            if part.best is None or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if best is None or best.fitness < part.fitness:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best)

        # Gather all the fitnesses in one list and print the stats
        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream)
    
    return best
    

