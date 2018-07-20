# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as st

import matplotlib.pyplot as plt
import seaborn as sns

rc = {'lines.linewidth': 2, 
      'axes.labelsize': 18, 
      'axes.titlesize': 18, 
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)

simple_update=np.array([[1,0],[-1,0],[0,1],[0,-1]],dtype=np.int)

def simple_propensity(parameters,population):
    beta_m,beta_p,gamma=parameters
    m,p=population
    return np.array([beta_m,m,beta_p*m,gamma*p])

def discrete_sample_draw(prob):
    return st.rv_discrete(values=(range(len(prob)),prob)).rvs()

def discrete_sample(prob):
    q=np.random.rand()
    i=0
    p_sum=0.0
    while p_sum<q:
        p_sum+=prob[i]
        i+=1   
    return i-1

def gillespie_draw(parameters,propensity_function,population):
    prop=propensity_function(parameters,population)
    prop_sum=prop.sum()
    time=np.random.exponential(1.0/prop_sum)
    rxn_prob=prop/prop_sum
    rxn=discrete_sample(rxn_prob)
    return rxn,time

def gillespie_ssa(parameters,propensity_function,update,population_0,time_points):
    pop_out = np.empty((len(time_points), update.shape[1]), dtype=np.int)
    i_time=1
    i=0
    population_prev=population_0
    t=time_points[0]
    population=population_0.copy()
    pop_out[0,:]=population
    while i<len(time_points):
        while t<time_points[i_time]:
            event,dt=gillespie_draw(parameters,propensity_function,population)
            population_prev=population.copy()
            population+=update[event,:]
            t+=dt
        i=np.searchsorted(time_points>t,True)
        pop_out[i_time:min(i,len(time_points))] = population_prev
        i_time=i
    return pop_out
# Specify parameters for calculation
params = np.array([3, 5, 0.4])
time_points = np.linspace(0, 50, 101)
population_0 = np.array([0, 0])
n_simulations = 50

# Seed random number generator for reproducibility
#np.random.seed(42)

# Initialize output array
pops = np.empty((n_simulations, len(time_points), 2))

# Run the calculations
for i in range(n_simulations):
    pops[i,:,:] = gillespie_ssa(params, simple_propensity, simple_update,
                                population_0, time_points)
    
fig, ax = plt.subplots(1, 2, figsize=(14, 5))

# Plot mRNA trajectories
for i in range(n_simulations):
    ax[0].plot(time_points, pops[i,:,0], '-', lw=0.3, alpha=0.2,
               color=sns.color_palette()[0])

# Plot mRNA mean
ax[0].plot(time_points, pops[:,:,0].mean(axis=0), '-', lw=6, 
           color=sns.color_palette()[2])

# Plot protein trajectories
for i in range(n_simulations):
    ax[1].plot(time_points, pops[i,:,1], 'k-', lw=0.3, alpha=0.2,
               color=sns.color_palette()[0])

# Plot protein mean
ax[1].plot(time_points, pops[:,:,1].mean(axis=0), 'r-', lw=6, 
           color=sns.color_palette()[1])

# Label axes
ax[0].set_xlabel('dimensionless time')
ax[1].set_xlabel('dimensionless time')
ax[0].set_ylabel('number of mRNAs')
ax[1].set_ylabel('number of proteins')
plt.tight_layout()