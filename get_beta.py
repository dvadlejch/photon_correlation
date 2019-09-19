# 19.9.2019
# code which is intended to calculate value of excess micromotion index beta given the photon-correlation signal

##
import numpy as np
from scipy.optimize import fsolve
from photon_correlation_functions import fl_signal

## parameters of the experiment
Omega = 2 * np.pi * 30.141e6  # RF drive freq
decay_rate = 132e6  # units 2*pi * decay freq
laser_detun = -1 / 2 * decay_rate  # laser detuning

## values from experiment
norm_mod_amp = 0.01  # normalized modulation amplitude determined by experiment


## solution

def root_func(beta, laser_detun, Omega, decay_rate, norm_mod_amp):
    return fl_signal(beta, laser_detun, Omega, decay_rate) - norm_mod_amp


sol = fsolve(root_func, np.array([0.5]), args=(laser_detun, Omega, decay_rate, norm_mod_amp), full_output=True)

print(sol)
