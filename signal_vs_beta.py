# 19.9.2019
# code which is intended to calculate value of excess micromotion index beta given the photon-correlation signal

##
import numpy as np
from scipy.optimize import fsolve
from photon_correlation_functions import fl_signal
import matplotlib.pyplot as plt
from scipy.special import j0, j1

## parameters of the experiment
Omega = 2 * np.pi * 14e6  # RF drive freq
decay_rate = 2 * np.pi * 10e6  # units 2*pi * decay freq
# laser_detun = -1/2 * decay_rate  # laser detuning
laser_detun = - 1 / 2 * decay_rate

## values from experiment
norm_mod_amp = 0.01  # normalized modulation amplitude determined by experiment

## signal vs beta plot

N = 300
beta = np.linspace(0, 0.7, N)
ph_signal = fl_signal(beta, laser_detun, Omega, decay_rate)

plt.plot(beta, ph_signal)
plt.show()
