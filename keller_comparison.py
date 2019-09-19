##
import numpy as np
from photon_correlation_functions import fl_signal
import matplotlib.pyplot as plt

## parameters of the experiment
Omega = 2 * np.pi * 25.42e6  # RF drive freq
decay_rate = 2*np.pi*19.6e6  # units 2*pi * decay freq
#decay_rate = 2*np.pi*25.42e6
beta = 0.085

## plot of DeltaS/S
N = 300 # number of plotted points
laser_detun = np.linspace(-3*decay_rate,0 , N)

phot_signal = fl_signal(beta, laser_detun, Omega, decay_rate)

plt.plot(laser_detun/decay_rate, phot_signal)
plt.show()