## imports
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from photon_correlation_functions import get_beta, fl_signal_second_order

## measurement parameters
Omega = 2 * np.pi * 30.03e6  # RF drive freq
decay_rate = 132e6  # units 2*pi * decay freq
# laser_detun = -1 / 2 * decay_rate  # laser detuning
laser_detun = -1/2  * decay_rate  # laser detuning

## calculating photon correlation histogram for given data
T = 2 * np.pi / Omega  # period of the RF drive
time_step = 128e-12  # time step in seconds
N_steps_per_rf_period = int(T / time_step)
t = time_step * np.arange(0, N_steps_per_rf_period)

S = fl_signal_second_order(0.2, laser_detun, Omega, decay_rate, t, -3.41e-2, 6.26e-1, 1e22)

plt.plot(t, S)
plt.show()