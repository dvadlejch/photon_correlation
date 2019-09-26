# 26.9.2019
# calculation of the modulation index beta for a real measured data

#from mat4py import loadmat
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
## data import
data_fluorescence = loadmat('micromotion_99_X_Y_110_bfiel2.mat')['namerene_hodnoty']
len_data = len(data_fluorescence)
## measurement parameters
Omega = 2 * np.pi * 30.151e6  # RF drive freq
decay_rate = 132e6  # units 2*pi * decay freq
#laser_detun = -1 / 2 * decay_rate  # laser detuning
laser_detun = -1/2 * decay_rate # laser detuning

## calculating photon correlation histogram for given data
T = 2*np.pi/Omega # period of the RF drive
time_step = 128e-12 # time step in seconds
N_steps_per_rf_period = int(T/time_step)
N_rf_periods = int(len_data/N_steps_per_rf_period)

# t_rf_cross = T * np.arange(0, int(len_data*time_step/T))  # times of rf zero crossings
# t = time_step * np.arange(0, len_data) # time of measurement

# now I want to divide measurement period into periods of RF and sum all of the fluorescence
cumulative_fluorescence = np.zeros(N_steps_per_rf_period)
for i in range(N_rf_periods):
    for j in range(N_steps_per_rf_period):
        cumulative_fluorescence[j] += data_fluorescence[j+i*N_steps_per_rf_period,1]


# plt.plot(t, data_fluorescence[:,1],'.')
plt.plot(cumulative_fluorescence,'.')
plt.show()
#
plt.figure()
plt.plot(data_fluorescence[:,1],'.')
plt.show()

## fit of the histogram
# S(t) = S_0 + Delta S * cos(Omega*t - phi)  # Keller2015

def fit_resid(x, Omega, S, time_step):
    # x = [S_0, Delta S, phi]
    len_S = len(S)
    S_fit = x[0] + x[1]*np.cos(Omega*time_step*np.arange(0,len_S) - x[2])
    return S - S_fit

fit = least_squares(fit_resid, [320000, 80000, 0], args=(Omega, cumulative_fluorescence, time_step))  # fit.x[0] = S_0, fit.x[1] = Delta S, fit.x[2] = phi

# -----------plot of the fit
S_fit = fit.x[0] + fit.x[1]*np.cos(Omega*time_step*np.arange(0,len(cumulative_fluorescence)) - fit.x[2])
plt.plot(S_fit)
plt.plot(cumulative_fluorescence, '.')
plt.show()
#------------