# 27.9.2019
# calculation of the modulation index beta for a real measured data
## import
# from mat4py import loadmat
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from photon_correlation_functions import get_beta_second, fl_signal_second_order
## data import
data_fluorescence = loadmat('micromotion_60_900_900_200.mat')['namerene_hodnoty']
#data_fluorescence = data_fluorescence[:,0:5]
len_data = len(data_fluorescence) # number of fluorescence data in each measurement of micromotion
N_of_data_sets = data_fluorescence.shape[1] # number of measurements of micromotion
## measurement parameters
Omega = 2 * np.pi * 30.03e6  # RF drive freq
decay_rate = 132e6  # units 2*pi * decay freq
# laser_detun = -1 / 2 * decay_rate  # laser detuning
laser_detun = -1/2  * decay_rate  # laser detuning

## calculating photon correlation histogram for given data
T = 2 * np.pi / Omega  # period of the RF drive
time_step = 128e-12  # time step in seconds
N_steps_per_rf_period = int(T / time_step)
N_rf_periods = int(len_data / N_steps_per_rf_period)

# t_rf_cross = T * np.arange(0, int(len_data*time_step/T))  # times of rf zero crossings
# t = time_step * np.arange(0, len_data) # time of measurement

# now I want to divide measurement period into periods of RF and sum all of the fluorescence
cumulative_fluorescence = np.zeros( (N_steps_per_rf_period, N_of_data_sets) )
S_0 = np.zeros(N_of_data_sets)
for k in range(N_of_data_sets):
    for i in range(N_rf_periods):
        for j in range(N_steps_per_rf_period):
            cumulative_fluorescence[j, k] += data_fluorescence[j + i * N_steps_per_rf_period, k]

    S_0[k] = cumulative_fluorescence[:, k].mean()
    cumulative_fluorescence[:, k] = cumulative_fluorescence[:, k] - S_0[k]
# plt.figure(1)
# plt.plot(data_fluorescence,'.')
# plt.show()
# plt.figure(2)
# plt.plot(cumulative_fluorescence, '.')
# plt.show()
# #
# plt.figure()
# plt.plot(data_fluorescence[:, 1], '.')
# plt.show()


## fit of the histogram
# S(t) = S_0 + Delta S * cos(Omega*t - phi)  # Keller2015

def fit_resid(x, Omega, S, time_step):
    # x = [Delta S, phi, Delta_S_2, phi_prime]
    len_S = len(S)
    S_fit = x[0] * np.cos(Omega * time_step * np.arange(0, len_S) + x[1]) + x[2] * np.cos(2*Omega * time_step * np.arange(0, len_S) - x[3])
    return S - S_fit

norm_mod_amp_second = np.zeros(N_of_data_sets)
for k in range(N_of_data_sets):
    # fit = least_squares(fit_resid, [320000, 80000, 0], args=(
    # Omega, cumulative_fluorescence[:,k], time_step), bounds=([0, 0, 0], [np.inf, np.inf, 2*np.pi]) )  # fit.x[0] = S_0, fit.x[1] = Delta S, fit.x[2] = phi
    fit = least_squares(fit_resid, [30000, 0, 1000, 0], args=(Omega, cumulative_fluorescence[:, k], time_step) )  # fit.x[0] = S_0, fit.x[1] = Delta S, fit.x[2] = phi
    print(fit)
    norm_mod_amp_second[k] = np.abs(fit.x[2] / fit.x[0]) # Delta S_2/Delta S

    S_fit = fit.x[0] * np.cos(Omega * time_step * np.arange(0, N_steps_per_rf_period) - fit.x[1]) + fit.x[2] * np.cos(2*Omega * time_step * np.arange(0, N_steps_per_rf_period) - fit.x[3])
    plt.figure(k)
    plt.plot(S_fit)
    plt.plot(cumulative_fluorescence[:,k], '.')
    plt.show()

# -----------plot of the fit
# S_fit = fit.x[0] + fit.x[1] * np.cos(Omega * time_step * np.arange(0, len(cumulative_fluorescence[:,1])) - fit.x[2])
# plt.plot(S_fit)
# plt.plot(cumulative_fluorescence[:,1], '.')
# plt.show()
# ------------

## calculating beta
m = 40 * 1.66053904e-27 # calcium mass
k_vec = 2*np.pi / 397e-9 # wave vector
e = 1.60217662e-19 # elem charge

beta = np.zeros(N_of_data_sets)
for k in range(N_of_data_sets):
    beta[k] = get_beta_second(Omega, decay_rate, laser_detun, norm_mod_amp_second[k])

E_rf = m*Omega**2 / (k_vec*e) * beta

# ---- final plot
matlab_plot_axes = loadmat('x_axis.mat')
delta_U = matlab_plot_axes['x']

# plt.plot(delta_U, E_rf, '.',delta_U, matlab_plot_axes['y'],'.')
# plt.show()
# plt.figure(1)
# plt.plot(delta_U, E_rf, '.')
# #
# plt.figure(2)
# plt.plot(delta_U, matlab_plot_axes['y'], '.')
# plt.show()
plt.figure()
plt.plot(beta, '.')
plt.show()