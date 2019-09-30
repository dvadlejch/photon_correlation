# 30.9.2019
# calculation of the modulation index beta for a real measured data

# from mat4py import loadmat
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from photon_correlation_functions import get_beta
# from matplotlib import rc
## data import
data_fluorescence = loadmat('micromotion_sym_20190228.mat')['namerene_hodnoty']
len_data = len(data_fluorescence) # number of fluorescence data in each measurement of micromotion
N_of_data_sets = data_fluorescence.shape[1] # number of measurements of micromotion
## measurement parameters
Omega = 2 * np.pi * 30.171e6  # RF drive freq
decay_rate = 132e6  # units 2*pi * decay freq
# laser_detun = -1 / 2 * decay_rate  # laser detuning
laser_detun = -2*np.pi*40e6  # laser detuning

## calculating photon correlation histogram for given data
T = 2 * np.pi / Omega  # period of the RF drive
time_step = 128e-12  # time step in seconds
N_steps_per_rf_period = int(T / time_step)
N_rf_periods = int(len_data / N_steps_per_rf_period)

# t_rf_cross = T * np.arange(0, int(len_data*time_step/T))  # times of rf zero crossings
# t = time_step * np.arange(0, len_data) # time of measurement

# now I want to divide measurement period into periods of RF and sum all of the fluorescence
cumulative_fluorescence = np.zeros( (N_steps_per_rf_period, N_of_data_sets) )
for k in range(N_of_data_sets):
    for i in range(N_rf_periods):
        for j in range(N_steps_per_rf_period):
            cumulative_fluorescence[j, k] += data_fluorescence[j + i * N_steps_per_rf_period, k]



## fit of the histogram
# S(t) = S_0 + Delta S * cos(Omega*t - phi)  # Keller2015

def fit_resid(x, Omega, S, time_step):
    # x = [S_0, Delta S, phi]
    len_S = len(S)
    S_fit = x[0] + x[1] * np.cos(Omega * time_step * np.arange(0, len_S) - x[2])
    return S - S_fit

norm_mod_amp = np.zeros(N_of_data_sets)
for k in range(N_of_data_sets):
    # fit = least_squares(fit_resid, [320000, 80000, 0], args=(
    # Omega, cumulative_fluorescence[:,k], time_step), bounds=([0, 0, 0], [np.inf, np.inf, 2*np.pi]) )  # fit.x[0] = S_0, fit.x[1] = Delta S, fit.x[2] = phi
    fit = least_squares(fit_resid, [320000, 80000, 0], args=(
        Omega, cumulative_fluorescence[:, k], time_step) )  # fit.x[0] = S_0, fit.x[1] = Delta S, fit.x[2] = phi
    norm_mod_amp[k] = np.abs(fit.x[1] / fit.x[0])  # Delta S/S_0

    # S_fit = fit.x[0] + fit.x[1] * np.cos(Omega * time_step * np.arange(0, len(cumulative_fluorescence)) - fit.x[2])
    # plt.figure(k)
    # plt.plot(S_fit)
    # plt.plot(cumulative_fluorescence[:,k], '.')
    # plt.show()


## calculating beta
m = 40 * 1.66053904e-27 # calcium mass
k_vec = 2*np.pi / 397.371e-9 # wave vector
e = 1.60217662e-19 # elem charge

beta = np.zeros(N_of_data_sets)
for k in range(N_of_data_sets):
    beta[k] = get_beta(Omega, decay_rate, laser_detun, norm_mod_amp[k])

E_rf = m*Omega**2 / (k_vec*e) * beta

## ---- final plot

U_1 = loadmat('voltages_28022019.mat')['U1']
U_2 = loadmat('voltages_28022019.mat')['U2']

U1 = np.zeros(len(U_1))
U2 = np.zeros(len(U_2))
# I must change data type of U1 and U2 to float
for k in range(len(U_1)):
    U1[k] = float(U_1[k,0])
    U2[k] = float(U_2[k, 0])

#----- Erf vs z position
alpha = (U1 - U2) / (U1+U2)
a=0.000357087
b=6.14272e-05
c=0.000214573
z_pos = (a*alpha + b*alpha**3 + c*alpha**5)*1e6

# final plot
#Direct input
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 30,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

plt.figure()
plt.rcParams["figure.figsize"] = (25,15)
plt.title(r'micromotion in axial axis with SYMETRIC driving')
plt.plot(z_pos, beta, '.', markersize=10)
plt.xlabel(r'$z_{\rm{pos}} [\rm{\mu m}]$')
plt.ylabel(r'$ \beta_z $ [-]')
plt.savefig('beta_vs_z_28022019.png', format='png', dpi=150)

plt.figure()
plt.rcParams["figure.figsize"] = (25,15)
plt.title(r'SYMETRIC drive')
plt.plot(z_pos, E_rf, '.', markersize=10)
plt.xlabel(r'$z_{\rm{pos}} [\rm{\mu m}]$')
plt.ylabel(r'$ E^{\rm{rf}}_{yz} $ \rm{[V/m]}')
plt.savefig('Erf_vs_z_28022019.png', format='png', dpi=150)