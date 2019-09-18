## 18.9.2019 photon correlation signal histogram

import numpy as np
import matplotlib.pyplot as plt

# parameters
S_0 = 10000  # photon/s
Delta_S = 100  # photon/s

Omega = 2 * np.pi * 30.141e6  # RF drive freq
RF_per = 2*np.pi/Omega
Phi = 0  # phase of RF drive

## photon-correlation measurement params
N_cyc = 100  # measurement is running through N_cyc cycles of RF drive
measure_time = N_cyc * RF_per

N = 1000  # number of measurements

## fluorescence signal

# creating distribution according to fluorescence prob.
t_det = []
while len(t_det) < N:
    t = np.random.uniform(0, measure_time, 1)  # random photon-arrival times
    S = S_0 + Delta_S * np.cos(Omega * t + Phi)  # photon detection prob. non scaled

    if S > np.random.uniform(S_0 - Delta_S, S_0 + Delta_S,1):
        t_det.append(t)

RF_cross = RF_per * np.linspace(0,N_cyc)  # times when RF crosses zero

# finding time to next rf cross for each element of t
tau = np.zeros(N)
for i in range(N):
    ind = np.where(RF_cross > t_det[i])
    #print(ind[0][0])
    tau[i] = RF_cross[ind[0][0]] - t_det[i]

## plot of histogram
plt.hist(tau)
plt.show()

