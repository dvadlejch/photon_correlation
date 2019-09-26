## 18.9.2019 photon correlation signal histogram

import numpy as np
import matplotlib.pyplot as plt

# parameters
S_0 = 1400  # photon/s
Delta_S = 1000  # photon/s

Omega = 2 * np.pi * 30.151e6  # RF drive freq
RF_per = 2*np.pi/Omega
Phi = 0  # phase of RF drive

## photon-correlation measurement params
N_cyc = 200  # measurement is running through N_cyc cycles of RF drive
measure_time = N_cyc * RF_per

N = 60000  # number of measurements

## fluorescence signal

# creating distribution according to fluorescence prob.
t_det = np.zeros(N)
k = 0
while k < N:
    t = np.random.uniform(0, measure_time, 1)  # random photon-arrival times
    S = S_0 + Delta_S * np.cos(Omega * t + Phi)  # photon detection prob. non scaled

    if S > np.random.uniform(0, S_0 + Delta_S,1):
        t_det[k] = t
        k += 1

RF_cross = RF_per * np.arange(N_cyc+1)  # times when RF crosses zero
#print(RF_cross)
# finding time to next rf cross for each element of t
tau = np.zeros(N)
for i in range(N):
    ind = np.where(RF_cross > t_det[i])
    #print(ind)
    tau[i] = RF_cross[ind[0][0]] - t_det[i]
    #print(RF_cross[ind[0][0]])

## plot of histogram
plt.hist(np.array(t_det)/RF_per, bins=int(N/500)) # time histogram
plt.hist(tau/RF_per, bins=int(N/500))
plt.show()


