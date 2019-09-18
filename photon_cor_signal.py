## 18.9.2019 photon correlation signal histogram

import numpy as np

# parameters
S_0 = 10000  # photon/s
Delta_S = 100  # photon/s

Omega = 2 * np.pi * 30.141e6  # RF drive freq
Phi = 0  # phase of RF drive

## photon-correlation measurement params
N_cyc = 100  # measurement is running through N_cyc cycles of RF drive
measure_time = N_cyc * 2 * np.pi / Omega


