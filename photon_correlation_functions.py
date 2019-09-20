# 19.9.2019
# file contains useful functions related to photon-correlation method

## packages
import numpy as np
from scipy.special import j0, j1
## defining a function calculating the fluorescence signal for given parameters

def get_A(decay_rate, detun):
    # returns driven oscillator amplitude in point given by detun param.
    return 1 / 2 * (decay_rate/2 - detun * 1j) / (detun ** 2 + (decay_rate / 2) ** 2)


def fl_signal(beta, laser_detun, Omega, decay_rate):
    # function calculates deltaS/S0 photon-correlation signal
    # input: beta, laser detuning, RF drive freq, decay rate
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    numer = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
    denom = j0(beta) ** 2 * np.absolute(A)**2 + j1(beta) ** 2 * (np.absolute(A_plus)**2 + np.absolute(A_minus)**2)

    return numer / denom

def signal_phase(laser_detun, Omega, decay_rate):

    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    return np.angle(np.conj(A) * A_plus - A * np.conj(A_minus))