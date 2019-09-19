# 19.9.2019
# file contains useful functions related to photon-correlation method

## packages
import numpy as np
from scipy.special import j0, j1
## defining a function calculating the fluorescence signal for given parameters

def get_A(decay_rate, detun):
    # returns driven oscillator amplitude in point given by detun param.
    return 1 / 2 * (decay_rate - detun * 1j) / (detun ** 2 + (decay_rate / 2) ** 2)


def fl_signal(beta, laser_detun, Omega, decay_rate):
    A_minus = get_A(decay_rate, laser_detun - Omega)
    A_plus = get_A(decay_rate, laser_detun + Omega)
    A = get_A(decay_rate, laser_detun)

    numer = 2 * j0(beta) * j1(beta) * np.absolute(np.conj(A) * A_plus - A * np.conj(A_minus))
    denom = j0(beta) ** 2 * A * np.conj(A) + j1(beta) ** 2 * (A_plus * np.conj(A_plus) + A_minus * np.conj(A_minus))

    return numer / denom