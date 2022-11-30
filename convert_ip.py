from astropy.modeling import models
from astropy import units as u
from astropy import constants as const
import numpy as np 
from scipy.integrate import simps
HEV = 4.13620e-15	# Planck's constant in eV 
C  = 2.997925e10
H = 6.6262e-27

def get_q_over_l(freq, spectrum):
    '''
    Calculate the ratio of ionizing photons to ionizing luminosity (up to 1000 Ryd)
    for a given spectrum

    Parameters:
        freq        array-like
                    frequency in Hz for spectrum (used for integrals)

        spectrum    array-like
                    spectrum shape in L_nu, F_nu or B_nu units. 

    Returns:
        Q/L         array-like
                    ratio of ionizing photons to ionizing luminosity (up to 1000 Ryd)
                    same shape as input
    '''
    rydberg_Hz = 13.6 / HEV 
    thousand_eV_Hz = 1000.0 / HEV
    select_Q = (freq>=rydberg_Hz)
    select_L = (freq >= rydberg_Hz) * (freq <= thousand_eV_Hz)
    
    # this is the integrand 
    q_integrand = spectrum / H / freq

    # integrate only over the selected ranged 
    Q = simps(q_integrand[select_Q], x=freq[select_Q])
    L = simps(spectrum[select_L], x=freq[select_L])
    
    return (Q/L)

def get_conversion_factor_blackbody(temp_eV = 50.0):
    '''
    Calculate the ratio of the two different ionization parameters (U/xi) 
    for a blackbody spectrum


    Parameters:
        temp_eV     float
                    temperature in eV

    Returns:
        factor      float
                    ratio of the two different ionization parameters
    '''
    bb = models.BlackBody(temperature=50.0*u.eV/const.k_B)
    freq = np.logspace(13,20,10000)
    spec = bb(freq)

    Q_over_L = get_q_over_l(freq, spec)
    factor = 1.0 / 4.0 / np.pi / C * Q_over_L
    return factor