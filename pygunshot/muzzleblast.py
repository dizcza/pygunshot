"""Muzzle blast component"""

import numpy as np

from pygunshot.domain import Gun


def seismicPulse(t_interval, ta, c, A, k, X, omega, theta, V, T):
    """
    Seismic pulse general model [1]_.

    References
    ----------
    1. Rabinovich, E. V., Filipenko, N. Y., & Shefel, G. S. (2018, May).
       Generalized model of seismic pulse. In Journal of Physics: Conference
       Series (Vol. 1015, No. 5, p. 052025). IOP Publishing.
    """
    t = t_interval - ta
    Pmb = c * A * np.cos(k * X - omega * t + theta) / np.cosh((X - V * t) / T)
    Pmb[t_interval < ta] = 0
    return Pmb


def friedlanderMW(t_interval, ta, amplitude, tau=0.05):
    """
    Friedlander model to calculate a muzzle blast wave.

    Parameters
    ----------
    t_interval -- Time, s (numpy array)
    ta -- Arrival time, s (float)
    tau -- Positive phase duration, s (float)
    amplitude -- Amplitude, Pa (float)

    Returns
    -------
    Pmb -- estimated muzzle blast pressure along the given time intervals
    """
    x = (t_interval - ta) / tau
    Pmb = amplitude * (1 - x) * np.exp(-x)
    Pmb[t_interval < ta] = 0
    return Pmb


def berlageMW(t_interval, ta, amplitude, nr=5, alpha=0.52, freq=20, phase=0):
    """
    Berlage model to calculate a muzzle blast wave.

    Parameters
    ----------
    t_interval -- equally spaced time array in s
    ta -- the time of arrival in s
    amplitude -- the pressure amplitude in Pa
    nr -- the rate of rising of the front edge of the MW
    alpha -- the attenuation rate of the MW
    freq -- the dominant frequency of the MW

    Returns
    -------
    Pmb -- estimated muzzle blast pressure along the given time intervals
    """
    t = t_interval - ta
    Pmb = amplitude * t ** nr * np.exp(-alpha * t) * np.sin(freq * t + phase)
    Pmb[t_interval < ta] = 0
    return Pmb


def getMuzzleBlastAtDistance(t_interval, gun: Gun, r, theta, csnd=341., gamma=1.24):
    """
    Calculate the muzzle blast component of the gunshot sound given ballistic parameters

    Parameters
    ----------
    t_interval -- Time array in seconds
    gun -- The barrel
    mu -- Momentum index (float)
    r -- Distance to the microphone from the gun (float)
    theta -- Angle between the boreline and the microphone position in radians (float)
    gamma -- Specific heat ratio (float)

    Returns
    -------
    Pmb -- Muzzle blast signal (numpy array)
    """
    l, lp = scalingLength(gun, theta=theta, csnd=csnd, gamma=gamma)
    ta = timeOfArrival(r, lp)
    tau = positivePhaseDuration(r, lp=lp, l=l, barrelLen=gun.barrelLen,
                                velocity=gun.velocity, csnd=csnd)
    Pb = peakOverpressure(r, lp)
    Pmb = friedlanderMW(t_interval, ta, amplitude=Pb * 101e3, tau=tau)
    # Pmb2 = berlageMW(t_interval, ta, amplitude=Pb * 101e3)
    return Pmb


def scalingLength(gun: Gun, theta, csnd=341., gamma=1.24, pinf=101e3):
    """
    Calculate the scaling length and the direction weighted scaling length
    
    Parameters
    ----------
    theta -- Angle between the boreline and the microphone position in radians (float)
    gamma -- Specific heat ratio (float)

    Returns
    -------
    l -- Scaling length (float)
    lp -- Weighted scaling length (float)
    """
    M = gun.mach_number(csnd)
    mu = gun.momentum_index(M, gamma=gamma)
    # TODO originally was pexit/pinf
    peb = gun.pexit
    # Energy deposition rate, eq. 2
    dEdt = (gamma * peb * gun.velocity) / (gamma - 1) * (
            1 + (gamma - 1) / 2 * M ** 2) * gun.bore_area
    l = np.sqrt(dEdt / (pinf * csnd))  # Eq. 3
    ratio = mu * np.cos(theta) + np.sqrt(1 - (mu * np.sin(theta)) ** 2)  # Eq.7
    lp = l * ratio
    return l, lp


def peakOverpressure(r, lp):
    """
    Calculate the peak overpressure of the muzzle blast
    
    Parameters
    ----------
    r -- Distance from the muzzle in m (float)
    lp -- Weighted scaling length (float)

    Returns
    -------
    Pb -- Peak overpressure in Pa (float)
    """
    rb = r / lp  # Def. 8
    if rb < 50:
        Pb = 0.89 * (lp / r) + 1.61 * (lp / r) ** 2   # Eq. 25
    else:
        Pb = 3.48975 / (rb * np.sqrt(np.log(33119.0 * rb)))  # Eq. 31

    return Pb


def timeOfArrival(r, lp, csnd=341):
    """
    Calculate the time of arrival of the muzzle blast
    
    Parameters
    ----------
    r -- Distance from the muzzle in m (float)
    lp -- Weighted scaling length (float)

    Returns
    -------
    ta --  Time of arrival in s (float)
    
    """
    rb = r / lp  # Def. 8
    X = np.sqrt(rb ** 2 + 1.04 * rb + 1.88)
    ta_norm = X - 0.52 * np.log(2 * X + 2 * rb + 1.04) - 0.56  # Eq. 27
    ta = ta_norm * lp / csnd  # Eq. 15
    return ta


def positivePhaseDuration(r, lp, l, barrelLen, velocity, csnd=341.):
    """
    Calculate the positive phase duration of the muzzle blast
    
    Parameters
    ----------
    r -- Distance from the muzzle in m (float)
    lp -- Weighted scaling length (float)
    l -- Scaling length (float)
    barrelLen -- Barrel length in m (float)
    velocity -- Exit speed of projectile in m/s (float)

    Returns
    -------
    tau -- Positive phase duration in s (float)
    
    """
    rb = r / lp  # Def. 8
    X = np.sqrt(rb ** 2 + 1.04 * rb + 1.88)
    delta = (barrelLen * csnd) / (velocity * l)  # Blow-down parameter, eq. 10
    G = 0.09 - 0.00379 * delta + 1.07 * (
                1 - 1.36 * np.exp(-0.049 * rb)) * l / lp  # Eq. 28
    if rb < 50:
        tau_norm = rb - X + 0.52 * np.log(2 * X + 2 * rb + 1.04) + 0.56 + G
    else:
        tau_norm = 2.99 * np.sqrt(np.log(33119.0 * rb)) - 8.534 + G
    tau = tau_norm * lp / csnd  # Eq. 15, 17
    return tau
