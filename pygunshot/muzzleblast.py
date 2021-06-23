"""Muzzle blast component"""

import numpy as np

import pygunshot.util as utl

from pygunshot.domain import Gun


def friedlander(time_arr, ta, tau, Pp=1.0, pinf=101e3, csnd=341):
    """
    Calculate Friedlander wave
    
    Parameters:
    ----------------
    t -- Time (numpy array)
    ta -- Arrival time (float)
    tau -- Positive phase duration (float)
    Pp -- Peak overpressure (float)
    pinf -- atmospheric pressure at t=inf (float)

    Returns:
    ----------------
    Ps -- Muzzle blast overpressure signal (numpy array)
    """
    x = (time_arr - ta) / tau
    Ps = np.where(time_arr >= ta,
                  Pp * pinf * (1 - x) * np.exp(-x),
                  np.zeros_like(time_arr, dtype=np.float32))
    return Ps


def calculate_muzzleblast(time_arr, gun: Gun, r, theta, csnd=341., gamma=1.24):
    """
    Calculate the muzzle blast component of the gunshot sound given ballistic parameters

    Parameters
    ----------
    time_arr -- Time array in seconds
    gun -- The barrel
    mu -- Momentum index (float)
    r -- Distance to the microphone from the gun (float)
    theta -- Angle between the boreline and the microphone position in radians (float)
    gamma -- Specific heat ratio (float)

    Returns
    -------
    Pmb -- Muzzle blast signal (numpy array)
    """
    l, lp = scaling_length(gun, theta=theta, csnd=csnd, gamma=gamma)
    ta = timeOfArrival(r, lp)
    tau = positivePhaseDuration(r, lp=lp, l=l, L=gun.barrelLength,
                                Vp=gun.uexit, csnd=csnd)
    Pb = peakOverpressure(r, lp)
    Pmb = friedlander(time_arr, ta, tau, Pb)
    return Pmb


def scaling_length(gun: Gun, theta, csnd=341., gamma=1.24, pinf=101e3):
    """
    Calculate the scaling length and the direction weighted scaling length
    
    Parameters:
    ----------------
    theta -- Angle between the boreline and the microphone position in radians (float)
    gamma -- Specific heat ratio (float)

    Returns:
    ----------------
    l -- Scaling length (float)
    lp -- Weighted scaling length (float)
    """
    M = gun.mach_number(csnd)
    mu = gun.momentum_index(M, gamma=gamma)
    peb = utl.convertPressureToPascals(gun.pexit) / pinf
    # Energy deposition rate, eq. 2
    dEdt = (gamma * peb * gun.uexit) / (gamma - 1) * (
            1 + (gamma - 1) / 2 * M ** 2) * gun.bore_area
    l = np.sqrt(dEdt / (pinf * csnd))  # Eq. 3
    ratio = mu * np.cos(theta) + np.sqrt(1 - (mu * np.sin(theta)) ** 2)  # Eq.7
    lp = l * ratio
    return l, lp


def peakOverpressure(r, lp):
    """
    Calculate the peak overpressure of the muzzle blast
    
    Parameters:
    ----------------
    r -- Distance from the muzzle in m (float)
    lp -- Weighted scaling length (float)

    Returns:
    ----------------
    Pb -- Peak overpressure in Pa (float)
    """
    rb = r / lp  # Def. 8
    if rb < 50:
        Pb = 0.89 * (lp / r) + 1.61 * (lp / r) ** 2   # Eq. 25
    else:
        Pb = 3.48975 / (rb * np.sqrt(np.log(33119.0 * rb)))  # Eq. 31

    return Pb / 100


def timeOfArrival(r, lp, csnd=341):
    """
    Calculate the time of arrival of the muzzle blast
    
    Parameters:
    ----------------
    r -- Distance from the muzzle in m (float)
    lsn -- Weighted scaling length (float)

    Returns:
    ----------------
    ta --  Time of arrival in s (float)
    
    """
    rb = r / lp  # Def. 8
    X = np.sqrt(rb ** 2 + 1.04 * rb + 1.88)
    ta_norm = X - 0.52 * np.log(2 * X + 2 * rb + 1.04) - 0.56  # Eq. 27
    ta = ta_norm * lp / csnd  # Eq. 15
    return ta


def positivePhaseDuration(r, lp, l, L, Vp, csnd=341.):
    """
    positivePhaseDuration(measurementDistance, scalingLength, unscaledScalingLength, boreLengthInMeters, exitSpeedOfProjectile):
    
    Calculate the positive phase duration of the muzzle blast
    
    Parameters:
    ----------------
    r -- Distance from the muzzle in m (float)
    lp -- Weighted scaling length (float)
    l -- Scaling length (float)
    L -- Barrel length in m (float)
    Vp -- Exit speed of projectile in m/s (float)

    Returns:
    ----------------
    tau -- Positive phase duration in s (float)
    
    """
    rb = r / lp  # Def. 8
    X = np.sqrt(rb ** 2 + 1.04 * rb + 1.88)
    delta = (L * csnd) / (Vp * l)  # Blow-down parameter, eq. 10
    G = 0.09 - 0.00379 * delta + 1.07 * (
                1 - 1.36 * np.exp(-0.049 * rb)) * l / lp  # Eq. 28
    if rb < 50:
        tau_norm = rb - X + 0.52 * np.log(2 * X + 2 * rb + 1.04) + 0.56 + G
    else:
        tau_norm = 2.99 * np.sqrt(np.log(33119.0 * rb)) - 8.534 + G
    tau = tau_norm * lp / csnd  # Eq. 15, 17
    return tau
