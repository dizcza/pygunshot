"""N-wave component"""

import numpy as np

from pygunshot.domain import Gun, Geometry


def nWaveAmplitude(M, bulletDiam, bulletLen, xmiss, patm=101e3):
    """
    Calculate N-wave overpressure amplitude

    Parameters
    ----------
    M -- Mach number
    bulletDiam -- bullet diameter in m
    bulletLen -- bullet length in m
    xmiss -- miss distance in m
    patm -- atmosperic pressure in Pa

    Returns
    -------
    pmax -- N-wave overpressure amplitude in Pa
    """
    pmax = 0.53 * patm * bulletDiam * ((M ** 2 - 1) ** 0.125) / (
            (xmiss ** 0.75) * (bulletLen ** 0.25))
    return pmax


def nWaveDuration(M, bulletDiam, bulletLen, xmiss, csnd=341):
    """
    Calculate N-wave period
    
    Parameters
    ----------
    M -- Mach number
    bulletDiam -- bullet diameter in m
    bulletLen -- bullet length in m
    xmiss -- miss distance in m
    csnd -- speed of sound in m/s

    Returns
    -------
    Td -- N-wave period in s
    """
    L = 1.82 * bulletDiam * M * (xmiss ** 0.25) / (
            ((M ** 2 - 1) ** 0.375) * (bulletLen ** 0.25))
    Td = L / csnd
    return Td


def nWaveTimeOfArrival(r, theta, cone_angle, velocity, csnd=341):
    """
    Calculate N-wave time of arrival

    Parameters
    ----------
    r -- the dist to the mic in m
    theta -- the angle between the gun look and the mic in rads
    cone_angle -- cone Mach angle in rads
    velocity -- projectile velocity in m/s
    csnd -- speed of sound in m/s

    Returns
    -------
    ta -- N-wave time of arrival in s
    """
    xmiss = r * np.sin(theta)
    bullet_travel = r * np.cos(theta)
    sound_travel = xmiss * np.cos(cone_angle)
    ta = bullet_travel / velocity + sound_travel / csnd
    return ta


def nWaveTimeOfArrivalReflected(geometry: Geometry, gun: Gun, csnd=341):
    """
    Calculate the time of arrival of the reflected N-wave

    Parameters
    ----------
    geometry -- Geometry object
    gun -- the firearm
    csnd -- the speed of sound in m/s

    Returns
    -------
    ta -- time of arrival in s
    xsource -- 3D coordinates of the imaginary source where the ray is started
    """
    gravity = np.array([0, 0, 1])
    r, theta = geometry.mic_coords_polar()
    xmiss = r * np.sin(theta)
    cone_angle = gun.cone_angle(gun.mach_number(csnd))
    bullet_travel = r * np.cos(theta) - xmiss * np.tan(cone_angle)
    xsource = geometry.xgun + bullet_travel * geometry.ngun
    xsource -= 2 * xsource[2] * gravity
    sound_travel = np.linalg.norm(xsource - geometry.xmic)
    ta = bullet_travel / gun.velocity + sound_travel / csnd
    return ta, xsource


def nWaveRiseTime(pmax, patm=101e3, csnd=341, lamb=6.8e-8):
    """
    Calculate N-wave rise time

    Parameters
    ----------
    pmax -- N-wave overpressure amplitude in Pa
    patm -- atmospheric pressure in Pa
    csnd -- speed of sound in m/s
    lamb -- air molecular mean free path

    Returns
    -------
    trise -- N-wave rise time in s
    """
    trise = (lamb / csnd) * (patm / pmax)
    return trise
    

def nWave(t_interval, gun: Gun, geometry: Geometry, patm=101e3, csnd=341):
    """
    Calculate the N-wave (shock wave) at the microphone position.

    If the projectile velocity speed is smaller than the speed of sound or the
    observer position is in the initial Mach cone, no sonic boom will be
    detected, and a None is returned.

    Parameters
    ----------
    xgun -- Gun position (3x1 numpy array)
    ngun -- Barrel look direction (3x1 numpy array)
    xmic -- Microphone position (3x1 numpy array)
    v -- Projectile velocity in m/s (float)
    d -- Projectile diameter im m (float)
    l -- Projectile length in m (float)
    x -- Miss distance in m (float)

    Returns
    -------
    Pnw -- N-wave pressure signal (numpy array)
    ta -- time of arrival in s
    """
    if gun.velocity <= csnd:
        # supersonic speed is required
        return None, None
    M = gun.mach_number(csnd)
    cone_angle = gun.cone_angle(M)
    r, theta = geometry.mic_coords_polar()
    if cone_angle > np.pi - theta:
        # the observer won't see the sonic boom
        return None, None

    t = t_interval
    Pnw = np.zeros_like(t)
    xmiss = r * np.sin(theta)  # shortest dist to the bullet trajectory
    pmax = nWaveAmplitude(M, bulletDiam=gun.bulletDiam,
                          bulletLen=gun.bulletLen, xmiss=xmiss, patm=patm)
    ta = nWaveTimeOfArrival(r, theta, cone_angle, gun.velocity, csnd=csnd)
    Td = nWaveDuration(M, bulletDiam=gun.bulletDiam, bulletLen=gun.bulletLen,
                       xmiss=xmiss, csnd=csnd)
    tr = nWaveRiseTime(pmax, patm=patm, csnd=csnd)

    mask_rise1 = (t > ta) & (t <= ta + tr)
    mask_fall = (t > ta + tr) & (t <= ta + Td - tr)
    mask_rise2 = (t > ta + Td - tr) & (t < ta + Td)
    Pnw[mask_rise1] = pmax * (t[mask_rise1] - ta) / tr
    Pnw[mask_fall] = pmax * (1 - 2 * (t[mask_fall] - (ta + tr)) / (
            Td - 2 * tr))
    Pnw[mask_rise2] = pmax * ((t[mask_rise2] - (ta + Td - tr)) / tr - 1)
    return Pnw, ta
