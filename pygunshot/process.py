import numpy as np

import pygunshot.muzzleblast  as mb
import pygunshot.nwave as nw
from pygunshot.domain import Gun, Geometry


def getAnechoicGunShot(geometry: Geometry, gun: Gun, duration, Fs=96000,
                       csnd=341., gamma=1.24):
    """
    Get anechoic gunshot
    
    Parameters:
    ----------------
    geomDict -- Dictionary containing the scene geometry (dict)
    ballistDict -- Dictionary containing the ballistic information (dict)
    duration -- Duration of the output signal in s (float)
    Fs -- Sampling rate in Hz (int)

    Returns:
    ----------------
    sig -- Anechoic gunshot sound signal 
    """
    r, theta = geometry.mic_coords_polar()
    time_arr = np.linspace(0, duration, num=int(duration * Fs))
    Pmb = mb.calculate_muzzleblast(time_arr, gun, r, theta, csnd=csnd,
                                   gamma=gamma)

    Pnw = None
    M = gun.mach_number(csnd)
    cone_angle = gun.cone_angle(M)
    if gun.uexit > csnd and (theta < np.pi - cone_angle):  # We have sonic boom
        Pnw = nw.calculateNWave(time_arr, gun, geometry)

    return Pmb, Pnw
