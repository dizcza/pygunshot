import numpy as np

from pygunshot.atmoshpere import atmosphericAttenuation
from pygunshot.domain import Gun, Geometry
from pygunshot.muzzleblast import muzzleBlast
from pygunshot.nwave import nWave


def getAnechoicGunShot(t_interval, geometry: Geometry, gun: Gun, Fs=96000,
                       patm=101e3, csnd=341., gamma=1.24):
    """
    Get anechoic gunshot
    
    Parameters
    ----------
    t_interval -- Time array in seconds
    geometry -- Geometry object (gun and mic coordinates)
    gun -- The firearm
    Fs -- Sampling rate in Hz (int)
    patm -- atmospheric pressure in Pa
    csnd -- speed of sound in m/s
    gamma -- specific heat ratio

    Returns
    -------
    sig -- Total anechoic gunshot sound signal (with atmospheric absorption)
    Pmb -- Muzzle blast component
    Pnw -- N-wave component
    """
    r, theta = geometry.mic_coords_polar()
    Pmb = muzzleBlast(t_interval, gun, distance=r, theta=theta, patm=patm,
                      csnd=csnd, gamma=gamma)
    sig = Pmb.copy()

    Pnw = nWave(t_interval, gun, geometry, patm=patm, csnd=csnd)
    if Pnw is not None:
        sig += Pnw

    sig = atmosphericAttenuation(sig, distance=r, Fs=Fs)

    return sig, Pmb, Pnw
