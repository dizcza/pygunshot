import numpy as np

import pygunshot.muzzleblast as mb
import pygunshot.nwave as nw
from pygunshot.atmoshpere import atmosphericAttenuation, getReflectedSignal
from pygunshot.domain import Gun, Geometry


def getAnechoicGunShot(t_interval, geometry: Geometry, gun: Gun, Fs=96000,
                       patm=101e3, csnd=341., gamma=1.24, reflection=True):
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
    reflection -- add reflected components (MW and NW) or not

    Returns
    -------
    sig -- Total anechoic gunshot sound signal
           with atmospheric absorption and reflections
    Pmb -- Muzzle blast component
    Pnw -- N-wave component
    """
    def timeShift(signal, ta, ta_refl):
        ta_id = np.where(t_interval >= ta)[0][0]
        ta_refl_id = np.where(t_interval >= ta_refl)[0][0]
        shift_size = ta_refl_id - ta_id
        signal_shifted = np.zeros_like(signal)
        signal_shifted[ta_refl_id:] = signal[ta_id:-shift_size]
        return signal_shifted

    def getReflectedMuzzleBlast():
        Pmb_refl = np.zeros_like(Pmb)
        ta_refl, xsource = mb.timeOfArrivalReflected(geometry, csnd=csnd)
        ray = geometry.xmic - xsource
        reflected_ray_dist = geometry.xmic[2] / ray[2] * np.linalg.norm(ray)
        incidence_angle = np.arctan2(ray[2], np.linalg.norm(ray[:2]))
        if ta_refl < t_interval[-1]:
            # FIXME: MB wave is spherical but the reflected component looks
            #  more realistic with the plane wave model
            Pmb_refl = timeShift(Pmb, ta=ta_mb, ta_refl=ta_refl)
            Pmb_refl = getReflectedSignal(
                Pmb_refl,
                Fs=Fs,
                angle=incidence_angle,
                reflected_ray_dist=reflected_ray_dist,
                reflection_model='plane'
            )
        return Pmb_refl

    def getReflectedNWave():
        Pnw_refl = np.zeros_like(Pnw)
        ta_refl, xsource = nw.nWaveTimeOfArrivalReflected(geometry, gun,
                                                          csnd=csnd)
        ray = geometry.xmic - xsource
        incidence_angle = np.arctan2(ray[2], np.linalg.norm(ray[:2]))
        if ta_refl < t_interval[-1]:
            Pnw_refl = timeShift(Pnw, ta=ta_nw, ta_refl=ta_refl)
            Pnw_refl = getReflectedSignal(
                Pnw_refl,
                Fs=Fs,
                angle=incidence_angle,
                reflection_model='plane'
            )
        return Pnw_refl

    r, theta = geometry.mic_coords_polar()
    Pmb, ta_mb = mb.muzzleBlast(t_interval, gun, distance=r, theta=theta,
                                patm=patm, csnd=csnd, gamma=gamma)
    sig = Pmb.copy()
    if reflection:
        sig += getReflectedMuzzleBlast()

    Pnw, ta_nw = nw.nWave(t_interval, gun, geometry, patm=patm, csnd=csnd)
    if Pnw is not None:
        sig += Pnw
        if reflection:
            sig += getReflectedNWave()

    sig = atmosphericAttenuation(sig, distance=r, Fs=Fs)

    return sig, Pmb, Pnw
