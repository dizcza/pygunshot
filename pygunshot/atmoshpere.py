import numpy as np

from scipy.signal import fftconvolve


def sound_speed(temperature=20):
    """
    Calculate the speed of sound adjusted by the air temperature

    Parameters
    ----------
    temperature -- air temperature in Celsius

    Returns
    -------
    csnd -- speed of sound in m/s
    """
    csnd = 331.3 * np.sqrt(1 + temperature / 273.15)
    return csnd


def atmosphericAttenuation(signal, distance, Fs, **kwargs):
    """
    Apply atmospheric absorption to the `signal` by convolving it with the
    impulse response obtained from :func:`Atmosphere.impulse_response`.

    It does not account for the geometrical attenuation.

    Parameters
    ----------
    signal -- pressure signal in Pa (time domain)
    distance -- the distance to the source point in m
    Fs -- sampling frequency of the `signal` in Hz
    kwargs -- passed to `Atmosphere` class

    Returns
    -------
    signal_attenuated - attenuated signal of the same shape in the time domain
    """
    # pip install acoustics
    from acoustics.atmosphere import Atmosphere

    atm = Atmosphere(**kwargs)
    ir = atm.impulse_response(distance=distance, fs=Fs, ntaps=len(signal))
    signal_attenuated = fftconvolve(signal, ir, mode='same')
    return signal_attenuated


def getReflectedSignal(signal, Fs, angle, flow_resistivity=2e5,
                       reflected_ray_dist=None, reflection_model='plane'):
    """
    Compute the reflected signal from the incidence `signal`.

    Parameters
    ----------
    signal -- pressure signal in Pa (time domain)
    Fs -- sampling frequency of the `signal` in Hz
    angle -- angle of incidence in rads
    flow_resistivity -- boundary flow resistivity. Defaults to grass
    reflected_ray_dist -- the distance travelled by the reflected ray
                          (needs only for the spherical reflection model)
    reflection_model -- either 'plane' or 'spherical'

    Returns
    -------
    signal_reflected -- reflected signal of the same shape in the time domain
    """
    from acoustics.signal import impulse_response_real_even
    from acoustics.reflection import Boundary

    freq = np.fft.rfftfreq(len(signal), d=1. / Fs)
    freq[0] = 1e-10  # avoid zero
    b = Boundary(freq, flow_resistivity=flow_resistivity, angle=angle,
                 distance=reflected_ray_dist,
                 reflection_model=reflection_model)
    ir = impulse_response_real_even(b.reflection_factor, ntaps=len(signal))
    signal_reflected = fftconvolve(signal, ir, mode='same')
    return signal_reflected


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from pygunshot.muzzleblast import friedlanderMW

    Fs = 96000
    duration = 0.5   # s
    distance = 1000  # m
    t_interval = np.linspace(0, duration, int(Fs * duration))
    pmb = friedlanderMW(t_interval, ta=0.1, amplitude=300, tau=0.05)
    plt.plot(t_interval, pmb, ls='--', label='Friedlander')
    pmb = atmosphericAttenuation(pmb, distance=distance, Fs=Fs)
    plt.plot(t_interval, pmb, label='Friedlander with atm absorption')
    plt.title(f"The influence of atmospheric absorption at {distance} m")
    plt.xlabel("Time, s")
    plt.ylabel("ΔP, Pa")
    plt.legend()
    plt.show()
