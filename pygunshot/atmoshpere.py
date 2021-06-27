import numpy as np


def atmosphericAttenuation(signal, distance, Fs, **kwargs):
    """
    Apply atmospheric absorption to the `signal` for all its FFT frequencies.
    It does not account for the geometrical attenuation.

    Parameters
    ----------
    signal - a pressure waveform (time domain)
    distance - the distance to the source point, m
    Fs - sampling frequency of the `signal`, Hz
    kwargs - passed to `Atmosphere` class

    Returns
    -------
    signal_attenuated - attenuated signal in the original time domain
    """
    # pip install acoustics
    from acoustics.atmosphere import Atmosphere

    atm = Atmosphere(**kwargs)
    signal_rfft = np.fft.rfft(signal)
    freq = np.fft.rfftfreq(n=len(signal), d=1. / Fs)
    a_coef = atm.attenuation_coefficient(freq)
    # signal_rfft *= 10 ** (-a_coef * distance / 20)
    signal_rfft *= np.exp(-a_coef * distance)
    signal_attenuated = np.fft.irfft(signal_rfft)
    return signal_attenuated
