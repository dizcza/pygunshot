def atmosphericAttenuation(signal, distance, Fs, **kwargs):
    """
    Apply atmospheric absorption to the `signal` by convolving it with the
    impulse response obtained from :func:`Atmosphere.impulse_response`.

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
    ir = atm.impulse_response(distance=distance, fs=Fs, ntaps=len(signal))
    signal_attenuated = np.convolve(signal, ir, mode='same')
    return signal_attenuated


if __name__ == '__main__':
    import numpy as np
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
