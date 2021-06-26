import json
from pathlib import Path

import scipy.io.wavfile as wav


def recordWave(filename, sig, Fs=192000):
    """
    Write the normalized signal as a wav file
    
    Parameters:
    ----------------
    filename -- File name including its path (string)
    sig -- Signal to be stored as a wav file (numpy array)
    Fs -- Sampling rate in Hz (float)    
    """
    Path(filename).parent.mkdir(exist_ok=True)
    siga = abs(sig)
    sig = sig / siga.max()
    wav.write(filename, Fs, sig)


def loadDict(filename):
    """
    Load a dictionary stored a JSON object
    
    Parameters:
    ----------------
    dictionary -- Dictionary of items (dict)
    filename -- Path and name of the file to be written (string)    
    """
    with open(filename, 'r') as f:
        return json.load(f)
