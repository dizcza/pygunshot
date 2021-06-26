#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

from pygunshot.process import getAnechoicGunShot
from pygunshot.util import loadDict, recordWave
from pygunshot.domain import Gun, Geometry

# Load the geometry and ballistic data
geom = Geometry(loadDict('Geometry/ExampleGeometry.json'))
gun = Gun(loadDict('Guns/300ShortMagnum.json'))

# Set duration and sampling rate
duration = 0.1
Fs = 96000

# Calculate anechoic signal
Pmb, Pnw = getAnechoicGunShot(geom, gun, duration, Fs)

plt.plot(np.arange(len(Pmb)) / Fs, Pmb)
plt.xlim([0.05, 0.07])
plt.xlabel("Time, s")
plt.ylabel("ΔP, Pa")
# plt.savefig("Output/example.png")
plt.show()

# np.savetxt("Output/Pmb.txt", Pmb)

# Save as a normalised WAVE file
# recordWave('Output/BrowningBDA380_anechoic.wav', sig, Fs)
