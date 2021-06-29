#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

from pygunshot.process import getAnechoicGunShot
from pygunshot.util import loadDict, recordWave
from pygunshot.domain import Gun, Geometry

# Load the geometry and ballistic data
geom = Geometry(loadDict('Geometry/ExampleGeometry.json'))
gun = Gun(loadDict('Guns/RossiMagnumR971.json'))

cone_angle = gun.cone_angle(gun.mach_number())
# geom.plot_geometry(cone_angle)
# plt.show()

# Set duration and sampling rate
duration = 0.1
Fs = 96000

# Calculate anechoic signal
t_interval = np.linspace(0, duration, num=int(duration * Fs))
sig, Pmb, Pnw = getAnechoicGunShot(t_interval, geom, gun, Fs=Fs)

plt.plot(t_interval, sig)
# plt.xlim([0.05, 0.07])
plt.xlabel("Time, s")
plt.ylabel("ΔP, Pa")
# plt.savefig("Output/example.png")
plt.show()

# np.savetxt("Output/Pmb.txt", Pmb)

# Save as a normalised WAVE file
# recordWave('Output/BrowningBDA380_anechoic.wav', sig, Fs)
