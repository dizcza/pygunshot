import math

import numpy as np

import pygunshot.util as utl


class Gun:
    """
    A gun.

    bulletDiam -- Bullet diameter in m (float)
    bulletLen -- Bullet length in m (float)
    barrelLength -- Barrel length of the gun (float)
    uexit -- Exit velocty of the bullet in m/s (float)
    """
    def __init__(self, ballistDict):
        self.bulletDiam = ballistDict['bulletDiam']
        self.bulletLen = ballistDict['bulletLen']
        self.barrelLength = ballistDict['barrelLength']
        self.pexit = ballistDict['pexit']
        self.uexit = ballistDict['uexit']
        self.gunlabel = ballistDict.get('gunlabel')
        self.ammolabel = ballistDict.get('ammolabel')

    def mach_number(self, csnd=341.):
        return self.uexit / csnd

    def momentum_index(self, M, gamma=1.24):
        """
        Calculate the momentum index.

        Parameters
        ----------
        M -- Mach number (float)
        gamma -- Specific heat ratio (float)

        Returns
        -------
        mu -- Momentum index (float)
        """
        pe = utl.convertPressureToPascals(self.pexit)
        pe /= 101e3  # normalize to [0, 1] range
        xmod = M * math.sqrt(gamma * pe / 2)
        # FIXME: mu = 0.83 - 0.0063 * xmod  # Eq. 26
        mu = 0.87 - 0.01 * xmod
        return mu

    def cone_angle(self, M):
        """
        Calculate the cone angle.

        Parameters:
        ----------------
        M -- Mach number (float)

        Returns:
        ----------------
        theta -- Mach cone angle in radians (float)

        """
        return np.arcsin(1.0 / M)

    @property
    def bore_area(self):
        """
        Returns
        -------
        Ae -- Bore area in m2 (float)
        """
        return math.pi * (self.bulletDiam / 2) ** 2


class Geometry:
    """
    The environment geometry.

    xgun -- Gun position (3x1 numpy array)
    xmic -- Microphone position (3x1 numpy array)
    ngun -- Gun look direction (3x1 numpy array)
    """
    def __init__(self, geomDict):
        self.xgun = np.array(geomDict['xgun'])
        self.xmic = np.array(geomDict['xmic'])
        ngun = np.array(geomDict['ngun'])
        self.ngun = ngun / np.linalg.norm(ngun)
        self.label = geomDict.get('label')

    def mic_coords_polar(self):
        mic_direction = self.xmic - self.xgun  # gun at the origin
        r = np.linalg.norm(mic_direction)  # distance to the mic
        theta = np.arccos(np.dot(mic_direction, self.ngun) / r)
        return r, theta
