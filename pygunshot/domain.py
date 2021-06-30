import math

import numpy as np
import matplotlib.pyplot as plt


class Gun:
    """
    A gun.

    bulletDiam -- Bullet diameter in m (float)
    bulletLen -- Bullet length in m (float)
    barrelLen -- Barrel length of the gun (float)
    velocity -- Exit velocity of the bullet in m/s (float)
    """
    def __init__(self, ballistDict):
        self.bulletDiam = ballistDict['bulletDiam']
        self.bulletLen = ballistDict['bulletLen']
        self.barrelLen = ballistDict['barrelLen']
        self.pexit = 98066.5 * ballistDict['pexit']
        self.velocity = ballistDict['velocity']
        self.gunlabel = ballistDict.get('gunlabel')
        self.ammolabel = ballistDict.get('ammolabel')

    def mach_number(self, csnd=341.):
        return self.velocity / csnd

    def momentum_index(self, M, gamma=1.24, patm=101e3):
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
        pe = self.pexit / patm
        xmod = M * math.sqrt(gamma * pe / 2)  # Eq. 11
        mu = 0.83 - 0.0063 * xmod  # Eq. 26
        return mu

    @staticmethod
    def cone_angle(M):
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
        """
        Returns
        -------
        r -- the dist to the mic in m
        theta -- the angle between the gun look and the mic in rads
        """
        mic_vec = self.xmic - self.xgun  # gun at the origin
        r = np.linalg.norm(mic_vec)  # distance to the mic
        theta = np.arccos(np.dot(mic_vec, self.ngun) / r)
        return r, theta

    @property
    def xmiss(self):
        """
        Returns
        -------
        xmiss -- shortest dist from the mic to the bullet trajectory in m
        """
        r, theta = self.mic_coords_polar()
        xmiss = r * np.sin(theta)
        return xmiss

    def plot_geometry(self, cone_angle=None):
        """
        Plot the gun and the mic (observer) in 2D.

        Parameters
        ----------
        cone_angle -- Mach cone angle in rads, optional
        """
        fig, ax = plt.subplots()
        ax.set_aspect('equal', 'box')
        ax.scatter(*self.xgun[:2], s=100, marker='H',
                   label=f'gun (z={self.xgun[2]} m)')
        ax.scatter(*self.xmic[:2], s=90, marker='D',
                   label=f'mic (z={self.xmic[2]} m)')
        mic_vec = self.xmic - self.xgun  # gun at the origin
        r = np.linalg.norm(mic_vec)  # distance to the mic
        ngun = 0.1 * r * self.ngun[:2]
        ax.arrow(*self.xgun[:2], *ngun, head_width=0.02 * r, color='black',
                 zorder=-1)

        # autoscale fix to mind the text positions.
        # See https://stackoverflow.com/questions/11545062/
        # matplotlib-autoscale-axes-to-include-annotations
        plt.get_current_fig_manager().canvas.draw()
        for handle in [ax.legend()]:
            bbox = handle.get_window_extent()
            bbox_data = bbox.transformed(ax.transData.inverted())
            ax.update_datalim(bbox_data.corners())
        ax.autoscale_view()

        if cone_angle is not None:
            for angle in (cone_angle, -cone_angle):
                s, c = np.sin(angle), np.cos(angle)
                mrot = np.array(((c, -s), (s, c)))
                mach_tail = mrot.dot(-ngun) + self.xgun[:2]
                ax.plot([mach_tail[0], self.xgun[0]],
                        [mach_tail[1], self.xgun[1]],
                        color='black', ls='--', lw=1)

        ax.grid()
        ax.set_xlabel('x, m')
        ax.set_ylabel('y, m')
        return ax
