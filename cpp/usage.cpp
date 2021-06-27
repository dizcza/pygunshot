#include <stdio.h>
#include <stdlib.h>
#include "gunshot.h"

#define SAMPLING_RATE 96000  // Hz

/* Minimal code to read and display the output signal in Python
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> pmb = np.fromfile("output.float32", dtype=np.float32)
>>> plt.plot(pmb)
>>> plt.show()
 */


/* Simulate the muzzle blast wave propagation at the mic location.*/
void test_getAnechoicGunShotAtDistance() {
    const float duration = 1.5f;  // in s
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float tInterval[nPoints];
    float pmb[nPoints];
    for (int i = 0; i < nPoints; i++) {
        tInterval[i] = i * 1.f / nPoints;
    }

    // Values are taken from ExampleGeometry.json
    float xgun[] = {4.25, 3.0, 1.25};
    float ngun[] = {0, 1, 0};
    float xmic[] = {8.05, 22.0, 1.5};
    GunShot::Geometry geom;
    GunShot::initGeometry(geom, xgun, xmic, ngun);

    // 300ShortMagnum.json
    GunShot::Gun gun(0.0091, 0.0182, 0.1016, 2646.8, 221, "300 Short Magnum", ".357 Magnum");

    // Call the main function
    getAnechoicGunShotAtDistance(pmb, tInterval, nPoints, geom, gun);

    // Read in Python with `pmb = np.fromfile(f, dtype=np.float32)`
    FILE *f = fopen("output_at_distance.float32", "w");
    fwrite(pmb, sizeof(float), nPoints, f);
    fclose(f);
}

/* Generate either a Friedlander or a Berlage wave at the source point (barrel
 * exit) of the specified length (duration) in s.
 */
void test_MuzzleBlastGeneration() {
    const float duration = 1.5f;  // in s
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float tInterval[nPoints];
    float pmb[nPoints];
    for (int i = 0; i < nPoints; i++) {
        tInterval[i] = i * 1.f / nPoints;
    }

//    GunShot::berlageMW(pmb, tInterval, nPoints, 0, 300, 5, 0.52f, 20);
    GunShot::friedlanderMW(pmb, tInterval, nPoints, 0, 300, 0.05f);

    // Read in Python with `pmb = np.fromfile(f, dtype=np.float32)`
    FILE *f = fopen("output_at_source_point.float32", "w");
    fwrite(pmb, sizeof(float), nPoints, f);
    fclose(f);
}


int main() {
    test_getAnechoicGunShotAtDistance();
    test_MuzzleBlastGeneration();
    return 0;
}
