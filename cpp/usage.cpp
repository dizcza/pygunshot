#include <stdio.h>
#include <stdlib.h>
#include "gunshot.h"

#define SAMPLING_RATE 96000  // Hz

/* Minimal code to read and display the output signal in Python
  >>> import matplotlib.pyplot as plt
  >>> import numpy as np
  >>> Pmb = np.fromfile("output.float32", dtype=np.float32)
  >>> plt.plot(Pmb)
  >>> plt.show()
 */


/* Simulate the muzzle blast wave propagation at the mic location.*/
void test_getAnechoicGunShotAtDistance() {
    const float duration = 0.1f;  // in s
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float tInterval[nPoints];
    float Pmb[nPoints];
    for (int32_t i = 0; i < nPoints; i++) {
        tInterval[i] = i * 1.f / SAMPLING_RATE;
    }

    // Values are taken from ExampleGeometry.json
    float xgun[] = {4.25, 3.0, 1.25};
    float ngun[] = {0, 1, 0};
    float xmic[] = {8.05, 22.0, 1.5};
    Gunshot::Geometry geom;
    Gunshot::initGeometry(geom, xgun, xmic, ngun);

    // RossiMagnumR971.json
    Gunshot::Gun gun(0.0091, 0.0182, 0.1016, 2817.0, 388, "Rossi Magnum R971", ".357 Magnum");

    // Call the main function
    Gunshot::getAnechoicGunshotAtDistance(Pmb, tInterval, nPoints, geom, gun);

    FILE *f = fopen("output_total.float32", "w");
    fwrite(Pmb, sizeof(float), nPoints, f);
    fclose(f);
}

/* Generate either a Friedlander or a Berlage wave at the source point (barrel
 * exit) of the specified length (duration) in s.
 */
void test_MuzzleBlastGeneration() {
    const float duration = 0.1f;  // in s
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float tInterval[nPoints];
    float Pmb[nPoints];
    for (int32_t i = 0; i < nPoints; i++) {
        tInterval[i] = i * 1.f / SAMPLING_RATE;
    }

    float ta = 0.05623617307421957;
    float tau = 0.0006669207333907665;
    float amplitude = 236.53990944414042;
    MuzzleBlast::friedlander(Pmb, tInterval, nPoints, ta, amplitude, tau);

    FILE *f = fopen("Pmb_example.float32", "w");
    fwrite(Pmb, sizeof(float), nPoints, f);
    fclose(f);
}


void test_NWaveGeneration() {
    const float duration = 0.1f;  // in s
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float tInterval[nPoints];
    float Pnw[nPoints];
    for (int32_t i = 0; i < nPoints; i++) {
        tInterval[i] = i * 1.f / SAMPLING_RATE;
    }

    // RossiMagnumR971 with the example geometry
    float pmax = 417.5844224926806;
    float ta = 0.05429687354014838;
    float theta = 0.19781125425388243;
    float Td = 0.0003323567661332675;
    float tr = 4.823159432796164e-08;
    NWave::generateNWave(Pnw, tInterval, nPoints, pmax, ta, Td, tr);

    FILE *f = fopen("Pnw_example.float32", "w");
    fwrite(Pnw, sizeof(float), nPoints, f);
    fclose(f);
}


int main() {
    test_getAnechoicGunShotAtDistance();
    test_MuzzleBlastGeneration();
    test_NWaveGeneration();
    return 0;
}
