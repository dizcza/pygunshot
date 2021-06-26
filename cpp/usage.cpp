#include <stdio.h>
#include <stdlib.h>
#include "gunshot.h"

#define SAMPLING_RATE 96000  // Hz


int main() {
    const float duration = 0.1f;
    const int32_t nPoints = (int32_t) (duration * SAMPLING_RATE);
    float *tInterval = (float*) malloc(sizeof(float) * nPoints);
    float *pmb = (float*) malloc(sizeof(float) * nPoints);

    GunShot::berlageMW(pmb, tInterval, nPoints, 0, 500, 5, 0.52f, 20);

    free(pmb);
    free(tInterval);

    return 0;
}
