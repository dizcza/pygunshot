/**
 * Gunshot Muzzle Blast waves generation.
 */

#include <math.h>
#include "gunshot.h"

namespace GunShot {

    /* Gun (barrel) constructor
     *
     * @param bulletDiam - the bullet diameter in m
     * @param bulletLen  - the bullet length in m
     * @param barrelLen  - the barrel length in m
     * @param pexit      - the muzzle pressure in kg/cm2
     * @param velocity   - the exit velocity of the bullet in m/s
     * @param gunLabel   - the gun label, optional
     * @param ammoLabel  - the ammo label, optional
     */
    Gun::Gun(float bulletDiam, float bulletLen, float barrelLen, float pexit, float velocity, std::string gunLabel, std::string ammoLabel) {
        this->bulletDiam = bulletDiam;
        this->bulletLen = bulletLen;
        this->barrelLen = barrelLen;
        this->pexit = 98066.5f * pexit;  // convert to Pa
        this->velocity = velocity;
        this->boreArea = M_PI * bulletDiam * bulletDiam / 4;
        this->gunLabel = gunLabel;
        this->ammoLabel = ammoLabel;
    }

    /* Calculates the Mach number.
     *
     * @param csnd - the speed of sound
     * @returns the Mach number at projectile ejection
     */
    float Gun::machNumber(float csnd) const {
        return this->velocity / csnd;
    }

    /* Calculates the momentum index defined as the ratio of the sound
     * pressure at the front and at the rear of the firing position.
     *
     * @param M     - the Mach number
     * @param gamma - specific heat ratio for exiting propellant gas
     * @param pinf  - atmospheric pressure at t=inf
     * @returns the momentum index
     */
    float Gun::momentumIndex(float M, float gamma, float pinf) const {
        float pe = this->pexit / pinf;
        float xmod = M * sqrt(gamma * pe / 2);
        float mu = 0.83f - 0.0063f * xmod;
        return mu;
    }

    /* Calculates the Mach cone angle.
     *
     * @param M - the Mach number
     * @returns the cone angle in radians
     */
    float Gun::coneAngle(float M) const {
        return asin(1.f / M);
    }

    /* Initialize the Geometry structure.
     *
     * @param geometry - Geometry object to be initialized
     * @param xgun - gun 3D coordinates
     * @param xmic - mic 3D coordinates
     * @param ngun - gun barren look 3D direction
     */
    void initGeometry(Geometry &geometry, const float *xgun, const float *xmic, const float *ngun) {
        float ngun_norm = ngun[0] * ngun[0] + ngun[1] * ngun[1] + ngun[2] * ngun[2];
        ngun_norm = sqrt(ngun_norm);
        for (int i = 0; i < 3; i++) {
            geometry.xgun[i] = xgun[i];
            geometry.xmic[i] = xmic[i];
            geometry.ngun[i] = ngun[i] / ngun_norm;
        }
    }

    /* Main function.
     * It estimates the total pressure signal given a geometry, a gun, and time intervals.
     * Currently, only the muzzle blast is calculated (the shock wave will be added later).
     *
     * @param pmb       - estimated muzzle blast pressure (a vec of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param geometry  - the environment geometry structure (coordinates)
     * @param gun       - the gun
     * @param csnd      - the speed of sound
     * @param gamma     - specific heat ratio for exiting propellant gas
     */
    void getAnechoicGunShotAtDistance(float *pmb, const float *tInterval, int32_t nPoints, const Geometry &geometry, const Gun &gun, float csnd, float gamma) {
        float r, theta;
        micPolarCoordsFromGeometry(&r, &theta, geometry);
        getMuzzleBlastAtDistance(pmb, tInterval, nPoints, gun, r, theta, csnd, gamma);
        // TODO: add shock wave
    }

    /***************** MUZZLE BLAST *****************/

    /* Calculates the muzzle blast wave at the microphone position
     *
     * @param pmb       - estimated muzzle blast pressure (a vec of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param gun       - the gun
     * @param r         - the distance between the gun and the microphone
     * @param theta     - the angle between the boreline and the microphone position in radians
     * @param csnd      - the speed of sound
     * @param gamma     - specific heat ratio for exiting propellant gas
     */
    void getMuzzleBlastAtDistance(float *pmb, const float *tInterval, int32_t nPoints, const Gun &gun, float r, float theta, float csnd, float gamma) {
        float l, lp;
        scalingLength(&l, &lp, gun, theta, csnd, gamma);
        float ta = timeOfArrival(r, lp, csnd);
        float tau = positivePhaseDuration(r, lp, l, gun.barrelLen, gun.velocity, csnd);
        float Pb = peakOverpressure(r, lp);
        float amplitude = Pb * ATMOSPHERIC_PRESSURE;
        friedlanderMW(pmb, tInterval, nPoints, ta, amplitude, tau);
        //berlageMW(pmb, tInterval, nPoints, ta, amplitude);
    }

    /* Friedlander model of a muzzle blast wave.
     *
     * @param pmb       - estimated muzzle blast pressure (a vec of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param tau       - positive phase duration in s
     */
    void friedlanderMW(float *pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float tau) {
        float x;
        for (int32_t i = 0; i < nPoints; i++) {
            if (tInterval[i] < ta) {
                pmb[i] = 0;
            } else {
                x = (tInterval[i] - ta) / tau;
                pmb[i] = amplitude * (1 - x) * exp(-x);
            }
        }
    }

    /* Berlage model of a muzzle blast wave.
     *
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param nr        - the rate of rising of the front edge of the MW
     * @param alpha     - the attenuation rate of the MW
     * @param freq      - the dominant frequency of the MW
     */
    void berlageMW(float *pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float nr, float alpha, float freq) {
        float t;
        for (int32_t i = 0; i < nPoints; i++) {
            if (tInterval[i] < ta) {
                pmb[i] = 0;
            } else {
                t = tInterval[i] - ta;
                pmb[i] = amplitude * pow(t, nr) * exp(-alpha * t) * sin(freq * t);
            }
        }
    }

    /* Calculates the scaling length `l` and the direction weighted scaling length `lp`.
     *
     * @param l     - (to be written) the scaling length in m
     * @param lp    - (to be written) the direction weighted scaling length in m
     * @param gun   - the gun
     * @param theta - the angle between the boreline and the microphone position in radians
     * @param csnd  - the speed of sound
     * @param gamma - specific heat ratio for exiting propellant gas
     * @param pinf  - atmospheric pressure at t=inf
     */
    void scalingLength(float *l, float *lp, const Gun &gun, float theta, float csnd, float gamma, float pinf) {
        float M = gun.machNumber(csnd);
        float mu = gun.momentumIndex(M, gamma, pinf);
        float peb = gun.pexit / pinf;
        // Energy deposition rate, eq. 2
        float dEdT = (gamma * peb * gun.velocity) / (gamma - 1) * (1 + (gamma - 1) / 2 * M * M) * gun.boreArea;
        *l = sqrt(dEdT / (pinf * csnd));
        float ratio = mu * cos(theta) + sqrt(1 - pow(mu * sin(theta), 2));
        *lp = (*l) * ratio;
    }

    /* Calculates the peak overpressure of the muzzle blast.
     *
     * @param r  - the distance between the gun and the microphone
     * @param lp - the direction weighted scaling length in m
     * @returns the peak overpressure in Pa
     */
    float peakOverpressure(float r, float lp) {
        const float rb = r / lp;
        const float rb_inv = lp / r;
        float Pb;
        if (rb < 50) {
            Pb = 0.89f * rb_inv + 1.61f * rb_inv * rb_inv;
        } else {
            Pb = 3.48975f / (rb * sqrt(log(33119 * rb)));
        }
        return Pb;
    }

    /* Calculates the time of arrival of the muzzle blast.
     *
     * @param r    - the distance between the gun and the microphone
     * @param lp   - the direction weighted scaling length in m
     * @param csnd - the speed of sound
     * @returns the time of arrival in s
     */
    float timeOfArrival(float r, float lp, float csnd) {
        const float rb = r / lp;
        const float X = sqrt(rb * rb + 1.04f * rb + 1.88f);
        float ta_norm = X - 0.52f * log(2 * X + 2 * rb + 1.04f) - 0.56f;  // Eq. 27
        float ta = ta_norm * lp / csnd;
        return ta;
    }

    /* Calculates the positive phase duration of the muzzle blast.
     *
     * @param r         - the distance between the gun and the microphone
     * @param lp        - the direction weighted scaling length in m
     * @param l         - the scaling length in m
     * @param barrelLen - the barrel length
     * @param velocity  - the exit velocity of the bullet in m/s
     * @param csnd      - the speed of sound
     * @returns the positive phase duration in s
     */
    float positivePhaseDuration(float r, float lp, float l, float barrelLen, float velocity, float csnd) {
        const float rb = r / lp;
        const float X = sqrt(rb * rb + 1.04f * rb + 1.88f);
        const float delta = (barrelLen * csnd) / (velocity * l);
        const float G = 0.09 - 0.00379 * delta + 1.07 * (1 - 1.36 * exp(-0.049 * rb)) * l / lp;  // Eq. 28
        float tau_norm;
        if (rb < 50) {
            tau_norm = rb - X + 0.52 * log(2 * X + 2 * rb + 1.04) + 0.56 + G;
        } else {
            tau_norm = 2.99 * sqrt(log(33119 * rb)) - 8.534 + G;
        }
        float tau = tau_norm * lp / csnd;
        return tau;
    }

    /***************** ENVIRONMENT *****************/

    /* Get polar coordinates (r and theta) of the mic w.r.t. the barrel.
     *
     * @param r        - (to be written) the distance between the gun and the microphone
     * @param theta    - (to be written) the angle between the boreline and the microphone position in radians
     * @param geometry - the environment geometry structure (coordinates)
     */
    void micPolarCoordsFromGeometry(float *r, float *theta, const Geometry &geometry) {
        float cos_anlge = 0;
        float dist = 0;
        float dx;
        for (int i = 0; i < 3; i++) {
            dx = geometry.xmic[i] - geometry.xgun[i];
            dist += dx * dx;
            cos_anlge += dx * geometry.ngun[i];
        }
        dist = sqrt(dist);
        cos_anlge /= dist;

        *r = dist;
        *theta = acos(cos_anlge);
    }

    /* Calculates the speed of sound adjusted by the air temperature.
     *
     * @param temperature - the air temperature in Celsius
     * @returns the speed of sound in m/s
     */
    float getSoundSpeed(float temperature) {
        return 331.3 * sqrt(1 + temperature / 273.15);
    }

}  /* namespace GunShot */
