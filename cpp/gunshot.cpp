/**
 * Gunshot Muzzle Blast waves generation.
 */

#include <math.h>
#include <cstring>
#include "gunshot.h"


namespace Gunshot {

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
     * @param patm  - atmospheric pressure in Pa
     * @returns the momentum index
     */
    float Gun::momentumIndex(float M, float gamma, float patm) const {
        float pe = this->pexit / patm;
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

    /* Main function.
     * It estimates the total pressure signal given a geometry, a gun, and time intervals.
     * Currently, only the muzzle blast is calculated (the shock wave will be added later).
     *
     * @param signal    - estimated total pressure signal in Pa (muzzle blast + N-Wave)
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `signal`
     * @param geometry  - the environment geometry structure (coordinates)
     * @param gun       - the gun
     * @param patm      - atmospheric pressure in Pa
     * @param csnd      - the speed of sound
     * @param gamma     - specific heat ratio for exiting propellant gas
     */
    void getAnechoicGunshotAtDistance(float *signal, const float *tInterval, int32_t nPoints, const Geometry &geometry, const Gun &gun, float patm, float csnd, float gamma) {
        float r, theta;
        micPolarCoordsFromGeometry(&r, &theta, geometry);
        MuzzleBlast::muzzleBlast(signal, tInterval, nPoints, gun, r, theta, csnd, gamma);
        if (NWave::isObserved(gun, geometry, csnd)) {
            float *Pnw = new float[nPoints];
            NWave::nWave(Pnw, tInterval, nPoints, gun, geometry, patm, csnd);
            for (int32_t i = 0; i < nPoints; i++) {
                signal[i] += Pnw[i];
            }
            delete [] Pnw;
        }
    }

    /***************** ENVIRONMENT *****************/

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

}  /* namespace Gunshot */


namespace MuzzleBlast {

    /* Calculates the muzzle blast wave at the microphone position
     *
     * @param Pmb       - estimated muzzle blast pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pmb`
     * @param gun       - the gun
     * @param r         - the distance between the gun and the microphone
     * @param theta     - the angle between the boreline and the microphone position in radians
     * @param csnd      - the speed of sound
     * @param gamma     - specific heat ratio for exiting propellant gas
     */
    void muzzleBlast(float *Pmb, const float *tInterval, int32_t nPoints, const Gunshot::Gun &gun, float r, float theta, float csnd, float gamma) {
        float l, lp;
        scalingLength(&l, &lp, gun, theta, csnd, gamma);
        float ta = timeOfArrival(r, lp, csnd);
        float tau = positivePhaseDuration(r, lp, l, gun.barrelLen, gun.velocity, csnd);
        float Pb = peakOverpressure(r, lp);
        float amplitude = Pb * ATMOSPHERIC_PRESSURE;
        friedlander(Pmb, tInterval, nPoints, ta, amplitude, tau);
    }

    /* Friedlander model of a muzzle blast wave.
     *
     * @param Pmb       - estimated muzzle blast pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param tau       - positive phase duration in s
     */
    void friedlander(float *Pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float tau) {
        float x;
        for (int32_t i = 0; i < nPoints; i++) {
            if (tInterval[i] < ta) {
                Pmb[i] = 0;
            } else {
                x = (tInterval[i] - ta) / tau;
                Pmb[i] = amplitude * (1 - x) * exp(-x);
            }
        }
    }

    /* Berlage model of a muzzle blast wave.
     *
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param nr        - the rate of rising of the front edge of the MW
     * @param alpha     - the attenuation rate of the MW
     * @param freq      - the dominant frequency of the MW
     */
    void berlage(float *Pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float nr, float alpha, float freq) {
        float t;
        for (int32_t i = 0; i < nPoints; i++) {
            if (tInterval[i] < ta) {
                Pmb[i] = 0;
            } else {
                t = tInterval[i] - ta;
                Pmb[i] = amplitude * pow(t, nr) * exp(-alpha * t) * sin(freq * t);
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
     * @param patm  - atmospheric pressure in Pa
     */
    void scalingLength(float *l, float *lp, const Gunshot::Gun &gun, float theta, float csnd, float gamma, float patm) {
        float M = gun.machNumber(csnd);
        float mu = gun.momentumIndex(M, gamma, patm);
        float peb = gun.pexit / patm;
        // Energy deposition rate, eq. 2
        float dEdT = (gamma * peb * gun.velocity) / (gamma - 1) * (1 + (gamma - 1) / 2 * M * M) * gun.boreArea;
        *l = sqrt(dEdT / (patm * csnd));
        // TODO find a better constant for the scaling length
        (*l) *= 10;
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

}  /* namespace MuzzleBlast */


namespace NWave {

    /* Generate an N-wave.
     *
     * @param Pnw       - estimated N-wave pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pnw`
     * @param pmax      - N-wave amplitude in Pa
     * @param ta        - N-wave time of arrival in s
     * @param Td        - N-wave total duration in s
     */
    void generateNWave(float *Pnw, const float *tInterval, int32_t nPoints, float pmax, float ta, float Td, float tr) {
        float t, p;
        for (int32_t i = 0; i < nPoints; i++) {
            t = tInterval[i];
            p = 0;
            if ((t > ta) & (t <= (ta + tr))) {
                p = pmax * (t - ta) / tr;
            } else if ((t > (ta + tr)) & (t <= (ta + Td - tr))) {
                p = pmax * (1 - 2 * (t - ta - tr) / (Td - 2 * tr));
            } else if ((t > (ta + Td - tr)) & (t < (ta + Td))) {
                p = pmax * ((t - (ta + Td - tr)) / tr - 1);
            }
            Pnw[i] = p;
        }
    }

    /* Calculate the N-wave component at the microphone position.
     *
     * @param Pnw       - estimated N-wave pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `signal`
     * @param gun       - the gun
     * @param geometry  - the environment geometry structure (coordinates)
     * @param patm      - atmospheric pressure in Pa
     * @param csnd      - the speed of sound
     * @returns
     *    1, if no N-wave has been observed and `Pnw` is set to zeros
     *    0, if an N-wave has been observed and `Pnw` contains correct values
     */
    int8_t nWave(float *Pnw, const float *tInterval, int32_t nPoints, const Gunshot::Gun &gun, const Gunshot::Geometry &geometry, float patm, float csnd) {
        if (!isObserved(gun, geometry, csnd)) {
            memset(Pnw, 0, nPoints * sizeof(float));
            return 1;
        }

        const float M = gun.machNumber(csnd);
        float r, theta;
        Gunshot::micPolarCoordsFromGeometry(&r, &theta, geometry);

        const float xmiss = r * sin(theta);
        const float pmax = amplitude(M, gun.bulletDiam, gun.bulletLen, xmiss, patm);
        const float ta = timeOfArrival(r, theta, gun.velocity, csnd);
        const float Td = duration(M, gun.bulletDiam, gun.bulletLen, xmiss, csnd);
        const float tr = riseTime(pmax, patm, csnd);

        generateNWave(Pnw, tInterval, nPoints, pmax, ta, Td, tr);

        return 0;
    }

    /* Check whether an N-wave is observed at the microphone position or not.
     *
     * @param gun       - the gun
     * @param geometry  - the environment geometry structure (coordinates)
     * @param csnd      - the speed of sound
     * @returns         - true if an N-wave is observed; false, otherwise
     */
    bool isObserved(const Gunshot::Gun &gun, const Gunshot::Geometry &geometry, float csnd) {
         if (gun.velocity <= csnd) {
              // supersonic speed is required
              return false;
          }
  
          float M = gun.machNumber(csnd);
          float coneAngle = gun.coneAngle(M);
          float r, theta;
          Gunshot::micPolarCoordsFromGeometry(&r, &theta, geometry);
          if (coneAngle > M_PI - theta) {
              // the observer won't see the sonic boom
              return false;
          }

          return true;
    }

    /* Calculate the N-wave rise time.
     *
     * @param pmax      - N-wave amplitude in Pa
     * @param patm      - atmospheric pressure in Pa
     * @param csnd      - the speed of sound
     * @param lamb      - air molecular mean free path in m
     * @returns tr      - the rise time in s
     */
    float riseTime(float pmax, float patm, float csnd, float lamb) {
        return (lamb / csnd) * (patm / pmax);
    }

    /* Calculate the N-wave arrival time.
     *
     * @param r         - the distance between the gun and the microphone
     * @param theta     - the angle between the boreline and the microphone position in radians
     * @param velocity  - the projectile velocity in m/s
     * @param csnd      - the speed of sound in m/s
     * @returns ta      - the time of arrival in s
     */
    float timeOfArrival(float r, float theta, float velocity, float csnd) {
        float xmiss = r * sin(theta);
        float M = velocity / csnd;
        float coneAngle = asin(1.f / M);
        float bulletTravel = r * cos(theta);
        float soundTravel = xmiss * cos(coneAngle);
        float ta = bulletTravel / velocity + soundTravel / csnd;
        return ta;
    }

    /* Calculate the N-wave total duration time.
     *
     * @param M          - the Mach number
     * @param bulletDiam - the bullet diameter in m
     * @param bulletLen  - the bullet length in m
     * @param xmiss      - the shortest dist from the mic to the bullet trajectory in m
     * @param csnd       - the speed of sound in m/s
     * @returns Td       - the total duration time in s
     */
    float duration(float M, float bulletDiam, float bulletLen, float xmiss, float csnd) {
        float L = 1.82 * bulletDiam * M * pow(xmiss, 0.25) / (pow(M * M - 1, 0.375) * pow(bulletLen, 0.25));
        float Td = L / csnd;
        return Td;
    }

    /* Calculate the N-wave amplitude.
     *
     * @param M          - the Mach number
     * @param bulletDiam - the bullet diameter in m
     * @param bulletLen  - the bullet length in m
     * @param patm       - the atmospheric pressure in Pa
     * @param csnd       - the speed of sound in m/s
     * @returns pmax     - the N-wave amplitude in Pa
     */
    float amplitude(float M, float bulletDiam, float bulletLen, float xmiss, float patm) {
        float pmax = 0.53 * patm * bulletDiam * pow(M * M - 1, 0.125) / (pow(xmiss, 0.75) * pow(bulletLen, 0.25));
        return pmax;
    }

}  /* namespace NWave */
