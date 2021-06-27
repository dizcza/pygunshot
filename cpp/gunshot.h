/**
 * Gunshot Muzzle Blast waves generation.
 */

#ifndef GUNSHOT_H
#define GUNSHOT_H

extern "C" {
#include <inttypes.h>
}

#include <string>

#define SOUND_SPEED           341
#define GAMMA                 1.24f
#define ATMOSPHERIC_PRESSURE  101000

namespace GunShot {
    class Gun {
    public:
        float bulletDiam;  // Bullet diameter in m
        float bulletLen;   // Bullet length in m
        float barrelLen;   // Barrel length in m
        float pexit;       // Muzzle pressure in Pa
        float velocity;    // Exit velocity of the bullet in m/s
        float boreArea;    // Bore area in m2
        std::string gunLabel;
        std::string ammoLabel;

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
        Gun(float bulletDiam, float bulletLen, float barrelLen, float pexit, float velocity, std::string gunLabel = NULL, std::string ammoLabel = NULL);

        /* Calculates the Mach number.
         *
         * @param csnd - the speed of sound
         * @returns the Mach number at projectile ejection
         */
        float machNumber(float csnd = SOUND_SPEED) const;

        /* Calculates the momentum index defined as the ratio of the sound
         * pressure at the front and at the rear of the firing position.
         *
         * @param M     - the Mach number
         * @param gamma - specific heat ratio for exiting propellant gas
         * @param pinf  - atmospheric pressure at t=inf
         * @returns the momentum index
         */
        float momentumIndex(float M, float gamma = GAMMA, float pinf = ATMOSPHERIC_PRESSURE) const;

        /* Calculates the Mach cone angle.
         *
         * @param M - the Mach number
         * @returns the cone angle in radians
         */
        float coneAngle(float M) const;
    };

    typedef struct Geometry {
        float xgun[3];  // Gun coordinates
        float xmic[3];  // Microphone coordinates
        float ngun[3];  // Gun look direction (a unit vector)
    } Geometry;

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
    void getAnechoicGunShotAtDistance(float *pmb, const float *tInterval, int32_t nPoints, const Geometry &geometry, const Gun &gun, float csnd = SOUND_SPEED, float gamma = GAMMA);

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
    void getMuzzleBlastAtDistance(float *pmb, const float *tInterval, int32_t nPoints, const Gun &gun, float r, float theta, float csnd = SOUND_SPEED, float gamma = GAMMA);

    /* Friedlander model of a muzzle blast wave.
     *
     * @param pmb       - estimated muzzle blast pressure (a vec of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param tau       - positive phase duration in s
     */
    void friedlanderMW(float *pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float tau = 0.05f);

    /* Berlage model of a muzzle blast wave.
     *
     * @param pmb       - estimated muzzle blast pressure (a vec of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param nr        - the rate of rising of the front edge of the MW
     * @param alpha     - the attenuation rate of the MW
     * @param freq      - the dominant frequency of the MW
     */
    void berlageMW(float *pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float nr = 5, float alpha = 0.52f, float freq = 20);

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
    void scalingLength(float *l, float *lp, const Gun &gun, float theta, float csnd = SOUND_SPEED, float gamma = GAMMA, float pinf = ATMOSPHERIC_PRESSURE);

    /* Calculates the peak overpressure of the muzzle blast.
     *
     * @param r  - the distance between the gun and the microphone
     * @param lp - the direction weighted scaling length in m
     * @returns the peak overpressure in Pa
     */
    float peakOverpressure(float r, float lp);

    /* Calculates the time of arrival of the muzzle blast.
     *
     * @param r    - the distance between the gun and the microphone
     * @param lp   - the direction weighted scaling length in m
     * @param csnd - the speed of sound
     * @returns the time of arrival in s
     */
    float timeOfArrival(float r, float lp, float csnd = SOUND_SPEED);

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
    float positivePhaseDuration(float r, float lp, float l, float barrelLen, float velocity, float csnd = SOUND_SPEED);

    /***************** ENVIRONMENT *****************/

    /* Initialize the Geometry structure.
     *
     * @param geometry - Geometry object to be initialized
     * @param xgun - gun 3D coordinates
     * @param xmic - mic 3D coordinates
     * @param ngun - gun barren look 3D direction
     */
    void initGeometry(Geometry &geometry, const float *xgun, const float *xmic, const float *ngun);

    /* Get polar coordinates (r and theta) of the mic w.r.t. the barrel.
     *
     * @param r        - (to be written) the distance between the gun and the microphone
     * @param theta    - (to be written) the angle between the boreline and the microphone position in radians
     * @param geometry - the environment geometry structure (coordinates)
     */
    void micPolarCoordsFromGeometry(float *r, float *theta, const Geometry &geometry);

    /* Calculates the speed of sound adjusted by the air temperature.
     *
     * @param temperature - the air temperature in Celsius
     * @returns the speed of sound in m/s
     */
    float getSoundSpeed(float temperature);

}  /* namespace GunShot */

#endif /* GUNSHOT_H */
