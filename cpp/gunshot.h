/**
 * Gunshot Muzzle Blast waves generation.
 */

#ifndef GUNSHOT_H
#define GUNSHOT_H

extern "C" {
#include <inttypes.h>
}

#include <string>

// Default parameters
#define SOUND_SPEED                   341
#define GAMMA                         1.24f
#define ATMOSPHERIC_PRESSURE          101000
#define AIR_MOLECULAR_MEAN_FREE_PATH  6.8e-8


namespace Gunshot {
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
         * @param patm  - atmospheric pressure in Pa
         * @returns the momentum index
         */
        float momentumIndex(float M, float gamma = GAMMA, float patm = ATMOSPHERIC_PRESSURE) const;

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
     * It estimates the total pressure signal (muzzle blast + N-wave) given
     * a geometry, a gun, and time intervals.
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
    void getAnechoicGunshotAtDistance(float *signal, const float *tInterval, int32_t nPoints, const Geometry &geometry, const Gun &gun, float patm = ATMOSPHERIC_PRESSURE, float csnd = SOUND_SPEED, float gamma = GAMMA);

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
    void muzzleBlast(float *Pmb, const float *tInterval, int32_t nPoints, const Gunshot::Gun &gun, float r, float theta, float csnd = SOUND_SPEED, float gamma = GAMMA);

    /* Friedlander model of a muzzle blast wave.
     *
     * @param Pmb       - estimated muzzle blast pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param tau       - positive phase duration in s
     */
    void friedlander(float *Pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float tau);

    /* Berlage model of a muzzle blast wave.
     *
     * @param Pmb       - estimated muzzle blast pressure (array of size nPoints) along the given time intervals
     * @param tInterval - equally spaced time array in s
     * @param nPoints   - the length of `tInterval` and `Pmb`
     * @param ta        - the time of arrival in s
     * @param amplitude - the pressure amplitude in Pa
     * @param nr        - the rate of rising of the front edge of the MW
     * @param alpha     - the attenuation rate of the MW
     * @param freq      - the dominant frequency of the MW
     */
    void berlage(float *Pmb, const float *tInterval, int32_t nPoints, float ta, float amplitude, float nr = 5, float alpha = 0.52f, float freq = 20);

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
    void scalingLength(float *l, float *lp, const Gunshot::Gun &gun, float theta, float csnd = SOUND_SPEED, float gamma = GAMMA, float patm = ATMOSPHERIC_PRESSURE);

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
    void generateNWave(float *Pnw, const float *tInterval, int32_t nPoints, float pmax, float ta, float Td, float tr);

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
    int8_t nWave(float *Pnw, const float *tInterval, int32_t nPoints, const Gunshot::Gun &gun, const Gunshot::Geometry &geometry, float patm = ATMOSPHERIC_PRESSURE, float csnd = SOUND_SPEED);

    /* Check whether an N-wave is observed at the microphone position or not.
     *
     * @param gun       - the gun
     * @param geometry  - the environment geometry structure (coordinates)
     * @param csnd      - the speed of sound
     * @returns         - true if an N-wave is observed; false, otherwise
     */
    bool isObserved(const Gunshot::Gun &gun, const Gunshot::Geometry &geometry, float csnd = SOUND_SPEED);

    /* Calculate the N-wave rise time.
     *
     * @param pmax      - N-wave amplitude in Pa
     * @param patm      - atmospheric pressure in Pa
     * @param csnd      - the speed of sound
     * @param lamb      - air molecular mean free path in m
     * @returns tr      - the rise time in s
     */
    float riseTime(float pmax, float patm = ATMOSPHERIC_PRESSURE, float csnd = SOUND_SPEED, float lamb=AIR_MOLECULAR_MEAN_FREE_PATH);

    /* Calculate the N-wave arrival time.
     *
     * @param r         - the distance between the gun and the microphone
     * @param theta     - the angle between the boreline and the microphone position in radians
     * @param velocity  - the projectile velocity in m/s
     * @param csnd      - the speed of sound in m/s
     * @returns ta      - the time of arrival in s
     */
    float timeOfArrival(float r, float theta, float velocity, float csnd = SOUND_SPEED);

    /* Calculate the N-wave total duration time.
     *
     * @param M          - the Mach number
     * @param bulletDiam - the bullet diameter in m
     * @param bulletLen  - the bullet length in m
     * @param xmiss      - the shortest dist from the mic to the bullet trajectory in m
     * @param csnd       - the speed of sound in m/s
     * @returns Td       - the total duration time in s
     */
    float duration(float M, float bulletDiam, float bulletLen, float xmiss, float csnd = SOUND_SPEED);

    /* Calculate the N-wave amplitude.
     *
     * @param M          - the Mach number
     * @param bulletDiam - the bullet diameter in m
     * @param bulletLen  - the bullet length in m
     * @param patm       - the atmospheric pressure in Pa
     * @param csnd       - the speed of sound in m/s
     * @returns pmax     - the N-wave amplitude in Pa
     */
    float amplitude(float M, float bulletDiam, float bulletLen, float xmiss, float patm = ATMOSPHERIC_PRESSURE);

}  /* namespace NWave */

#endif /* GUNSHOT_H */
