/**
 * Implements a radial profile of the form
 *
 *   f(a) = 1 - a^b
 * 
 * where a is the normalized radius (range is [0,1])
 * and b is the exponent parameter.
 */

#include <cmath>
#include <softlib/config.h>
#include <softlib/DistributionFunction/PowerRadialProfile.h>
#include <softlib/SOFTLibException.h>

/**
 * Constructor.
 *
 * rmin: Minimum allowed radius.
 * rmax: Maximum allowed radius.
 * b:    Exponent of the radial profile.
 *
 * This radial profile is identically zero outside
 * the bounds rmin and rmax.
 */
PowerRadialProfile::PowerRadialProfile(const slibreal_t rmin, const slibreal_t rmax, const slibreal_t b) {
    if (rmin >= rmax)
        throw SOFTLibException("Power radial profile: rmin must be less than rmax. rmin = %e, rmax = %e.", rmin, rmax);

    this->rmin = rmin;
    this->rmax = rmax;
    this->b    = b;
    this->delr = 1.0 / (rmax-rmin);
}

/**
 * Evaluate the radial profile.
 *
 * r:           Radial position of particle.
 * drift_shift: Orbit drift shift.
 */
slibreal_t PowerRadialProfile::Eval(const slibreal_t r, const slibreal_t drift_shift) {
    slibreal_t rho = r-drift_shift;

    if (rho >= this->rmax || rho < this->rmin)
        return 0.0;
    else
        return 1.0 - pow((rho-this->rmin)*this->delr, this->b);
}

/**
 * Clone this radial profile.
 */
PowerRadialProfile *PowerRadialProfile::MinClone() {
    return new PowerRadialProfile(this->rmin, this->rmax, this->b);
}
