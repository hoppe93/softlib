/**
 * Implementation of the linear radial profile.
 */

#include <softlib/config.h>
#include <softlib/DistributionFunction/LinearRadialProfile.h>
#include <softlib/SOFTLibException.h>

/**
 * Constructor.
 *
 * rmin: Innermost radius.
 * rmax: Outermost radius.
 *
 * Outside the interval [rmin, rmax] the radial profile
 * will be identically zero.
 */
LinearRadialProfile::LinearRadialProfile(const slibreal_t rmin, const slibreal_t rmax) {
    if (rmin >= rmax)
        throw SOFTLibException("Linear radial profile: rmin must be less than rmax. rmin = %e, rmax = %e.", rmin, rmax);

    this->rmin = rmin;
    this->rmax = rmax;
    this->delr = 1.0 / (rmax-rmin);
}

/**
 * Evaluate the radial profile in the given point.
 *
 * r:           Radial position of particle.
 * drift_shift: Orbit drift shift.
 */
slibreal_t LinearRadialProfile::Eval(const slibreal_t r, const slibreal_t drift_shift) {
    slibreal_t rho = r - drift_shift;

    if (rho >= this->rmax || rho < this->rmin)
        return 0.0;
    else
        return 1.0 - (rho - this->rmin) * this->delr;
}

/**
 * Clone this radial profile.
 */
LinearRadialProfile *LinearRadialProfile::MinClone() {
    return new LinearRadialProfile(this->rmin, this->rmax);
}
