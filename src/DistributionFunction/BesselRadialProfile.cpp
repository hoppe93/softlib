/**
 * Implements a radial profile of the form
 *
 *   f(r/a) = J_0(r/a)
 * 
 * where r/a is the normalized minor radius (range is [0,1]).
 */

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <softlib/config.h>
#include <softlib/DistributionFunction/BesselRadialProfile.h>
#include <softlib/SOFTLibException.h>

/**
 * Constructor.
 *
 * rmin: Minimum allowed radius.
 * rmax: Maximum allowed radius.
 *
 * This radial profile is identically zero outside
 * the bounds rmin and rmax.
 */
BesselRadialProfile::BesselRadialProfile(const slibreal_t rmin, const slibreal_t rmax) {
    if (rmin >= rmax)
        throw SOFTLibException("Bessel radial profile: rmin must be less than rmax. rmin = %e, rmax = %e.", rmin, rmax);

    this->rmin = rmin;
    this->rmax = rmax;
}

/**
 * Evaluate the radial profile.
 *
 * r:           Radial position of particle.
 * drift_shift: Orbit drift shift.
 */
slibreal_t BesselRadialProfile::Eval(const slibreal_t r, const slibreal_t drift_shift) {
    slibreal_t rho = r-drift_shift;

    if (rho >= this->rmax || rho < this->rmin)
        return 0.0;
    else {
        slibreal_t x = (rho-this->rmin)/(this->rmax-this->rmin) * J0_ZERO1;
        return gsl_sf_bessel_J0(x);
    }
}

/**
 * Clone this radial profile.
 */
BesselRadialProfile *BesselRadialProfile::MinClone() {
    return new BesselRadialProfile(this->rmin, this->rmax);
}

