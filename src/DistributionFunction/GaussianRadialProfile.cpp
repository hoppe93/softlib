/**
 * Implements a radial profile of the form
 *
 *   f(x) = a * exp(-(x-b/c)**2)
 *
 * where x is the normalized radius (range is [0,1])
 * and a,b,c are the gaussian parameters.
 */

#include <cmath>
#include <softlib/config.h>
#include <softlib/DistributionFunction/GaussianRadialProfile.h>
#include <softlib/SOFTLibException.h>

/**
 * Constructor.
 *
 * a,b,c:    Gaussian Parameters.
 *
 */
GaussianRadialProfile::GaussianRadialProfile(const slibreal_t rmin, const slibreal_t rmax,const slibreal_t a, const slibreal_t b, const slibreal_t c) {
    if (rmin >= rmax)
        throw SOFTLibException("Gaussian radial profile: rmin must be less than rmax. rmin = %e, rmax = %e.", rmin, rmax);
    if (a <= 0)
        throw SOFTLibException("Gaussian radial profile: a must be positive. a = %e.", a);
    if (b <= 0)
        throw SOFTLibException("Gaussian radial profile: b must be positive. b = %e.", b);
    if (c <= 0)
        throw SOFTLibException("Gaussian radial profile: c must be positive. c = %e.", c);

    this->rmin = rmin;
    this->rmax = rmax;
    this->a = a;
    this->b = b;
    this->c = c;
}

/**
 * Evaluate the radial profile.
 *
 * r:           Radial position of particle.
 * drift_shift: Orbit drift shift.
 */
slibreal_t GaussianRadialProfile::Eval(const slibreal_t r, const slibreal_t drift_shift) {
    slibreal_t rho = r-drift_shift;

    if (rho >= this->rmax || rho < this->rmin)
        return 0.0;
    else
        return (this->a) * exp(-(rho-this->b)*(rho-this->b)/ (this->c * this->c));
}

/**
 * Clone this radial profile.
 */
GaussianRadialProfile *GaussianRadialProfile::MinClone() {
    return new GaussianRadialProfile(this->rmin, this->rmax,this->a, this->b, this->c);
}
