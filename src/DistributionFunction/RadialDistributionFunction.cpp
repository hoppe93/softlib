/**
 * Implementation of the radial distribution function.
 */

#include <softlib/config.h>
#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/RadialProfile.h>

/**
 * Constructor.
 */
RadialDistributionFunction::RadialDistributionFunction() { }
RadialDistributionFunction::RadialDistributionFunction(
    RadialProfile *rp, MomentumSpaceDistributionFunction *msdf
) {
    this->Initialize(rp, msdf);
}

/**
 * Evaluate the radial distribution function.
 *
 * rho:         Radial location to evaluate the function in.
 * p:           Momentum of the particle.
 * xi:          Particle pitch.
 * drift_shift: Orbit drift shift of the particle.
 */
slibreal_t RadialDistributionFunction::Eval(const slibreal_t rho, const slibreal_t p, const slibreal_t xi, const slibreal_t drift_shift) {
    return this->rp->Eval(rho, drift_shift) * this->msdf->Eval(p, xi);
}

/**
 * Initialize the radial distribution function with the
 * given radial profile and momentum-space distribution
 * function.
 *
 * rp:   Function to use as radial profile.
 * msdf: Momentum-space distribution function to use.
 */
void RadialDistributionFunction::Initialize(RadialProfile *rp, MomentumSpaceDistributionFunction *msdf) {
    this->rp = rp;
    this->msdf = msdf;
}

/**
 * Clone this distribution function while re-using
 * as many as possible of allocated objects.
 */
RadialDistributionFunction *RadialDistributionFunction::MinClone() {
    RadialDistributionFunction *rdist = new RadialDistributionFunction();
    RadialProfile *r = rp->MinClone();
    MomentumSpaceDistributionFunction *m = msdf->MinClone();

    rdist->Initialize(r, m);
    return rdist;
}

