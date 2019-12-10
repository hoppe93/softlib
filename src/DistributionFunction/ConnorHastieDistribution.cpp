/**
 * Implementation of Eq. (60) in [Connor and Hastie, Nucl. Fusion 15, 415 (1975)]
 * which gives the tail of the distribution function in the limit
 * pperp << ppar:
 *
 *   f = 1/ppar * exp[ -(alpha+1)/(2*(Zeff+1)) * pperp^2 / ppar ]
 */

#include <cmath>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/ConnorHastieDistribution.h>
#include <softlib/DistributionFunction/DistributionFunction.h>

/**
 * Constructor.
 */
ConnorHastieDistribution::ConnorHastieDistribution() { }
ConnorHastieDistribution::ConnorHastieDistribution(
    const slibreal_t EHat, const slibreal_t Zeff
) {
    this->Initialize(EHat, Zeff);
}

/**
 * Evaluate the distribution function the
 * given phase-space point.
 *
 * r:           Radial point to evaluate at
 *    (not used in this implementation).
 * p:           Momentum point to evaluate at (in units of mc).
 * xi:          Pitch point to evaluate at.
 * drift_shift: Amount by which the magnetic
 *    drift surfaces are shifted (not used in
 *    this implementation).
 */
slibreal_t ConnorHastieDistribution::Eval(const slibreal_t p, const slibreal_t xi) {
    return (1.0/(p*abs(xi)) * exp(-this->D * p * (1-xi*xi)/abs(xi)));
}

/**
 * Initialize the distribution function object.
 *
 * EHat:     Electric field strength, in units of the
 *    threshold electric field introduced by
 *    [Connor & Hastie, Nucl. Fusion 15, 415 (1975)].
 * Zeff:     Effective charge of the plasma.
 */
void ConnorHastieDistribution::Initialize(
    const slibreal_t EHat, const slibreal_t Zeff
) {
    this->EHat = EHat;
    this->Zeff = Zeff;
    this->D    = (this->EHat + 1) / (2*(this->Zeff + 1));
}

/**
 * Clone this analytic avalanche distribution function.
 */
ConnorHastieDistribution *ConnorHastieDistribution::MinClone() {
    return new ConnorHastieDistribution(this->EHat, this->Zeff);
}

