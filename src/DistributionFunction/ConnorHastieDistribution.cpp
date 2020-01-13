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
    const slibreal_t EHat, const slibreal_t Zeff,
    const slibreal_t pMax, const slibreal_t deltaP
) {
    this->Initialize(EHat, Zeff, pMax, deltaP);
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
    if (deltaP == 0 || p < this->pMax)
        return (1.0/(p*abs(xi)) * exp(-this->D * p * (1-xi*xi)/abs(xi)));
    else {  // Apply Gaussian cut-off
        slibreal_t del = (p-this->pMax) / this->deltaP;
        return (1.0/(p*abs(xi)) * exp(-this->D * p * (1-xi*xi)/abs(xi) - del*del));
    }
}

/**
 * Initialize the distribution function object.
 *
 * EHat:     Electric field strength, in units of the
 *    threshold electric field introduced by
 *    [Connor & Hastie, Nucl. Fusion 15, 415 (1975)].
 * Zeff:     Effective charge of the plasma.
 * pMax:     Momentum at which Gaussian cut-off should
 *           be applied.
 * deltaP:   Rate at which the Gaussian cut-off decays.
 */
void ConnorHastieDistribution::Initialize(
    const slibreal_t EHat, const slibreal_t Zeff,
    const slibreal_t pMax, const slibreal_t deltaP
) {
    this->EHat   = EHat;
    this->Zeff   = Zeff;
    this->pMax   = pMax;
    this->deltaP = deltaP;
    this->D    = (this->EHat + 1) / (2*(this->Zeff + 1));
}

/**
 * Clone this analytic avalanche distribution function.
 */
ConnorHastieDistribution *ConnorHastieDistribution::MinClone() {
    return new ConnorHastieDistribution(this->EHat, this->Zeff, this->pMax, this->deltaP);
}

