/**
 * Implementation of the quasi-steady-state spatially homogeneous
 * runaway-electron distribution function originally given in
 * [F\"ul\"op et al, Phys. Plasmas 13, 062506 (2006)] and modified
 * in [Embr\'eus et al, in preparation].
 *
 * The formula implemented is:
 *
 *   f(p, \xi) = m_e*c*A(p)/(2\pi*\gamma_0*p^2)
 *        exp(-\gamma/\gamma_0 - A(p)*(1 - |\xi|)) *
 *        1 / (1 - exp(-2*A(p))),
 * 
 * where m_e is the electron mass, c the speed of light,
 *
 *   A(p)     = (E/E_c + 1) / (Z_eff + 1) * \gamma,
 *   \gamma_0 = log(\Lambda) sqrt( Z_eff + 5 ).
 */

#include <cmath>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/DistributionFunction/DistributionFunction.h>

/**
 * Constructor.
 */
AnalyticalAvalanche::AnalyticalAvalanche() { }
AnalyticalAvalanche::AnalyticalAvalanche(
    const slibreal_t EHat, const slibreal_t lnLambda, const slibreal_t Zeff
) {
    this->Initialize(EHat, lnLambda, Zeff);
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
slibreal_t AnalyticalAvalanche::Eval(const slibreal_t p, const slibreal_t xi) {
    slibreal_t A, gmm, f, p2;

    p2 = p*p;
    gmm = sqrt(p2 + 1.0);
    A = Ac * gmm;

    f  = ELECTRON_MASS*LIGHTSPEED * A / (2.0 * PI * gmm0 * p2);
    f *= exp(-gmm/gmm0 - A * (1.0 - abs(xi))) / (1.0 - exp(-2.0*A));

    return f;
}

/**
 * Initialize the analytical avalanche
 * distribution function object.
 *
 * EHat:     Electric field strength, in units of the
 *    threshold electric field introduced by
 *    [Connor & Hastie, Nucl. Fusion 15, 415 (1975)].
 * lnLambda: Coulomb logarithm.
 * Zeff:     Effective charge of the plasma.
 */
void AnalyticalAvalanche::Initialize(
    const slibreal_t EHat, const slibreal_t lnLambda, const slibreal_t Zeff
) {
    this->EHat = EHat;
    this->lnLambda = lnLambda;
    this->Zeff = Zeff;

    this->Ac = (EHat + 1.0) / (Zeff + 1.0);
    this->gmm0 = lnLambda * sqrt(Zeff + 5.0);
}

/**
 * Clone this analytic avalanche distribution function.
 */
AnalyticalAvalanche *AnalyticalAvalanche::MinClone() {
    AnalyticalAvalanche *da = new AnalyticalAvalanche();
    da->Initialize(this->EHat, this->lnLambda, this->Zeff);

    return da;
}

