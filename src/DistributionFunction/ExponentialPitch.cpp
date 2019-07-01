/**
 * Implementation of a pitch angle distribution function that
 * takes the form
 *
 *   f(\xi) = exp(C\xi)
 */

#include <cmath>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/ExponentialPitch.h>
#include <softlib/DistributionFunction/DistributionFunction.h>

/**
 * Constructor.
 */
ExponentialPitch::ExponentialPitch() { }
ExponentialPitch::ExponentialPitch(const slibreal_t C) : C(C) {}

/**
 * Evaluate the distribution function the
 * given phase-space point.
 *
 * r:           Radial point to evaluate at
 *              (not used in this implementation).
 * p:           Momentum point to evaluate at (in units of mc).
 *              (not used in the implementation)
 * xi:          Pitch point to evaluate at.
 * drift_shift: Amount by which the magnetic
 *              drift surfaces are shifted (not used in
 *              this implementation).
 */
slibreal_t ExponentialPitch::Eval(const slibreal_t, const slibreal_t xi) {
    return exp(C * abs(xi));
}

/**
 * Initialize the distribution function.
 *
 * C: Pitch angle parameter.
 */
void ExponentialPitch::Initialize(const slibreal_t C) { this->C = C; }

/**
 * Clone this analytic avalanche distribution function.
 */
ExponentialPitch *ExponentialPitch::MinClone() {
    ExponentialPitch *da = new ExponentialPitch();
    da->Initialize(this->C);

    return da;
}

