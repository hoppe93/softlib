/**
 * Implementation of safety factor methods for the analytical
 * 2D magnetic field. These methods evaluate and differentiate
 * the safety factor q(r).
 */

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>

/**
 * Compute the safety factor in the given point
 *
 * r: Minor radius at which to evaluate the
 *    safety factor.
 */
slibreal_t MagneticFieldAnalytical2D::__GetSafetyFactor(slibreal_t r) {
	slibreal_t rn = r/rminor;

	switch (safety_factor_type) {
		case MFASF_CONSTANT:
			return safety_factor_param1;
		case MFASF_LINEAR:
			return safety_factor_param1*rn + 1.0;
		case MFASF_QUADRATIC:
			return safety_factor_param1*rn*rn + 1.0;
		case MFASF_EXPONENTIAL:
			return exp(rn * safety_factor_param1);
		default: return 1.0;
	}
}
/**
 * Compute the derivative of q(r) with respect
 * to minor radius r, and multiply it by r.
 *
 * r: Minor radius at which to evaluate.
 */
slibreal_t MagneticFieldAnalytical2D::__GetrDqDr(slibreal_t r) {
	slibreal_t rn = r / rminor;

	switch (safety_factor_type) {
		case MFASF_CONSTANT:
			return 0;
		case MFASF_LINEAR:
			return safety_factor_param1 * rn;
		case MFASF_QUADRATIC:
			return 2*safety_factor_param1 * rn * rn;
		case MFASF_EXPONENTIAL:
			return exp(rn * safety_factor_param1) * safety_factor_param1 * rn;
		default: return 0;
	}
}

