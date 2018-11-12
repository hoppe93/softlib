/**
 * Implementation of a general integrator
 * equation class.
 */

#include <softlib/config.h>
#include <softlib/GeneralEquation.h>
#include <softlib/IntegratorEquation.h>

template<unsigned int N>
GeneralEquation<N>::GeneralEquation(Vector<N>& (*eqn)(const slibreal_t, const Vector<N>&, Vector<N>&)) {
	this->equation = eqn;
}

/**
 * Evaluate the equation given 
 * to this class.
 *
 * t: Time at requested point.
 * y: Function value in t.
 * dydxout: On return, contains function value
 *    in the requested point (also returned).
 */
template<unsigned int N>
Vector<N>& GeneralEquation<N>::Evaluate(const slibreal_t t, const Vector<N>& y, Vector<N>& dydxout) {
	return equation(t, y, dydxout);
}

/**
 * Sets which function to use when
 * evaluating the derivative (RHS)
 * of the ODE to solve.
 */
template<unsigned int N>
void GeneralEquation<N>::SetDerivative(Vector<N>& (*eqn)(const slibreal_t, const Vector<N>&, Vector<N>&)) {
	this->equation = eqn;
}

