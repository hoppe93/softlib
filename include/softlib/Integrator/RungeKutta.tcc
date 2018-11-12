/**
 * Implementation of a Dormand-Prince fifth order Runge-Kutta method.
 */

#include <cmath>
#include <softlib/config.h>
#include <softlib/RungeKutta.h>
#include <softlib/Vector.h>
#include <softlib/SOFTLibException.h>

#ifndef SQR
#	define SQR(x) ((x)*(x))
#endif

/**
 * Constructor
 *
 * tol: Relative tolerance to use (by default,
 *    the absolute tolerance is set to the same value)
 */
template<unsigned int N>
RungeKutta<N>::RungeKutta(const slibreal_t tol) {
	this->atol = tol;
	this->rtol = tol;
	this->errold = 1.0e-4;
	this->optimal_step = 1.0;
}

/**
 * Estimate the error made in
 * the current timestep.
 *
 * y1: Solution in "previous" timestep.
 * y2: Solution in "current" timestep.
 */
template<unsigned int N>
slibreal_t RungeKutta<N>::EstimateError(const Vector<N> &y1, const Vector<N> &y2) {
	unsigned int i;
	slibreal_t err = 0, sk;

	for (i = 0; i < N; i++) {
		sk = atol + rtol*std::max(fabs(y1[i]), fabs(y2[i]));
		err += SQR(yerr[i]/sk);
	}

	return sqrt(err / N);
}

/**
 * Returns true if coefficients for
 * dense output is stored by this
 * RungeKutta solver.
 */
template<unsigned int N>
bool RungeKutta<N>::HasDenseOutput() { return store_dense_coeffs; }

/**
 * Sets the initial value of the
 * integrator. If the integrator already
 * has taken steps, this resets the
 * integrator and discards any solutions.
 *
 * (This function overrides the identical-looking
 *  function of the same name in 'Integrator' so
 *  that the 'RungeKutta::StoreSolution()' function
 *  also stores dense output if necessary)
 *
 * y0: Initial value of the function at t = 0, x = x0.
 */
template<unsigned int N>
void RungeKutta<N>::InitialValue(const Vector<N>& y0) {
	Vector<N> yp0;
	yp0 = this->eqn->Evaluate(0.0, y0, yp0);

    this->ntimesteps = 0;

	StoreSolution(0.0, y0, yp0);
}

/**
 * Advance the integrator one step. The new
 * step is stored in the 'solutions' vector.
 * The time of the new step is returned.
 */
template<unsigned int N>
slibreal_t RungeKutta<N>::Step() {
	slibreal_t T, Tnew, err;
	y = this->LastSolution();
	dydx = this->LastDerivative();

	h = optimal_step;
	T = this->LastTime();

	while (true) {
		yout = InnerStep(T, h, y, dydx);
		err = EstimateError(y, yout);
		if (this->UpdateTimestep(err)) break;

		h = optimal_step;
		if (fabs(h) <= fabs(T)*REAL_EPSILON)
			throw SOFTLibException("Step-size underflow in Runge-Kutta integrator.");
	}

	yold = y;
	dydxold = dydx;

	Tnew = T + h;
	this->timestep = h;
	this->StoreSolution(Tnew, yout, dydxout);

	return Tnew;
}

/**
 * Override of the 'StoreSolution' method
 * of 'Integrator'. This override aims to
 * also build dense output coefficients
 * when output is stored.
 *
 * t:  Time at current point.
 * y:  Function value at current time.
 * yp: Derivative value at current time.
 */
template<unsigned int N>
void RungeKutta<N>::StoreSolution(const slibreal_t t, const Vector<N>& y, const Vector<N>& yp) {
	Integrator<N>::StoreSolution(t, y, yp);

	/* Dense output coefficients should only be stored
	 * when requested, and only after the first timestep
	 * has been taken. */
	if (HasDenseOutput() && this->ntimesteps > 1)
		PrepareDense(this->timestep, yold, y, dydxold, yp);
}

/**
 * Updates the length of the timestep.
 * Returns 'true' if the previously calculated
 * solution satisfies the accuracy requirements,
 * returns false if the step should be redone.
 * 
 * err: Error estimated for previous step.
 */
template<unsigned int N>
bool RungeKutta<N>::UpdateTimestep(const slibreal_t err) {
	slibreal_t scale;

	if (err <= 1.0) {
		if (err == 0.0)
			scale = maxscale;
		else {
			scale = safe*pow(err, -alpha) * pow(errold, beta);
			if (scale < minscale) scale = minscale;
			if (scale > maxscale) scale = maxscale;
		}

		/* If previous step was rejected,
		 * don't allow step to increase */
		if (previous_step_rejected)
			optimal_step = h*std::min(scale, 1.0);
		else
			optimal_step = h*scale;

		errold = std::max(err, 1e-4);
		previous_step_rejected = false;
		return true;
	} else {
		scale = std::max(safe*pow(err, -alpha), minscale);
		optimal_step = h*scale;
		previous_step_rejected = true;
		return false;
	}
}

/***********
 * SETTERS *
 ***********/
template<unsigned int N>
void RungeKutta<N>::SetExponentAlpha(const slibreal_t a) { this->alpha = a; }
template<unsigned int N>
void RungeKutta<N>::SetExponentBeta(const slibreal_t b) { this->beta = b; }
template<unsigned int N>
void RungeKutta<N>::SetMaxScaleFactor(const slibreal_t m) { this->maxscale = m; }
template<unsigned int N>
void RungeKutta<N>::SetMinScaleFactor(const slibreal_t m) { this->minscale = m; }
template<unsigned int N>
void RungeKutta<N>::SetSafetyFactor(const slibreal_t s) { this->safe = s; }
template<unsigned int N>
void RungeKutta<N>::ToggleDenseOutput(const bool d) { this->store_dense_coeffs = d; }

