/**
 * Implementation of abstract class 'Integrator'.
 *
 * The 'Integrator' class is used to construct objects that are able to
 * integrate a set of differential equations. This class is used as a
 * superclass for, for example, the RKF45 integrator (RKF45.cpp).
 */

#include <cstdlib>
#include <vector>
#include <softlib/Integrator.h>
#include <softlib/Vector.h>

/**
 * Constructor
 */
template<unsigned int N>
Integrator<N>::Integrator() {
	this->timestep = 1;

	this->ntsallocated = 0;
	this->ntimesteps = 0;
	this->time = nullptr;
	this->solution = nullptr;
	this->derivative = nullptr;
}
template<unsigned int N>
Integrator<N>::Integrator(const Vector<N>& y0) {
	this->ntsallocated = 0;
	this->ntimesteps = 0;
	this->time = nullptr;
	this->solution = nullptr;
	this->derivative = nullptr;

	InitialValue(y0);
}

/**
 * Sets the initial value of the
 * integrator. If the integrator already
 * has taken steps, this resets the
 * integrator and discards any solutions.
 *
 * y0: Initial value of the function at t = 0, x = x0.
 */
template<unsigned int N>
void Integrator<N>::InitialValue(const Vector<N>& y0) {
	Vector<N> yp0;
	yp0 = this->eqn->Evaluate(0.0, y0, yp0);

    // Reset timestep counter
	this->ntimesteps = 0;

	StoreSolution(0.0, y0, yp0);
}

/**
 * Get the last time value
 * which was solved for.
 */
template<unsigned int N>
slibreal_t Integrator<N>::LastTime() {
	if (ntimesteps <= 0) return 0.0;
	else return time[ntimesteps-1];
}

/**
 * Get the solution from the
 * last completed timestep.
 */
template<unsigned int N>
slibreal_t *Integrator<N>::LastSolution() {
	if (ntimesteps <= 0) return nullptr;
	else return solution+((ntimesteps-1)*N);
}

/**
 * Get the derivative evaluated
 * in the most recently solved
 * for timestep.
 */
template<unsigned int N>
slibreal_t *Integrator<N>::LastDerivative() {
	if (ntimesteps <= 0) return nullptr;
	else return derivative+((ntimesteps-1)*N);
}

/**
 * Generate dense output with a uniform timestep
 * between times 'tmin' and 'tmax'.
 *
 * nt:   Number of timesteps to return.
 * tmin: Time in first step to return.
 * tmax: Time in last step to return.
 * out:  Array containing solution on return. Must have nt elements.
 * outt: Array containing time points on return. Must have nt*N elements.
 */
template<unsigned int N>
void Integrator<N>::OutputDense(const unsigned int nt, const slibreal_t tmin, const slibreal_t tmax, slibreal_t *out, slibreal_t *outt) {
    unsigned int i;
    slibreal_t dt = 1;

    if (nt > 1)
        dt = (tmax-tmin) / (nt-1.0);
    else if (nt == 1) {
        if (tmin != tmax)
            throw SOFTLibException("Only one timestep was requested but tmin != tmax.");
    } else
        throw SOFTLibException("Invalid number of timesteps requested: %u.", nt);

    // Generate time array (outt)
    for (i = 0; i < nt; i++) {
        // Use multiplication to avoid
        // accumulating round-off error
        outt[i] = tmin + i*dt;
    }

    // Explicitly set last timestep to prevent round-off errors
    outt[nt-1] = tmax;

    OutputDense(nt, outt, out);
}

/**
 * Sets the equation to solve.
 *
 * eq: IntegratorEquation object representing
 *   the equation ODE solve.
 */
template<unsigned int N>
void Integrator<N>::SetEquation(IntegratorEquation<N> *eq) { this->eqn = eq; }

/**
 * Returns the number of timesteps
 * taken so far.
 */
template<unsigned int N>
unsigned int Integrator<N>::StepsTaken() { return ntimesteps; }

/**
 * Store a solution.
 *
 * t:  Time of the solution.
 * y:  Solution value at t.
 * yp: Function value/derivative at t.
 */
template<unsigned int N>
void Integrator<N>::StoreSolution(const slibreal_t t, const Vector<N>& y, const Vector<N>& yp) {
	if (ntsallocated <= ntimesteps) {
		time = (slibreal_t*)realloc(time, sizeof(slibreal_t)*(ntsallocated+INTEGRATOR_ALLOC_SIZE));
		solution = (slibreal_t*)realloc(solution, sizeof(slibreal_t)*N*(ntsallocated+INTEGRATOR_ALLOC_SIZE));
		derivative = (slibreal_t*)realloc(derivative, sizeof(slibreal_t)*N*(ntsallocated+INTEGRATOR_ALLOC_SIZE));

		ntsallocated += INTEGRATOR_ALLOC_SIZE;
	}

	time[ntimesteps] = t;
	y.ToArray(solution+(ntimesteps*N));
	yp.ToArray(derivative+(ntimesteps*N));

	ntimesteps++;
}

/**
 * Returns the time value corresponding
 * to integrator step i.
 *
 * i: Index of time step to return.
 */
template<unsigned int N>
slibreal_t Integrator<N>::TimeAt(const unsigned int i) {
	return time[i];
}
/**
 * Returns a pointer to the ODE solution
 * corresponding to time index i.
 *
 * i: Index of time step to return.
 */
template<unsigned int N>
slibreal_t *Integrator<N>::SolutionAt(const unsigned int i) {
	return solution+(i*N);
}

/**
 * Returns a pointer to the derivative
 * /equation value corresponding to
 * time index i.
 *
 * i: Index of time step to return.
 */
template<unsigned int N>
slibreal_t *Integrator<N>::DerivativeAt(const unsigned int i) {
	return derivative+(i*N);
}

