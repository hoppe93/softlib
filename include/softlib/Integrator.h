#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include <vector>
#include <softlib/config.h>
#include <softlib/IntegratorEquation.h>
#include <softlib/Vector.h>

/* Number of timesteps to allocate when
 * the storage arrays 'solution', 'derivative',
 * and 'time' run out of memory. */
#define INTEGRATOR_ALLOC_SIZE 1000

template<unsigned int N>
class Integrator {
	protected:
		slibreal_t timestep;

		/* Solution list (y) */
		slibreal_t *solution;
		//Vector<N> *solution;
		/* List of derivatives of y_n,
		 * evaluated at x_n (n = 0,1,2,...) */
		slibreal_t *derivative;
		slibreal_t *time;
		unsigned int ntimesteps, ntsallocated;

		IntegratorEquation<N> *eqn;

	public:
		Integrator();
		Integrator(const Vector<N>&);

		/* Initialize the integrator */
		virtual void InitialValue(const Vector<N>&);
		/* Get the last solution */
		slibreal_t LastTime();
		slibreal_t *LastSolution();
		slibreal_t *LastDerivative();
		unsigned int StepsTaken();

        /* Returns the solution in the given points */
		void OutputDense(const unsigned int, const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);		/* nt, tmin, tmax, out, outt */
        virtual void OutputDense(const unsigned int, slibreal_t*, slibreal_t*) = 0;     /* nt, out, outt */

		/* Advances the integrator one step */
		virtual slibreal_t Step() = 0;

		/* Set the equation to integrate */
		void SetEquation(IntegratorEquation<N>*);
		/* Pushes a timestep to the internal step stack */
		void StoreSolution(const slibreal_t, const Vector<N>&, const Vector<N>&);

		slibreal_t TimeAt(const unsigned int);
		slibreal_t *SolutionAt(const unsigned int);
		slibreal_t *DerivativeAt(const unsigned int);
};

// Implementation
#include <softlib/Integrator/Integrator.tcc>

#endif/*_INTEGRATOR_H*/
