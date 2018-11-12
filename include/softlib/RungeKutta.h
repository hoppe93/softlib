#ifndef _RUNGE_KUTTA_H
#define _RUNGE_KUTTA_H

#include <softlib/config.h>
#include <softlib/Integrator.h>

template<unsigned int N>
class RungeKutta : public Integrator<N> {
	private:
		slibreal_t errold;
		bool previous_step_rejected = false,
			 store_dense_coeffs = true;
	protected:
		/* Vectors in which we store the output
		   from the 'InnerStep' routine. */
		Vector<N> yout, yerr, dydxout, y, dydx,
			yold, dydxold;
		/* Timestep to take in 'InnerStep' */
		slibreal_t h;

		/* Parameters used for rescaling
		   the timestep. */
		slibreal_t
			safe=0.9,				/* Safety factor */
			beta = 0.0,				/* Exponent for PI control (0 = no PI control; 0.04-0.08 is otherwise a good default) */
			alpha = 0.2-beta*0.75,	/* Exponent of rescaling */
			minscale = 0.2,			/* Minimum value to decrease step with */
			maxscale = 10.0;		/* Maximum value to increase step with */

		slibreal_t optimal_step, atol, rtol;
	public:
		RungeKutta(const slibreal_t=RK_DEFAULT_TOLERANCE);

		slibreal_t EstimateError(const Vector<N>&, const Vector<N>&);
		void InitialValue(const Vector<N>&);
		slibreal_t Step();
		void StoreSolution(const slibreal_t, const Vector<N>&, const Vector<N>&);
		bool UpdateTimestep(const slibreal_t);

		/* Virtual methods */
		virtual Vector<N>& InnerStep(const slibreal_t, const slibreal_t, const Vector<N>&, const Vector<N>&) = 0;
		virtual void PrepareDense(const slibreal_t, const Vector<N>&, const Vector<N>&, const Vector<N>&, const Vector<N>&) = 0;

		bool HasDenseOutput();

		/* Setters */
		void SetExponentAlpha(const slibreal_t);
		void SetExponentBeta(const slibreal_t);
		void SetMaxScaleFactor(const slibreal_t);
		void SetMinScaleFactor(const slibreal_t);
		void SetSafetyFactor(const slibreal_t);
		void ToggleDenseOutput(const bool);
};

// Implementation
#include <softlib/Integrator/RungeKutta.tcc>

#endif/*_RUNGE_KUTTA_H*/
