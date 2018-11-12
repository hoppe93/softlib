#ifndef _INTEGRATOR_EQUATION_H
#define _INTEGRATOR_EQUATION_H

#include <softlib/Vector.h>

template<unsigned int N>
class IntegratorEquation {
	protected:
		/* Vector in which to store
		 * the return value of 'Evaluate' */
		Vector<N> evaluateRetval;
	public:
		virtual Vector<N>& Evaluate(const slibreal_t, const Vector<N>&, Vector<N>&) = 0;
};

#endif/*_INTEGRATOR_EQUATION_H*/
