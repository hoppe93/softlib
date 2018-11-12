#ifndef _GENERAL_EQUATION_H
#define _GENERAL_EQUATION_H

#include <softlib/IntegratorEquation.h>

template<unsigned int N>
class GeneralEquation : public IntegratorEquation<N> {
	protected:
		Vector<N>& (*equation)(const slibreal_t, const Vector<N>&, Vector<N>&);
	
	public:
		GeneralEquation(Vector<N>& (*equation)(const slibreal_t, const Vector<N>&, Vector<N>&));
		Vector<N>& Evaluate(const slibreal_t, const Vector<N>&, Vector<N>&);

		void SetDerivative(Vector<N>& (*eqn)(const slibreal_t, const Vector<N>&, Vector<N>&));
};

// Implementation
#include <softlib/Integrator/GeneralEquation.tcc>

#endif/*_GENERAL_EQUATION_H*/
