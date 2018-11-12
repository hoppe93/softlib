#ifndef _TEST_INTEGRATOR_H
#define _TEST_INTEGRATOR_H

#include <softlib/config.h>
#include <softlib/GeneralEquation.h>
#include <softlib/Integrator.h>
#include <softlib/Vector.h>

#include "runtest.h"

class Test_Integrator : public UnitTest {
	protected:
		static slibreal_t ode1result;
		static Vector<1>& ode1eqn(const slibreal_t, const Vector<1>&, Vector<1>&);
		static slibreal_t *ode1val(const slibreal_t, const Vector<1>&);

		static slibreal_t ode2result[2];
		static Vector<2>& ode2eqn(const slibreal_t, const Vector<2>&, Vector<2>&);
		static slibreal_t *ode2val(const slibreal_t, const Vector<2>&);

		slibreal_t tolerance;
	public:
		Test_Integrator(const string&);

		GeneralEquation<1> *GetODE1Equation();
		GeneralEquation<2> *GetODE2Equation();

		bool RunODE1Test(Integrator<1>&);
		bool RunODE2Test(Integrator<2>&);
		template<unsigned int N>
		bool RunODETest(Integrator<N>&, const slibreal_t, const Vector<N>&, slibreal_t* (*odeval)(slibreal_t, const Vector<N>&));

		void SetTolerance(const slibreal_t);
};

#endif/*_TEST_INTEGRATOR_H*/
