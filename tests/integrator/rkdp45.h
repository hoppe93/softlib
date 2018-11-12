#ifndef _TEST_RKDP45_H
#define _TEST_RKDP45_H

#include <softlib/RKDP45.h>

#include "integrator.h"
#include "runtest.h"

class Test_RKDP45 : public Test_Integrator {
	private:
		
	public:
		Test_RKDP45();
		Test_RKDP45(const string&);
		bool Run();

		bool RunODE1TestDense(RKDP45<1>&);
		bool RunODE2TestDense(RKDP45<2>&);
		template<unsigned int N>
		bool RunODETestDense(RKDP45<N>&, const slibreal_t, const Vector<N>&, slibreal_t* (*odeval)(slibreal_t, const Vector<N>&));
};

#endif/*_TEST_RKDP45_H*/
