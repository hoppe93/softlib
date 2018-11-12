/**
 * Unit test for the Dormand-Prince fifth-order
 * Runge-Kutta integrator.
 */

#include <cmath>
#include <string>

#include <softlib/RKDP45.h>
#include <softlib/Vector.h>

#include "rkdp45.h"
#include "runtest.h"

using namespace std;

/**
 * Constructor
 */
Test_RKDP45::Test_RKDP45() : Test_Integrator("int_rkdp45") {}
Test_RKDP45::Test_RKDP45(const string& name) : Test_Integrator(name) {}

/**
 * Run the tests of the RKDP45 implementation.
 */
bool Test_RKDP45::Run() {
	bool success = true;
	RKDP45<1> ode1;
	RKDP45<2> ode2;

	/* Test 1D ODE */
	if ((success=RunODE1Test(ode1)))
		this->PrintOK("1-D ODE successfully solved.");
	else this->PrintError("[int_rkdp45]: Failed to solve 1-D ODE.");
	if ((success=RunODE1TestDense(ode1)))
		this->PrintOK("1-D ODE successfully solved with dense output.");
	else this->PrintError("[int_rkdp45]: Failed to solve 1-D ODE with dense output.");
	
	/* Test 2D ODE */
	if ((success=RunODE2Test(ode2)))
		this->PrintOK("2-D ODE successfully solved.");
	else this->PrintError("[int_rkdp45]: Failed to solve 2-D ODE.");
	if ((success=RunODE2TestDense(ode2)))
		this->PrintOK("2-D ODE successfully solved with dense output.");
	else this->PrintError("[int_rkdp45]: Failed to solve 2-D ODE with dense output.");
	
	return success;
}

/**
 * Runs through the test procedure for the ODEs.
 * Returns true if the solution is within the
 * desired tolerance, and false otherwise.
 *
 * solver: 'Integrator' object to use for solving the ODE.
 */
bool Test_RKDP45::RunODE1TestDense(RKDP45<1>& solver) {
	slibreal_t tmax = 10.0;
	Vector<1> x0;
	GeneralEquation<1> *eqn;
	x0[0] = 1.0;

	eqn = GetODE1Equation();
	solver.SetEquation(eqn);
	solver.InitialValue(x0);

	return RunODETestDense(solver, tmax, x0, ode1val);
}
bool Test_RKDP45::RunODE2TestDense(RKDP45<2>& solver) {
	slibreal_t tmax = 4.0;
	Vector<2> x0;
	GeneralEquation<2> *eqn;
	x0[0] = 1.0;
	x0[1] = 0.0;

	eqn = GetODE2Equation();
	solver.SetEquation(eqn);
	solver.InitialValue(x0);

	return RunODETestDense(solver, tmax, x0, ode2val);
}

/**
 * Generic routine for solving an ODE
 * and comparing the dense result to a
 * previously known solution.
 *
 * solver: Integrator object to use in solving the ODE.
 * tmax: Time to run until.
 * x0: Initial condition for the solution.
 * odeval: Correct solution to the ODE.
 */
template<unsigned int N>
bool Test_RKDP45::RunODETestDense(
	RKDP45<N>& solver, const slibreal_t tmax, const Vector<N>& x0,
	slibreal_t* (*odeval)(slibreal_t, const Vector<N>&)
) {
	unsigned int i, j, nsteps;
	slibreal_t *ssol, *rsol, *tp, absv;

	/* Solve equation */
	while (solver.Step() < tmax);

	/* Get dense solution */
	nsteps = 200;
	ssol = new slibreal_t[nsteps*N];
	tp = new slibreal_t[nsteps];

	solver.OutputDense(nsteps, 0.0, tmax, ssol, tp);

	for (i = 0; i < nsteps; i++) {
		rsol = odeval(tp[i], x0);
		absv = 0.0;

		for (j = 0; j < N; j++)
			absv += rsol[j]*rsol[j];

		absv = sqrt(absv);

		for (j = 0; j < N; j++) {
			if (fabs((ssol[i*N+j]-rsol[j])/absv) > tolerance) {
				printf("failing at i = %d, where error = %e\n", i, fabs((ssol[i*N+j]-rsol[j])/absv));
				delete [] ssol;
				delete [] tp;
				return false;
			}
		}
	}

	delete [] ssol;
	delete [] tp;
	return true;
}

