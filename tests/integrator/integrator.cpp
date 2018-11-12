/**
 * Implements the general ODE integrator test object.
 */

#include <cmath>

#include <softlib/GeneralEquation.h>
#include <softlib/Integrator.h>
#include <softlib/Vector.h>

#include "integrator.h"

slibreal_t Test_Integrator::ode1result = 0;
slibreal_t Test_Integrator::ode2result[2] = {0.0,0.0};

/*****************
 * ODE FUNCTIONS *
 *****************/
/**
 * Evaluates the RHS of the equation
 *
 *   dx/dt = t/x,
 *
 * which is solved by
 *   
 *   x(t) = x(0) * sqrt(1 + t^2).
 *
 */
Vector<1>& Test_Integrator::ode1eqn(const slibreal_t T, const Vector<1>& y, Vector<1>& dydxout) {
	dydxout[0] = T / y[0];
	return dydxout;
}
/**
 * Evaluates the solution x(t) corresponding
 * to the above equation.
 *
 * T:  Time point at which to evaluate the solution.
 * x0: Initial value, x(t=0).
 */
slibreal_t *Test_Integrator::ode1val(const slibreal_t T, const Vector<1>& x0) {
	Test_Integrator::ode1result = x0[0] * sqrt(1.0 + T*T);
	return &(Test_Integrator::ode1result);
}

/**
 * Evaluates the RHSs of the
 * equation system
 *
 *   dx/dt =-a*y,
 *   dy/dt = a*x,
 *
 * which is solved by
 *   
 *   x(t) = x(0)*cos(a*t) + y(0)*sin(a*t),
 *   y(t) = x(0)*sin(a*t) - y(0)*cos(a*t).
 */
Vector<2>& Test_Integrator::ode2eqn(const slibreal_t T __attribute__((unused)), const Vector<2>& y, Vector<2>& dydxout) {
	static slibreal_t a = 3.14159265359;

	dydxout[0] =-a*y[1];
	dydxout[1] = a*y[0];

	return dydxout;
}
/**
 * Evaluates the solution (x(t), y(t)) corresponding
 * to the above equation.
 *
 * T:  Time point at which to evaluate the solution.
 * x0: Initial value, (x(t=0), y(t=0)).
 */
slibreal_t *Test_Integrator::ode2val(const slibreal_t T, const Vector<2>& x0) {
	static slibreal_t a = 3.14159265359;

	Test_Integrator::ode2result[0] = x0[0] * cos(a*T) + x0[1] * sin(a*T);
	Test_Integrator::ode2result[1] = x0[0] * sin(a*T) - x0[1] * cos(a*T);

	return Test_Integrator::ode2result;
}

/**
 * Returns an 'IntegratorEquation' object
 * representing the test ODEs.
 */
GeneralEquation<1> *Test_Integrator::GetODE1Equation() { return new GeneralEquation<1>(ode1eqn); }
GeneralEquation<2> *Test_Integrator::GetODE2Equation() { return new GeneralEquation<2>(ode2eqn); }

Test_Integrator::Test_Integrator(const string& name) : UnitTest(name) {
	tolerance = 20*RK_DEFAULT_TOLERANCE;
}

/**
 * Runs through the test procedure for the ODEs.
 * Returns true if the solution is within the
 * desired tolerance, and false otherwise.
 *
 * solver: 'Integrator' object to use for solving the ODE.
 */
bool Test_Integrator::RunODE1Test(Integrator<1>& solver) {
	slibreal_t tmax = 10.0;
	Vector<1> x0;
	GeneralEquation<1> *eqn;
	x0[0] = 1.0;

	eqn = GetODE1Equation();
	solver.SetEquation(eqn);
	solver.InitialValue(x0);

	return RunODETest(solver, tmax, x0, ode1val);
}
bool Test_Integrator::RunODE2Test(Integrator<2>& solver) {
	slibreal_t tmax = 4.0;
	Vector<2> x0;
	GeneralEquation<2> *eqn;
	x0[0] = 1.0;
	x0[1] = 0.0;

	eqn = GetODE2Equation();
	solver.SetEquation(eqn);
	solver.InitialValue(x0);

	return RunODETest(solver, tmax, x0, ode2val);
}

/**
 * Generic routine for solving an ODE
 * and comparing the result to a
 * previously known solution.
 *
 * solver: Integrator object to use in solving the ODE.
 * tmax: Time to run until.
 * x0: Initial condition for the solution.
 * odeval: Correct solution to the ODE.
 */
template<unsigned int N>
bool Test_Integrator::RunODETest(
	Integrator<N>& solver, const slibreal_t tmax, const Vector<N>& x0,
	slibreal_t* (*odeval)(slibreal_t, const Vector<N>&)
) {
	unsigned int i, j, nsteps;
	slibreal_t *ssol, *rsol, T, absv;

	/* Solve equation */
	while (solver.Step() < tmax);

	nsteps = solver.StepsTaken();
	for (i = 0; i < nsteps; i++) {
		T = solver.TimeAt(i);
		ssol = solver.SolutionAt(i);
		rsol = odeval(T, x0);
		absv = 0.0;

		for (j = 0; j < N; j++)
			absv += rsol[j]*rsol[j];

		absv = sqrt(absv);

		for (j = 0; j < N; j++) {
			if (fabs((ssol[j]-rsol[j])/absv) > tolerance)
				return false;
		}
	}

	return true;
}

/**
 * Sets the tolerance of the similarity
 * test. For a test to pass, the relative
 * difference between the calculated solution
 * and the correct solution must be less
 * than this tolerance.
 *
 * tol: Tolerance value to use.
 */
void Test_Integrator::SetTolerance(const slibreal_t tol) {
	tolerance = tol;
}

