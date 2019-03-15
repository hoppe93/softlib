/**
 * Implementation of general unit tests for the
 * numerical distribution function types.
 */

#include <cmath>
#include <ctime>

#include <softlib/config.h>
#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/LinearRadialProfile.h>
#include "Test_NumericDistributionFunction.h"

/**
 * Generate an analytical distribution function consisting
 * of two parts:
 *
 *   Momentum-space: Analytical avalanche distribution
 *   Real-space:     Quadratic function
 *
 * RETURNS a new RadialDistributionFunction object.
 */
RadialDistributionFunction *Test_NumericDistributionFunction::GeneratePhaseSpaceDistribution(
    slibreal_t EHat, slibreal_t logLambda, slibreal_t Zeff,
    slibreal_t rmin, slibreal_t rmax
) {
    RadialDistributionFunction *rdf = new RadialDistributionFunction();
    LinearRadialProfile *rp = new LinearRadialProfile(rmin, rmax);
    AnalyticalAvalanche *aa = new AnalyticalAvalanche();

    aa->Initialize(EHat, logLambda, Zeff);
    rdf->Initialize(rp, aa);

    return rdf;
}

/**
 * Test a set of 'n' randomised and compare the values
 * of two distribution functions in those points. Note
 * that 'df1' is assumed to be the "correct" distribution
 * function.
 *
 * df1:  Distribution function one (assumed to be the one to check against).
 * df2:  Distribution function two (assumed to be the one to check).
 * n:    Number of test points to generate.
 * tol:  Tolerance (counted as fail if relative error is above this value).
 * rmin: Minimum allowed radius.
 * rmax: Maximum allowed radius.
 * pmin: Minimum allowed momentum.
 * pmax: Maximum allowed momentum.
 *
 * RETURNS the number of failed tests.
 */
unsigned int Test_NumericDistributionFunction::RandomPoints(
    DistributionFunction *df1, DistributionFunction *df2,
    const unsigned int n, const slibreal_t tol,
    const slibreal_t rmin, const slibreal_t rmax,
    const slibreal_t pmin, const slibreal_t pmax
) {
    slibreal_t *r, *p, *xi, *f;
    unsigned int i, failed=0;

    srand(time(NULL));

    // Generate test points
    r  = new slibreal_t[n];
    p  = new slibreal_t[n];
    xi = new slibreal_t[n];
    f  = new slibreal_t[n];

    for (i = 0; i < n; i++) {
        r[i]  = rmin + (rmax-rmin) * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
        p[i]  = pmin + (pmax-pmin) * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
        xi[i] = 1.0 - 2.0 * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);

        f[i]  = df1->Eval(r[i], p[i], xi[i]);
    }

    failed = this->CompareToTable(df2, r, p, xi, f, n, tol);

    delete [] r;
    delete [] p;
    delete [] xi;
    delete [] f;

    return failed;
}

