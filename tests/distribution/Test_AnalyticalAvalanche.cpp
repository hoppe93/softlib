/**
 * Unit test for the analytic avalanche distribution function.
 */

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/constants.h>
#include "Test_AnalyticalAvalanche.h"

/**
 * Destructor
 */
Test_AnalyticalAvalanche::~Test_AnalyticalAvalanche() {
    if (f != nullptr) {
        delete [] r;
        delete [] p;
        delete [] xi;
        delete [] f;
    }
}

/**
 * Reference implementation of the analytical avalanche
 * distribution function.
 *
 * p:   Particle momentum.
 * xi:  Particle pitch.
 * E:   Electric field, normalized to Ec.
 * lnL: Coulomb logarithm.
 * Z:   Effective plasma charge.
 */
slibreal_t Test_AnalyticalAvalanche::AA(
    slibreal_t p, slibreal_t xi,
    slibreal_t E, slibreal_t lnL, slibreal_t Z
) {
    slibreal_t A, g, g0, p2, prefac;

    p2 = p*p;
    g = sqrt(p2 + 1.0);
    g0 = lnL * sqrt(Z + 5.0);

    A = (E+1.0) / (Z+1.0) * g;

    prefac = ELECTRON_MASS * LIGHTSPEED * A/(2 * PI * g0 * p2);
    return prefac * exp(-g/g0 - A*(1+xi));
}
/**
 * Run the set of tests.
 */
bool Test_AnalyticalAvalanche::Run() {
    unsigned int i, j, k, failed=0;
    slibreal_t EHat, lnLambda, Zeff;
    slibreal_t *r, *p, *xi, *f;
    AnalyticalAvalanche *af = new AnalyticalAvalanche();

    srand(time(NULL));

    // Generate p and xi
    r  = new slibreal_t[NPOINTS];
    p  = new slibreal_t[NPOINTS];
    xi = new slibreal_t[NPOINTS];
    f  = new slibreal_t[NPOINTS];
    for (i = 0; i < NPOINTS; i++) {
        r[i] = 0.0;
        p[i] = PMAX * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
        xi[i] = 1.0 - 2.0 * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
    }

    for (i = 0; i < NEHAT; i++) {
        EHat = EHAT0 + (EHAT1-EHAT0) * ((slibreal_t)i) / ((slibreal_t)NEHAT-1.0);
        for (j = 0; j < NLOGLAMBDA; j++) {
            lnLambda = LOGLAMBDA0 + (LOGLAMBDA1-LOGLAMBDA0) * ((slibreal_t)j) / ((slibreal_t)NLOGLAMBDA-1.0);
            for (k = 0; k < NZEFF; k++) {
                Zeff = ZEFF0 + (ZEFF1-ZEFF0) * ((slibreal_t)k) / ((slibreal_t)NZEFF-1.0);

                af->Initialize(EHat, lnLambda, Zeff);
                // Generate f
                for (i = 0; i < NPOINTS; i++) {
                    f[i] = AA(p[i], xi[i], EHat, lnLambda, Zeff);
                }

                failed += CompareToTable(af, r, p, xi, f, NPOINTS, TOLERANCE);
            }
        }
    }

    if (failed > 0) {
        unsigned int ntot = NPOINTS*NZEFF*NLOGLAMBDA*NEHAT;
        slibreal_t perc = ((slibreal_t)failed) / ((slibreal_t)ntot) * 100.0;
        this->PrintError("Analytical avalanche distribution: %u of %u (%.2f%%) tests had an error larger than %d eps.", failed, ntot, perc, (int)(TOLERANCE/REAL_EPSILON));
        delete af;
        return false;
    } else {
        delete af;
        return true;
    }
}

