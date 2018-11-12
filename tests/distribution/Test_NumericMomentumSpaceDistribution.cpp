/**
 * Unit test for the momentum-space distribution function class.
 */

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>

#include "Test_DistributionFunction.h"
#include "Test_NumericMomentumSpaceDistribution.h"

/**
 * Generate a numerical distribution function from
 * the analytical avalanche distribution, which can
 * be used to test the interpolation routines.
 */
struct dfTablePair Test_NumericMomentumSpaceDistribution::Generate(slibreal_t E, slibreal_t lnL, slibreal_t Z) {
    struct dfTablePair retval;
    unsigned int i, j;
    slibreal_t *p, *xi, *f,
               *rtest, *ptest, *xitest, *ftest;
    AnalyticalAvalanche *af = new AnalyticalAvalanche();
    af->Initialize(E, lnL, Z);
    
    srand(time(NULL));

    // Generate p, xi and f
    p  = new slibreal_t[NP];
    xi = new slibreal_t[NXI];
    f  = new slibreal_t[NP*NXI];

    for (i = 0; i < NP; i++)
        p[i] = PMIN + (PMAX-PMIN) * ((slibreal_t)i) / ((slibreal_t)NP-1.0);

    for (j = 0; j < NXI; j++) {
        xi[j] = -1.0 + 2.0*((slibreal_t)j) / ((slibreal_t)NXI-1.0);
        for (i = 0; i < NP; i++)
            f[i+j*NP] = af->Eval(p[i], xi[j]);
    }

    // Create the momentum-space distribution
    NumericMomentumSpaceDistributionFunction *msdf = new NumericMomentumSpaceDistributionFunction();
    msdf->InitializeLog(NP, NXI, p, xi, f);

    // Generate test points
    rtest  = new slibreal_t[NTESTPOINTS];
    ptest  = new slibreal_t[NTESTPOINTS];
    xitest = new slibreal_t[NTESTPOINTS];
    ftest  = new slibreal_t[NTESTPOINTS];
    for (i = 0; i < NTESTPOINTS; i++) {
        rtest[i]  = 0.0;
        ptest[i]  = PMIN + (PMAX-PMIN) * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX-1.0);
        xitest[i] = -1.0 + 2.0*((slibreal_t)rand()) / ((slibreal_t)RAND_MAX-1.0);

        ftest[i] = af->Eval(ptest[i], xitest[i]);
    }
    
    retval.df = msdf;
    retval.r  = rtest;
    retval.p  = ptest;
    retval.xi = xitest;
    retval.f  = ftest;
    retval.n  = NTESTPOINTS;

    return retval;
}

/**
 * Run the test.
 */
bool Test_NumericMomentumSpaceDistribution::Run() {
    struct dfTablePair dftp;
    unsigned int failed;

    dftp = Generate(EHat, lnLambda, Zeff);
    failed = CompareToTable(dftp.df, dftp.r, dftp.p, dftp.xi, dftp.f, dftp.n, TOLERANCE);

    delete [] dftp.r;
    delete [] dftp.p;
    delete [] dftp.xi;
    delete [] dftp.f;
    delete dftp.df;

    if (failed > 0)
        this->PrintWarning("--> Numeric momentum-space distribution: %u of %u (%.2f%%) tests had an error larger than %.2f%%.", failed, dftp.n, ((slibreal_t)failed)/((slibreal_t)dftp.n)*100.0, TOLERANCE*100.0);
    
    if (failed > dftp.n*ALLOWED_FAILURES) {
        this->PrintError("--> Too many failures for numeric momentum-space distribution.");
        return false;
    } else
        return true;
}

