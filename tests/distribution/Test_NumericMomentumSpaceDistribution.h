#ifndef _TEST_NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_H
#define _TEST_NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_H

#include "runtest.h"
#include "Test_DistributionFunction.h"

class Test_NumericMomentumSpaceDistribution : public Test_DistributionFunction {
    private:
        const unsigned int
            NP              = 200,  // Number of points in p grid
            NXI             = 100,  // Number of points in xi grid
            NTESTPOINTS     = 300;  // Number of test points to evaluate both distributions in
        
        // Input parameters to analytical avalanche distribution
        const slibreal_t
            EHat        = 2.0,
            lnLambda    = 17.0,
            Zeff        = 3.0,

            PMAX        = 100.0,
            PMIN        = 0.1,      // p threshold to avoid f(p=0) = inf
            TOLERANCE   = 0.05,

            // Since the test points are random (which is
            // good to avoid testing interpolation in only
            // well-behaved points) we allow a certain percentage
            // of points to be "failures", i.e. to have an error
            // larger than TOLERANCE.
            ALLOWED_FAILURES = 0.1;
    public:
        Test_NumericMomentumSpaceDistribution(const string& name) : Test_DistributionFunction(name) {}
        struct dfTablePair Generate(slibreal_t, slibreal_t, slibreal_t);
        bool Run();
};

#endif/*_TEST_NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_H*/
