#ifndef _TEST_SOFT_DISTRIBUTION_FUNCTION_H
#define _TEST_SOFT_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include "Test_DistributionFunction.h"
#include "runtest.h"

class Test_SOFTDistributionFunction : public Test_DistributionFunction {
    private:
        const string inputfilename = "../tests/distribution/soft_test_distribution.mat";
        
        const unsigned int
            NPOINTS1        = 300,
            NPOINTS2        = 300;
        const slibreal_t
            EHAT            = 4.0,
            LOGLAMBDA       = 17.0,
            ZEFF            = 2.0,

            RMIN            = 1.6,
            RMAX            = 2.2,
            PMIN            = 0.1,
            PMAX            = 100.0,
            TOLERANCE1      = 0.2,
            TOLERANCE2      = 0.05,
            ALLOWED_FAILURES= 0.1;

    public:
        Test_SOFTDistributionFunction(const string& name) : Test_DistributionFunction(name) {}
        RadialDistributionFunction *GeneratePhaseSpaceDistribution(
            slibreal_t, slibreal_t, slibreal_t,
            slibreal_t, slibreal_t
        );
        unsigned int RandomPoints(
            DistributionFunction*, DistributionFunction*,
            const unsigned int, const slibreal_t,
            const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t
        );
        bool Run();
};

#endif/*_TEST_SOFT_DISTRIBUTION_FUNCTION_H*/
