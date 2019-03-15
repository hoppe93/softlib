#ifndef _TEST_SOFT_DISTRIBUTION_FUNCTION_H
#define _TEST_SOFT_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include "Test_NumericDistributionFunction.h"
#include "runtest.h"

class Test_SOFTDistributionFunction : public Test_NumericDistributionFunction {
    private:
        const string inputfilename = "../tests/distribution/soft_test_distribution.mat";
    public:
        Test_SOFTDistributionFunction(const string& name) : Test_NumericDistributionFunction(name) {}
        bool Run();
};

#endif/*_TEST_SOFT_DISTRIBUTION_FUNCTION_H*/
