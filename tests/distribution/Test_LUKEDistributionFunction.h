#ifndef _TEST_LUKE_DISTRIBUTION_FUNCTION_H
#define _TEST_LUKE_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include "Test_NumericDistributionFunction.h"
#include "runtest.h"

class Test_LUKEDistributionFunction : public Test_NumericDistributionFunction {
    private:
        const string inputfilename = "../tests/distribution/luke_test_distribution.mat";
    public:
        Test_LUKEDistributionFunction(const string& name) : Test_NumericDistributionFunction(name) {}
        bool Run();
};

#endif/*_TEST_LUKE_DISTRIBUTION_FUNCTION_H*/
