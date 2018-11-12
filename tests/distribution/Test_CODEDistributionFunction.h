#ifndef _TEST_CODE_DISTRIBUTION_FUNCTION_H
#define _TEST_CODE_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/CODEDistributionFunction.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include "Test_DistributionFunction.h"
#include "runtest.h"

class Test_CODEDistributionFunction : public Test_DistributionFunction {
    private:
        const string distfilename = "../tests/distribution/code_test_distribution.mat";
        const string corrfilename = "../tests/distribution/code_corr_values.mat";
        
        const unsigned int
            NPOINTS         = 300;
        const slibreal_t
            TOLERANCE1      = 0.1,
            TOLERANCE2      = 0.2,
            ALLOWED_FAILURES= 0.1;

    public:
        Test_CODEDistributionFunction(const string& name) : Test_DistributionFunction(name) {}

        bool CreateNumericAndTest(CODEDistributionFunction*);
        bool LoadAndTestCorrect(const string&, CODEDistributionFunction*);
        bool Run();
};

#endif/*_TEST_CODE_DISTRIBUTION_FUNCTION_H*/
