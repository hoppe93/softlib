#ifndef _TEST_DISTRIBUTION_FUNCTION_H
#define _TEST_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/DistributionFunction.h>
#include "runtest.h"

struct dfTablePair {
    DistributionFunction *df;
    slibreal_t *r, *p, *xi, *f;
    unsigned int n;
};

class Test_DistributionFunction : public UnitTest {
    public:
        Test_DistributionFunction(const string &name) : UnitTest(name) {}
        unsigned int CompareToTable(
            DistributionFunction*, const slibreal_t*,
            const slibreal_t*, const slibreal_t*,
            const slibreal_t*, const unsigned int,
            const slibreal_t
        );
};

#endif/*_TEST_DISTRIBUTION_FUNCTION_H*/
