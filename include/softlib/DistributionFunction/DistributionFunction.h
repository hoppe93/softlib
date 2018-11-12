#ifndef _DISTRIBUTION_FUNCTION_H
#define _DISTRIBUTION_FUNCTION_H

#include <softlib/config.h>

class DistributionFunction {
    public:
        virtual ~DistributionFunction() {}
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t drift_shift=0.0) = 0;
        virtual DistributionFunction *MinClone() = 0;
};

#endif/*_DISTRIBUTION_FUNCTION_H*/
