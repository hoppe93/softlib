#ifndef _ANALYTIC_AVALANCHE_H
#define _ANALYTIC_AVALANCHE_H

#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>

class AnalyticalAvalanche : public MomentumSpaceDistributionFunction {
    protected:
        slibreal_t EHat, lnLambda, Zeff;
        slibreal_t Ac, gmm0;
    public:
        AnalyticalAvalanche();
        AnalyticalAvalanche(const slibreal_t, const slibreal_t, const slibreal_t);
        slibreal_t Eval(const slibreal_t, const slibreal_t);
        void Initialize(const slibreal_t, const slibreal_t, const slibreal_t);
        AnalyticalAvalanche *MinClone();
};

#endif/*_ANALYTIC_AVALANCHE_H*/
