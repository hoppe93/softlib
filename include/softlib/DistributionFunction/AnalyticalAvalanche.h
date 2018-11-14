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
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) override;
        void Initialize(const slibreal_t, const slibreal_t, const slibreal_t);
        virtual AnalyticalAvalanche *MinClone() override;
};

#endif/*_ANALYTIC_AVALANCHE_H*/
