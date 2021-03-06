#ifndef _RADIAL_DISTRIBUTION_FUNCTION_H
#define _RADIAL_DISTRIBUTION_FUNCTION_H

class RadialDistributionFunction;

#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>
#include <softlib/DistributionFunction/RadialProfile.h>

class RadialDistributionFunction : public DistributionFunction {
    protected:
        MomentumSpaceDistributionFunction *msdf;
        RadialProfile *rp;
    public:
        RadialDistributionFunction();
        RadialDistributionFunction(RadialProfile*, MomentumSpaceDistributionFunction*);

        virtual slibreal_t Eval(const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t drift_shift=0.0) override;
        void Initialize(RadialProfile*, MomentumSpaceDistributionFunction*);
        virtual RadialDistributionFunction *MinClone() override;
};

#endif/*_RADIAL_DISTRIBUTION_FUNCTION_H*/
