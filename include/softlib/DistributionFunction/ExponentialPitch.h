#ifndef _EXPONENTIAL_PITCH_H
#define _EXPONENTIAL_PITCH_H

#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>

class ExponentialPitch : public MomentumSpaceDistributionFunction {
    protected:
        slibreal_t C0, C1;
    public:
        ExponentialPitch();
        ExponentialPitch(const slibreal_t, const slibreal_t C1=0);
        using MomentumSpaceDistributionFunction::Eval;
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) override;
        void Initialize(const slibreal_t, const slibreal_t C1=0);
        virtual ExponentialPitch *MinClone() override;
};

#endif/*_EXPONENTIAL_PITCH_H*/
