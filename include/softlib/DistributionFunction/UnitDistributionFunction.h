#ifndef _UNIT_DISTRIBUTION_FUNCTION_H
#define _UNIT_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>

class UnitDistributionFunction : public MomentumSpaceDistributionFunction {
	public:
        using MomentumSpaceDistributionFunction::Eval;
		virtual slibreal_t Eval(
            const slibreal_t p __attribute__((unused)),
            const slibreal_t xi __attribute__((unused))
        ) override { return 1.0; }
        virtual UnitDistributionFunction *MinClone() override;
};

#endif/*_UNIT_DISTRIBUTION_FUNCTION_H*/
