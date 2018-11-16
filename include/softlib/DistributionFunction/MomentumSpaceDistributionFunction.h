#ifndef _MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H
#define _MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H

#include <cstdlib>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <softlib/config.h>

class MomentumSpaceDistributionFunction;

#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/RadialProfile.h>
#include <softlib/DistributionFunction/UniformRadialProfile.h>

class MomentumSpaceDistributionFunction : public DistributionFunction {
    public:
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) = 0;
        virtual slibreal_t Eval(
            const slibreal_t rho __attribute__((unused)), const slibreal_t p,
            const slibreal_t xi, const slibreal_t drift_shift __attribute__((unused))
        ) override { return this->Eval(p, xi); }

        virtual MomentumSpaceDistributionFunction *MinClone() = 0;

        RadialDistributionFunction *ToRadialDistribution();
        RadialDistributionFunction *ToRadialDistribution(RadialProfile*);
};

#endif/*_MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H*/
