#ifndef _UNIFORM_RADIAL_PROFILE_H
#define _UNIFORM_RADIAL_PROFILE_H

#include <softlib/config.h>
#include <softlib/DistributionFunction/RadialProfile.h>

class UniformRadialProfile : public RadialProfile {
    public:
        slibreal_t Eval(const slibreal_t rho __attribute__((unused)), const slibreal_t drift_shift __attribute__((unused))) { return 1.0; }
        UniformRadialProfile *MinClone() { return new UniformRadialProfile(); }
};

#endif/*_UNIFORM_RADIAL_PROFILE_H*/
