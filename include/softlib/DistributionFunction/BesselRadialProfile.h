#ifndef _BESSEL_RADIAL_PROFILE_H
#define _BESSEL_RADIAL_PROFILE_H

#include <softlib/config.h>
#include <softlib/DistributionFunction/RadialProfile.h>

class BesselRadialProfile : public RadialProfile {
    private:
        // Bounds on radial interval on which this
        // profile is defined. Outside the interval,
        // this profile is identically zero.
        slibreal_t rmin, rmax;
        
        // First zero of the 0th Bessel function of the first kind J0
        const slibreal_t J0_ZERO1 = 2.40482555769577;
    public:
        BesselRadialProfile(const slibreal_t, const slibreal_t);
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t drift_shift=0.0) override;
        virtual BesselRadialProfile *MinClone() override;
};

#endif/*_BESSEL_RADIAL_PROFILE_H*/
