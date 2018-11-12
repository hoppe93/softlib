#ifndef _LINEAR_RADIAL_PROFILE_H
#define _LINEAR_RADIAL_PROFILE_H

#include <softlib/config.h>
#include <softlib/DistributionFunction/RadialProfile.h>

class LinearRadialProfile : public RadialProfile {
    private:
        // Bounds on radial interval on which this
        // profile is defined. Outside the interval,
        // this profile is identically zero.
        slibreal_t rmin, rmax, delr;
    public:
        LinearRadialProfile(const slibreal_t, const slibreal_t);
        slibreal_t Eval(const slibreal_t, const slibreal_t drift_shift=0.0);
        LinearRadialProfile *MinClone();
};

#endif/*_LINEAR_RADIAL_PROFILE_H*/
