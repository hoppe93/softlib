#ifndef _GAUSSIAN_RADIAL_PROFILE_H
#define _GAUSSIAN_RADIAL_PROFILE_H

#include <softlib/config.h>
#include <softlib/DistributionFunction/RadialProfile.h>

class GaussianRadialProfile : public RadialProfile {
    private:
      // Bounds on radial interval on which this
      // profile is defined. Outside the interval,
      // this profile is identically zero.
        slibreal_t rmin, rmax, delr;
        slibreal_t a,b,c;// Gaussian Parameters
    public:
        GaussianRadialProfile(const slibreal_t, const slibreal_t,const slibreal_t, const slibreal_t, const slibreal_t);
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t drift_shift=0.0) override;
        virtual GaussianRadialProfile *MinClone() override;
};

#endif/*_GAUSSIAN_RADIAL_PROFILE_H*/
