#ifndef _RADIAL_PROFILE_H
#define _RADIAL_PROFILE_H

class RadialProfile {
    public:
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) = 0;
        virtual RadialProfile *MinClone() = 0;
};

#endif/*_RADIAL_PROFILE_H*/
