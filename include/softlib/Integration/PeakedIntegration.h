#ifndef _SOFTLIB_PEAKED_INTEGRATION_H
#define _SOFTLIB_PEAKED_INTEGRATION_H

#include <functional>
#include <softlib/config.h>

class PeakedIntegration {
    private:
    public:
        PeakedIntegration();

        slibreal_t Evaluate(
            std::function<slibreal_t(slibreal_t)>&,
            const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t);
        slibreal_t Evaluate(
            std::function<slibreal_t(slibreal_t)>&,
            const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t
        );
        slibreal_t FindMaximum(
            std::function<slibreal_t(slibreal_t)>&, const slibreal_t tol=1e-2
        );
};

#endif/*_SOFTLIB_PEAKED_INTEGRATION_H*/
