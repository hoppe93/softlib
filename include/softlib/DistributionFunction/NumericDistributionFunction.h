#ifndef _NUMERIC_DISTRIBUTION_FUNCTION_H
#define _NUMERIC_DISTRIBUTION_FUNCTION_H

#include <cstdlib>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <softlib/config.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>

#ifdef INTERP_SPLINTER
#   error "SPLINTER interpolation has not been implemented in the NumericDistributionFunction class yet."
#else
#   include <gsl/gsl_interp2d.h>
#   include <gsl/gsl_spline2d.h>
#endif

class NumericDistributionFunction : public DistributionFunction {
    private:
        slibreal_t **f=nullptr,     // Distribution function (of size nr*np*nxi)
                   // Grid points
                   *r=nullptr,      // Radial coordinate
                   *p=nullptr,      // Magnitude of momentum
                   *xi=nullptr;     // Cosine of pitch angle
        unsigned int nr, np, nxi;   // Length of r, p and xi respectively

        /**
         * For simplicity, we consider the 3-D distribution
         * function as a set of 'nr' momentum-space
         * distribution functions f_r(p,xi). Since the radial
         * distribution function generally varies as a polynomial
         * (in contrast to the momentum-space part, which
         * generally varies exponentially) we choose to do
         * proper interpolation in the momentum-space part
         * and interpolate linearly in radius. We therefore
         * have that
         *
         *   f(R0,p,xi) = (1-dr)*f_{r}(p,xi) + dr*f_{r+1}(p,xi)
         *
         * where r is the radial index corresponding to
         * R(r) < R0 < R(r+1), and dr = R0-R(r).
         */
        unsigned int nmsdf;
        NumericMomentumSpaceDistributionFunction *msdf;

        slibreal_t rmin, rmax;
        slibreal_t pmin, pmax;
        slibreal_t ximin, ximax;

        bool logarithmic = false;
        int interptype = INTERPOLATION_CUBIC;       // Type of interpolation done
    public:
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t drift_shift=0.0);
        int __FindNearestR(const slibreal_t);
        void Initialize(
            const unsigned int, const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t**,
            int interp=INTERPOLATION_CUBIC
        );
        void InitializeLog(
            const unsigned int, const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t**,
            int interp=INTERPOLATION_CUBIC, bool alloc=false
        );
        bool IsLogarithmic() { return this->logarithmic; }
        virtual NumericDistributionFunction *MinClone();
        void SetLogarithmic(bool l) { this->logarithmic = l; }

        // Available interpolation methods
        enum {
            INTERPOLATION_LINEAR,
            INTERPOLATION_CUBIC
        };
};

#endif/*_NUMERIC_DISTRIBUTION_FUNCTION_H*/
