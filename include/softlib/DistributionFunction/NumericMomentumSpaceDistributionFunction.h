#ifndef _NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H
#define _NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H

#include <cstdlib>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <softlib/config.h>

#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>

#ifdef INTERP_SPLINTER
#   error "SPLINTER interpolation has not been implemented in the NumericMomentumSpaceDistributionFunction class yet."
#else
#   include <gsl/gsl_interp2d.h>
#   include <gsl/gsl_spline2d.h>
#endif

class NumericMomentumSpaceDistributionFunction : public MomentumSpaceDistributionFunction {
    private:
        slibreal_t *f=NULL,           // Distribution function (of size np*nxi)
                   *p=NULL, *xi=NULL; // Momentum-space grid points; p = magnitude of momentum, xi = cosine of pitch angle
        unsigned int np, nxi;         // Length of p and xi respectively

        /* Interpolation properties */
#ifdef INTERP_SPLINTER
#else   /* INTERP_GSL */
        gsl_interp_accel *pa, *xia;
        gsl_spline2d *fspline;
#endif

        slibreal_t pmin, pmax;
        slibreal_t ximin, ximax;

        bool flipPitchSign = false;
        bool logarithmic = false;
        int interptype = INTERPOLATION_LINEAR;
    public:
        using MomentumSpaceDistributionFunction::Eval;
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) override;
        void Initialize(
            const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*,
            int interp=INTERPOLATION_LINEAR
        );
        void InitializeLog(
            const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*,
            int interp=INTERPOLATION_LINEAR
        );
        bool IsLogarithmic() { return this->logarithmic; }
        virtual NumericMomentumSpaceDistributionFunction *MinClone() override;
        void SetLogarithmic(bool l) { this->logarithmic = l; }

        void FlipPitchSign(bool v=true) { this->flipPitchSign = v; }

        // Available interpolation methods
        enum {
            INTERPOLATION_LINEAR,
            INTERPOLATION_CUBIC
        };
};

#endif/*_NUMERIC_MOMENTUM_SPACE_DISTRIBUTION_FUNCTION_H*/
