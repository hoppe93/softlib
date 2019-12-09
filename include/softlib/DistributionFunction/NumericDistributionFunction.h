#ifndef _NUMERIC_DISTRIBUTION_FUNCTION_H
#define _NUMERIC_DISTRIBUTION_FUNCTION_H

#include <cstdlib>
#include <vector>
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
        std::vector<slibreal_t> r;      // Radial grid

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
         * R(r) < R0 < R(r+1), and dr = (R0-R(r)) / (R(r+1)-R(r)).
         */
        //unsigned int nmsdf;
        //NumericMomentumSpaceDistributionFunction *msdf;
        std::vector<NumericMomentumSpaceDistributionFunction*> msdf;

        slibreal_t rmin, rmax;

		bool allowExtrapolation = true;
        bool logarithmic = false;
        int interptype = INTERPOLATION_CUBIC;       // Type of interpolation done
    public:
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t drift_shift=0.0) override;
        int __FindNearestR(const slibreal_t);
        void Initialize(
            const unsigned int, const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t*,
            int interp=INTERPOLATION_CUBIC
        );
        void InitializeLog(
            const unsigned int, const unsigned int, const unsigned int,
            slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t*,
            int interp=INTERPOLATION_CUBIC, bool alloc=false
        );
        int InsertMomentumSpaceDistribution(
            slibreal_t, NumericMomentumSpaceDistributionFunction*
        );
        int InsertMomentumSpaceDistribution(
            const unsigned int, const unsigned int, slibreal_t,
            slibreal_t*, slibreal_t*, slibreal_t*, int
        );
        int InsertMomentumSpaceDistributionLog(
            const unsigned int, const unsigned int, slibreal_t,
            slibreal_t*, slibreal_t*, slibreal_t*, int, bool alloc=false
        );
        bool IsLogarithmic() { return this->logarithmic; }
        virtual NumericDistributionFunction *MinClone() override;
		void SetExtrapolationPermission(bool e) { this->allowExtrapolation = e; }
        void SetLogarithmic(bool l) { this->logarithmic = l; }

        void FlipPitchSign(bool v) {
            std::vector<NumericMomentumSpaceDistributionFunction*>::iterator it;
            for (it = this->msdf.begin(); it != this->msdf.end(); it++)
                (*it)->FlipPitchSign(v);
        }

        // Available interpolation methods
        enum {
            INTERPOLATION_LINEAR,
            INTERPOLATION_CUBIC
        };
};

#endif/*_NUMERIC_DISTRIBUTION_FUNCTION_H*/
