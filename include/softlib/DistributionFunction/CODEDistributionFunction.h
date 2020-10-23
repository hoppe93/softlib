#ifndef _CODE_DISTRIBUTION_FUNCTION_H
#define _CODE_DISTRIBUTION_FUNCTION_H

#include <string>
#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>
#include <softlib/SFile.h>

#ifdef INTERP_SPLINTER
#   error "SPLINTER interpolation has not been implemented in the CODEDistributionFunction class yet."
#else
#   include <gsl/gsl_interp.h>
#endif

class CODEDistributionFunction : public MomentumSpaceDistributionFunction {
    private:
        int interptype=INTERPOLATION_LINEAR;
        unsigned int timestep, np, nleg;

        slibreal_t *f=nullptr, *p=nullptr;
        slibreal_t *legendre=nullptr, xileg;
        unsigned int lmodes=0;

#ifdef INTERP_SPLINTER
#else
        gsl_interp_accel *interp_accel=nullptr;
        gsl_interp **interp=nullptr;
#endif
    public:
        CODEDistributionFunction();
        CODEDistributionFunction(const std::string&, int time=-1, int interptype=INTERPOLATION_LINEAR, const std::string& path="");
        CODEDistributionFunction(SFile*, int time=-1, int interptype=INTERPOLATION_LINEAR, const std::string& path="");
        ~CODEDistributionFunction();

        void ComputeLegendrePolynomials(slibreal_t*, slibreal_t, int);
        void Destroy();
        using MomentumSpaceDistributionFunction::Eval;
        virtual slibreal_t Eval(const slibreal_t, const slibreal_t) override;
        void Load(const std::string&, int time=-1, int interptype=INTERPOLATION_LINEAR, const std::string& path="");
        void Load(SFile*, int time=-1, int interptype=INTERPOLATION_LINEAR, const std::string& path="");
        void Initialize(const unsigned int, const unsigned int, slibreal_t*, slibreal_t*, int interptype=INTERPOLATION_LINEAR);
        void InitInterpolation(int interptype=INTERPOLATION_LINEAR);
        virtual CODEDistributionFunction *MinClone() override;
        NumericMomentumSpaceDistributionFunction *ToMomentumSpace(
            bool logarithmic=true, bool uniformTheta=false, unsigned int nxi=200,
            int interptype=NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR
        );

        slibreal_t GetPMax();
        slibreal_t GetPMin();

        enum {
            INTERPOLATION_LINEAR,
            INTERPOLATION_POLYNOMIAL,
            INTERPOLATION_CSPLINE,
            INTERPOLATION_CSPLINE_PERIODIC,
            INTERPOLATION_AKIMA,
            INTERPOLATION_AKIMA_PERIODIC,
            INTERPOLATION_STEFFEN
        };
};

#endif/*_CODE_DISTRIBUTION_FUNCTION_H*/
