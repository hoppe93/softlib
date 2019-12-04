#ifndef _GOCODE_DISTRIBUTION_FUNCTION_H
#define _GOCODE_DISTRIBUTION_FUNCTION_H

#include <string>
#include <softlib/DistributionFunction/CODEDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>


class GOCODEDistributionFunction : public DistributionFunction {
    private:
        unsigned int nr;
        slibreal_t *rho;
        CODEDistributionFunction **code;

        slibreal_t rhomin, rhomax;

        bool allowExtrapolation = true;

        int __FindNearestR(const slibreal_t);

    public:
        static const std::string MAGIC;

        GOCODEDistributionFunction() {};
        GOCODEDistributionFunction(const std::string&, MagneticField2D*, int, int interptype=CODEDistributionFunction::INTERPOLATION_LINEAR);
        ~GOCODEDistributionFunction();

        void Initialize(unsigned int, slibreal_t*, CODEDistributionFunction**);
        void Load(const std::string&, MagneticField2D*, int, int interptype=CODEDistributionFunction::INTERPOLATION_LINEAR);

        virtual slibreal_t Eval(const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t drift_shift=0.0) override;
        virtual GOCODEDistributionFunction *MinClone() override;
};

#endif/*_GOCODE_DISTRIBUTION_FUNCTION_H*/
