#ifndef _SOFT_DISTRIBUTION_FUNCTION_H
#define _SOFT_DISTRIBUTION_FUNCTION_H

#include <string>
#include <softlib/DistributionFunction/NumericDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SFile.h>

#define ELECTRON_MASS_EV        5.109989461e5       // Electron mass in units of eV/mc^2

class SOFTDistributionFunction : public NumericDistributionFunction {
    public:
        static const std::string MAGIC;

        struct softmdf_data {
            unsigned int np, nxi;
            slibreal_t *p, *xi;
            slibreal_t *f;
        };
        struct softdf_data {
            unsigned int nr;
            slibreal_t *r;
            struct softmdf_data *fr;
        };

        SOFTDistributionFunction(const std::string&, MagneticField2D *mf=nullptr, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);

        void Load(const std::string&, MagneticField2D*, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);
        static struct softdf_data *__Load(const std::string&, MagneticField2D *mf=nullptr);
        static void __LoadRadius(SFile*, const std::string, struct softmdf_data*);
        static NumericMomentumSpaceDistributionFunction *LoadMomentumSpace(
            const std::string&, MagneticField2D *mf=nullptr,
			bool logarithmic=true, unsigned int radindex=0,
            int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC
        );
        static bool Verify(const std::string&, double**, sfilesize_t, sfilesize_t, SFile*);
};

#endif/*_SOFT_DISTRIBUTION_FUNCTION_H*/
