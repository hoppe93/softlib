#ifndef _DEPR_SOFT_DISTRIBUTION_FUNCTION_H
#define _DEPR_SOFT_DISTRIBUTION_FUNCTION_H

#include <string>
#include <softlib/DistributionFunction/NumericDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SFile.h>

#define ELECTRON_MASS_EV        5.109989461e5       // Electron mass in units of eV/mc^2

struct deprsoftdf_data {
    unsigned int nr, np, nxi;
    slibreal_t *r, *p, *xi, *f;
};

class deprSOFTDistributionFunction : public NumericDistributionFunction {
    public:
        deprSOFTDistributionFunction(const std::string&, MagneticField2D *mf=nullptr, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);

        void Load(const std::string&, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);
        static struct deprsoftdf_data* __Load(const std::string&, MagneticField2D *mf=nullptr);
        static NumericMomentumSpaceDistributionFunction *LoadMomentumSpace(
            const std::string&, MagneticField2D *mf=nullptr,
			bool logarithmic=true, unsigned int radindex=0,
            int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC
        );
        static bool Verify(const std::string&, double**, sfilesize_t, sfilesize_t, SFile*);
};

#endif/*_DEPR_SOFT_DISTRIBUTION_FUNCTION_H*/
