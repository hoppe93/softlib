#ifndef _LUKE_DISTRIBUTION_FUNCTION_H
#define _LUKE_DISTRIBUTION_FUNCTION_H

#include <string>
#include <softlib/DistributionFunction/NumericDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SFile.h>

#define ELECTRON_MASS_EV        5.109989461e5       // Electron mass in units of eV/mc^2

class LUKEDistributionFunction : public NumericDistributionFunction {
	private:
		MagneticField2D *magfield;
    public:
		struct lukedf_data {
			unsigned int nr, np, nxi;
			slibreal_t *r, *p, *xi, *f;
		};

        LUKEDistributionFunction(const std::string&, MagneticField2D*, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);

        void Load(const std::string&, bool logarithmic=true, int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC);
        static struct lukedf_data* __Load(const std::string&, MagneticField2D*);
        static NumericMomentumSpaceDistributionFunction *LoadMomentumSpace(
            const std::string&, MagneticField2D*,
			bool logarithmic=true, unsigned int radindex=0,
            int interptype=NumericDistributionFunction::INTERPOLATION_CUBIC
        );
};

#endif/*_LUKE_DISTRIBUTION_FUNCTION_H*/
