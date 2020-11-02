#ifndef _DREAM_DISTRIBUTION_FUNCTION_H
#define _DREAM_DISTRIBUTION_FUNCTION_H

#include <string>
#include "softlib/config.h"
#include "softlib/DistributionFunction/DistributionFunction.h"
#include "softlib/DistributionFunction/NumericDistributionFunction.h"
#include "softlib/MagneticField/MagneticField2D.h"

class DREAMDistributionFunction : public NumericDistributionFunction {
private:
    struct dreamdf_data {
        std::string distname;
        unsigned int nr, nxi, np;
        slibreal_t *r, *xi, *p;
        slibreal_t *f;
    };
    static struct dreamdf_data *__Load(const std::string&, MagneticField2D *mf=nullptr, const std::string& distname="");

    struct dreamdf_data *rawdata;

public:
    DREAMDistributionFunction(const std::string&, MagneticField2D *mf=nullptr, const std::string& distname="", bool logarithmic=false, int interptype=NumericDistributionFunction::INTERPOLATION_LINEAR);
    ~DREAMDistributionFunction();

    void Load(const std::string&, MagneticField2D*, const std::string& distname="", bool logarithmic=false, int interptype=NumericDistributionFunction::INTERPOLATION_LINEAR);
};

#endif/*_DREAM_DISTRIBUTION_FUNCTION_H*/
