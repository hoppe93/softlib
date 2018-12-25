#ifndef _MAGNETIC_FIELD_ANALYTICAL2D_H
#define _MAGNETIC_FIELD_ANALYTICAL2D_H

#include <string>
#include <softlib/config.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>

enum MFASafetyFactorType {
	MFASF_CONSTANT,
	MFASF_LINEAR,
	MFASF_QUADRATIC,
	MFASF_EXPONENTIAL
};

enum MFAToroidalFieldSign {
    MFATFS_CW,      // Toroidal field pointing clock-wise (seen from above tokamak)
    MFATFS_CCW      // Toroidal field pointing counter-clock-wise (seen from above tokamak)
};

#define MAGNETIC_FIELD_ANALYTICAL2D_NDOMAINPOINTS 50

class MagneticFieldAnalytical2D : public MagneticField2D {
	private:
		slibreal_t retval[3];
		slibreal_t retval_magnitude;
        struct flux_diff flux_retval;

		/* Safety factor fields */
		slibreal_t safety_factor_param1,
				   safety_factor_param2;
		enum MFASafetyFactorType safety_factor_type;

		slibreal_t B0, Rm, rminor, sigmaB, **jacobian;
	public:
		MagneticFieldAnalytical2D(slibreal_t, slibreal_t, slibreal_t, enum MFAToroidalFieldSign, enum MFASafetyFactorType, slibreal_t, slibreal_t);
		MagneticFieldAnalytical2D(slibreal_t, slibreal_t, slibreal_t, enum MFAToroidalFieldSign, enum MFASafetyFactorType, slibreal_t, slibreal_t, const std::string&, const std::string&);
		virtual ~MagneticFieldAnalytical2D() {}

        MagneticFieldAnalytical2D *Clone();

        void ConstructDomain(const slibreal_t, const slibreal_t);
		slibreal_t *Eval(slibreal_t*);
		slibreal_t *Eval(slibreal_t, slibreal_t, slibreal_t);
		struct magnetic_field_data& EvalDerivatives(slibreal_t*);
		struct magnetic_field_data& EvalDerivatives(slibreal_t, slibreal_t, slibreal_t);
        virtual slibreal_t EvalFlux(slibreal_t, slibreal_t, slibreal_t) override;
        virtual struct flux_diff *EvalFluxDerivatives(slibreal_t, slibreal_t, slibreal_t) override;
        MagneticFieldNumeric2D *ToNumeric2D(const unsigned int nr=120, const unsigned int nz=120);

        slibreal_t GetRMajor() { return this->Rm; }
        slibreal_t GetRMinor() { return this->rminor; }
        virtual bool HasMagneticFlux() override { return true; }

		slibreal_t __GetrDqDr(slibreal_t);
		slibreal_t __GetSafetyFactor(slibreal_t);

        slibreal_t GetB0() { return this->B0; }
        void SetB0(const slibreal_t B0) { this->B0 = B0; }
};

#endif/*_MAGNETIC_FIELD_ANALYTICAL2D_H*/
