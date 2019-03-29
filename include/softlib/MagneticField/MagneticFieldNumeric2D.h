#ifndef _MAGNETIC_FIELD_NUMERIC_2D_H
#define _MAGNETIC_FIELD_NUMERIC_2D_H

#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SFile.h>

#ifdef INTERP_SPLINTER
#	error "SPLINTER interpolation has not been implemented in the MagneticFieldNumeric2D class yet."
#else
#	include <gsl/gsl_interp2d.h>
#	include <gsl/gsl_spline2d.h>
#endif

class MagneticFieldNumeric2D : public MagneticField2D {
	protected:
		slibreal_t *Br=nullptr, *Bphi=nullptr, *Bz=nullptr;
        slibreal_t *Psi=nullptr;
		slibreal_t *R=nullptr, *Z=nullptr;
		slibreal_t *rsep=nullptr, *zsep=nullptr;
		slibreal_t *rwall=nullptr, *zwall=nullptr;
		unsigned int nr, nz, nwall, nsep;

        /* Interpolated value storage */
        slibreal_t *interpval;
        slibreal_t **jacobian;
        struct flux_diff flux_data;

		slibreal_t rmin, rmax, zmin, zmax;

		slibreal_t **__EvalJacobian(slibreal_t, slibreal_t, slibreal_t);

        bool hasFluxCoordinates = false;

		/* Interpolation properties */
#	ifdef INTERP_SPLINTER
#	else	/* INTERP_GSL */
		gsl_interp_accel *ra, *za;
		gsl_spline2d *sBr, *sBphi, *sBz;
        gsl_spline2d *sPsi;
#	endif
	public:
		MagneticFieldNumeric2D();
		MagneticFieldNumeric2D(const std::string&);
		MagneticFieldNumeric2D(const std::string&, enum sfile_type);
		MagneticFieldNumeric2D(
			const std::string& name, const std::string& description,
			slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
			slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz, slibreal_t *Psi,
            slibreal_t raxis, slibreal_t zaxis,
			slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
			slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
		);
		virtual ~MagneticFieldNumeric2D();

        MagneticFieldNumeric2D *Clone();

		void BaseInit();
		slibreal_t *Eval(slibreal_t*);
		slibreal_t *Eval(slibreal_t, slibreal_t, slibreal_t);
		struct magnetic_field_data& EvalDerivatives(slibreal_t*);
		struct magnetic_field_data& EvalDerivatives(slibreal_t, slibreal_t, slibreal_t);
        virtual slibreal_t EvalFlux(slibreal_t, slibreal_t, slibreal_t) override;
        virtual struct flux_diff *EvalFluxDerivatives(slibreal_t, slibreal_t, slibreal_t) override;
		slibreal_t FindMaxRadius() override;
		slibreal_t FindMinRadius() override;
		std::string& GetDescription();
		std::string& GetName();
        slibreal_t *GetMagneticAxis();
		void Init(
			const std::string& name, const std::string& description,
			slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
			slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz, slibreal_t *Psi,
            slibreal_t raxis, slibreal_t zaxis,
			slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
			slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
		);
		void InitInterpolation();
		virtual void Load(const std::string&);
		virtual void Load(const std::string&, enum sfile_type);
		void Save(const std::string&);
		void Save(const std::string&, enum sfile_type);
		/* Transpose a matrix from SFile */
		double **Transpose(double**, sfilesize_t, sfilesize_t);

        virtual bool HasMagneticFlux() override { return hasFluxCoordinates; }
};

#endif/*_MAGNETIC_FIELD_NUMERIC_2D_H*/
