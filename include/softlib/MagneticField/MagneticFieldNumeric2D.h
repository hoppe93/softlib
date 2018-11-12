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
	private:
		slibreal_t *Br, *Bphi, *Bz;
		slibreal_t *R, *Z;
		slibreal_t *rsep, *zsep;
		slibreal_t *rwall, *zwall;
		unsigned int nr, nz, nwall, nsep;

        /* Interpolated value storage */
        slibreal_t *interpval;
        slibreal_t **jacobian;

		slibreal_t rmin, rmax, zmin, zmax;

		slibreal_t **__EvalJacobian(slibreal_t, slibreal_t, slibreal_t);

		/* Interpolation properties */
#	ifdef INTERP_SPLINTER
#	else	/* INTERP_GSL */
		gsl_interp_accel *ra, *za;
		gsl_spline2d *sBr, *sBphi, *sBz;
#	endif
	public:
		MagneticFieldNumeric2D(const std::string&);
		MagneticFieldNumeric2D(const std::string&, enum sfile_type);
		MagneticFieldNumeric2D(
			const std::string& name, const std::string& description,
			slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
			slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz,
            slibreal_t raxis, slibreal_t zaxis,
			slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
			slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
		);
		virtual ~MagneticFieldNumeric2D();

        MagneticFieldNumeric2D *Clone();

		slibreal_t *Eval(slibreal_t*);
		slibreal_t *Eval(slibreal_t, slibreal_t, slibreal_t);
		struct magnetic_field_data& EvalDerivatives(slibreal_t*);
		struct magnetic_field_data& EvalDerivatives(slibreal_t, slibreal_t, slibreal_t);
		slibreal_t FindMaxRadius() override;
		slibreal_t FindMinRadius() override;
		std::string& GetDescription();
		std::string& GetName();
        slibreal_t *GetMagneticAxis();
		void Init(
			const std::string& name, const std::string& description,
			slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
			slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz,
            slibreal_t raxis, slibreal_t zaxis,
			slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
			slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
		);
		void InitInterpolation();
		void Load(const std::string&);
		void Load(const std::string&, enum sfile_type);
		/* Transpose a matrix from SFile */
		double **Transpose(double**, sfilesize_t, sfilesize_t);
};

#endif/*_MAGNETIC_FIELD_NUMERIC_2D_H*/
