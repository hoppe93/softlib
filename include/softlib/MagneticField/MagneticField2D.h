#ifndef _MAGNETICFIELD2D_H
#define _MAGNETICFIELD2D_H

#include <softlib/config.h>
#include <softlib/SFile.h>
#include <softlib/Vector.h>

struct magnetic_field_data {
	slibreal_t B[3];
	slibreal_t gradB[3];
	slibreal_t curlB[3];
	slibreal_t Babs;
    slibreal_t **J;     // Jacobian (3x3 matrix)
};
struct flux_diff {
    slibreal_t psi;
    slibreal_t dpsi_dR, dpsi_dZ;
};

class MagneticField2D {
	private:
		slibreal_t *rdomain, *zdomain;
        unsigned int ndomain;
		slibreal_t maxradius=-1, minradius=-1;		// Maximum and minimum radius from which particles can be dropped

        slibreal_t *rseparatrix=nullptr, *zseparatrix=nullptr;
        unsigned int nseparatrix;

        bool bounds_set = false;
        struct {
            slibreal_t rmin, rmax, zmin, zmax;
        } bounds = {0,0,0,0};

		bool __IntersectsDomain3D(
			slibreal_t, slibreal_t, slibreal_t, slibreal_t,
			slibreal_t, slibreal_t, slibreal_t,
			slibreal_t, slibreal_t, slibreal_t
		);
	protected:
		slibreal_t magnetic_axis[2];
	public:
		std::string name, description;
		struct magnetic_field_data magdata;

		MagneticField2D();
        virtual ~MagneticField2D();

        virtual MagneticField2D *Clone() = 0;

		slibreal_t *Eval(Vector<3>&);
		slibreal_t *Eval(slibreal_t*);
		virtual slibreal_t *Eval(slibreal_t, slibreal_t, slibreal_t) = 0;
		struct magnetic_field_data& EvalDerivatives(Vector<3>&);
		struct magnetic_field_data& EvalDerivatives(slibreal_t*);
		virtual struct magnetic_field_data& EvalDerivatives(slibreal_t, slibreal_t, slibreal_t) = 0;
        slibreal_t EvalFlux(Vector<3>&);
        slibreal_t EvalFlux(slibreal_t*);
        virtual slibreal_t EvalFlux(slibreal_t, slibreal_t, slibreal_t) = 0;
        struct flux_diff *EvalFluxDerivatives(Vector<3>&);
        struct flux_diff *EvalFluxDerivatives(slibreal_t*);
        virtual struct flux_diff *EvalFluxDerivatives(slibreal_t, slibreal_t, slibreal_t) = 0;

		std::string& GetDescription() { return this->description; }
        slibreal_t *GetMagneticAxis();
		slibreal_t GetMagneticAxisR();
		slibreal_t GetMagneticAxisZ();
		slibreal_t GetMaxRadius();
		slibreal_t GetMinRadius();
		std::string& GetName() { return this->name; }
        virtual bool HasMagneticFlux() = 0;

        slibreal_t *GetRDomain() { return this->rdomain; }
        slibreal_t *GetZDomain() { return this->zdomain; }
        unsigned int GetNDomain() { return this->ndomain; }

        slibreal_t *GetRSeparatrix() { return this->rseparatrix; }
        slibreal_t *GetZSeparatrix() { return this->zseparatrix; }
        unsigned int GetNSeparatrix() { return this->nseparatrix; }

        void GetDomainBounds(slibreal_t&, slibreal_t&, slibreal_t&, slibreal_t&);

		bool CrossesDomain(slibreal_t*);
		bool CrossesDomain(slibreal_t*, slibreal_t*);
		bool CrossesDomain(slibreal_t, slibreal_t, slibreal_t);
		bool CrossesDomain(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
		virtual slibreal_t FindMaxRadius();
		virtual slibreal_t FindMinRadius();
        virtual void FindMaxZ(const slibreal_t, slibreal_t&, slibreal_t&);
		slibreal_t FindMaxRadius(unsigned int, slibreal_t*, slibreal_t*);
		slibreal_t FindMinRadius(unsigned int, slibreal_t*, slibreal_t*);
        void FindMaxZ(unsigned int, const slibreal_t*, const slibreal_t*, const slibreal_t, slibreal_t&, slibreal_t&);
        void FindMaxMin(unsigned int, const slibreal_t*, const slibreal_t*, const slibreal_t, slibreal_t&, slibreal_t&);
		bool IntersectsDomain3D(slibreal_t*, slibreal_t*, bool has_outer_wall=true);
		bool IntersectsDomain3D(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, bool has_outer_wall=true);
        void SetDomain(slibreal_t*, slibreal_t*, unsigned int);
        void SetSeparatrix(slibreal_t*, slibreal_t*, unsigned int);
};

#endif/*_MAGNETICFIELD2D_H*/
