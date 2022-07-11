/**
 * Implementation of the 'MagneticFieldAnalytical2D'
 * abstract class, from which various analytical
 * magnetic fields can be derived.
 */

#include <cmath>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>

using namespace std;

/**
 * Constructor
 *
 * B0: Magnetic field strength on-axis.
 * Rm: Major radius of tokamak.
 * rminor: Minor radius of tokamak (used to define domain).
 * sigmaB: Toroidal field sign.
 * sigmaI: Plasma current sign.
 * qt: Type of q-factor.
 * qa1: q-factor param1.
 * qa2: q-factor param2 (may be optional).
 *
 * name: Name of magnetic field.
 * description: Description of the magnetic field.
 */
MagneticFieldAnalytical2D::MagneticFieldAnalytical2D(
	slibreal_t B0, slibreal_t Rm, slibreal_t maxis_z, slibreal_t rminor,
    enum MFAFieldSign sigmaB, enum MFAFieldSign sigmaI, enum MFASafetyFactorType qt,
	slibreal_t qa1, slibreal_t qa2
) : MagneticFieldAnalytical2D(
    B0, Rm, maxis_z, rminor, sigmaB, sigmaI, qt, qa1, qa2,
    "Analytical Magnetic Field",
    "A magnetic field with a simple analytical description"
) {}
MagneticFieldAnalytical2D::MagneticFieldAnalytical2D(
	slibreal_t B0, slibreal_t Rm, slibreal_t maxis_z, slibreal_t rminor,
    enum MFAFieldSign sigmaB, enum MFAFieldSign sigmaI, enum MFASafetyFactorType qt,
	slibreal_t qa1, slibreal_t qa2, const string &name, const string &description
) {
	this->B0 	 	  = B0;
	this->Rm 	 	  = Rm;
    this->maxis_z     = maxis_z;
	this->rminor 	  = rminor;
    this->sigmaB      = (sigmaB==MFAFS_CW?(-1.0):(+1.0));
	this->sigmaI      = (sigmaI==MFAFS_CW?(-1.0):(+1.0));
	this->name   	  = name;
	this->description = description;

    this->magnetic_axis[0] = Rm;
    this->magnetic_axis[1] = maxis_z;
    this->jacobian    = new slibreal_t*[3];
    this->jacobian[0] = new slibreal_t[3*3];
    this->jacobian[1] = this->jacobian[0] + 3;
    this->jacobian[2] = this->jacobian[1] + 3;

	this->safety_factor_type   = qt;
	this->safety_factor_param1 = qa1;
	this->safety_factor_param2 = qa2;

    ConstructDomain(Rm, this->maxis_z, rminor);
}

/**
 * Deep-clones this magnetic field object.
 */
MagneticFieldAnalytical2D *MagneticFieldAnalytical2D::Clone() {
    return new MagneticFieldAnalytical2D(
        this->B0, this->Rm, this->maxis_z, this->rminor,
		((this->sigmaB>0)?MFAFS_CCW:MFAFS_CW),
		((this->sigmaI>0)?MFAFS_CCW:MFAFS_CW),
        this->safety_factor_type, this->safety_factor_param1,
        this->safety_factor_param2, this->name, this->description
    );
}

/**
 * Constructs a domain for the 2D analytical magnetic
 * field. The domain is a circle centered at (Rm, zm)
 * with radius rminor.
 *
 * Rm:     Major radius of device (center of domain).
 * zm:     Vertical position of magnetic axis.
 * rminor: Minor radius of device (radius of domain).
 */
void MagneticFieldAnalytical2D::ConstructDomain(const slibreal_t Rm, const slibreal_t zm, const slibreal_t rminor) {
    unsigned int i;
    const unsigned int N = MAGNETIC_FIELD_ANALYTICAL2D_NDOMAINPOINTS;
    slibreal_t *r, *z, theta;

    r = new slibreal_t[N];
    z = new slibreal_t[N];

    for (i = 0; i < N-1; i++) {
        theta = 2.0*M_PI * (((slibreal_t)i) / ((slibreal_t)(N-1)));
        r[i] = Rm + rminor*cos(theta);
        z[i] = zm + rminor*sin(theta);
    }

    r[N-1] = r[0];
    z[N-1] = z[0];

    SetDomain(r, z, N);
    SetSeparatrix(r, z, N);
}

/**
 * Evaluate the magnetic field in the given point.
 *
 * xyz/x,y,z: Cartesian coordinates of the point at
 *   which to evaluate the magnetic field.
 */
slibreal_t *MagneticFieldAnalytical2D::Eval(slibreal_t x, slibreal_t y, slibreal_t z) {
	slibreal_t R = hypot(x, y), r = hypot(Rm-R, z-this->maxis_z),
		sintheta, costheta, sinphi, cosphi, pf,
		iR = 1.0 / R, t1;
	
	if (r != 0) sintheta = (z-this->maxis_z) / r, costheta = (Rm-R) / r;
	else sintheta = 0, costheta = 0;

	sinphi = y * iR, cosphi = x * iR;
	pf = B0*Rm * iR;
	t1 = r / (Rm*__GetSafetyFactor(r));

	retval[0] = pf * (sigmaI * t1 * (sintheta*cosphi) - sigmaB*sinphi);
	retval[1] = pf * (sigmaI * t1 * (sintheta*sinphi) + sigmaB*cosphi);
	retval[2] = pf * (sigmaI * t1 * costheta);

	return retval;
}

/**
 * Evaluate derivatives (grad |B| & curl B) of the
 * magnetic field in the given point.
 *
 * xyz/x,y,z: Cartesian coordinates of the point at
 *   which to evaluate the magnetic field Jacobian.
 *
 * RETURNS a struct containing the magnetic field
 *   vector, its magnitude, the gradient of the
 *   magnitude as well as the curl of the magnetic
 *   field.
 */
struct magnetic_field_data& MagneticFieldAnalytical2D::EvalDerivatives(
    slibreal_t x, slibreal_t y, slibreal_t z
) {
	slibreal_t R = hypot(x, y), Z = z-this->maxis_z,
        r = hypot(Rm-R, Z),
		sintheta, costheta, sinphi, cosphi,
		r_qRm, sqr, q, rdq_dr,
		DBr, DB0, ddr_rB0, iR = 1.0 / R;
	
	if (r != 0) sintheta = Z / r, costheta = (Rm-R) / r;
	else sintheta = 0, costheta = 0;

	sinphi = y * iR, cosphi = x * iR;
	q = __GetSafetyFactor(r);
	rdq_dr = __GetrDqDr(r);
	r_qRm = r / (Rm * q);
	sqr = sqrt(r_qRm*r_qRm + 1);

	DBr = B0*Rm * iR * (r_qRm/(Rm * q)*(1 - rdq_dr)/sqr + costheta * iR * sqr);
	DB0 =-B0*Rm*sintheta * iR*iR * sqr;
	ddr_rB0 = sigmaI * B0 * iR / q * (Rm*iR + 1 - rdq_dr/q);

	magdata.B[0] = B0*Rm*iR * (sigmaI * r_qRm * (sintheta*cosphi) - sigmaB*sinphi);
	magdata.B[1] = B0*Rm*iR * (sigmaI * r_qRm * (sintheta*sinphi) + sigmaB*cosphi);
	magdata.B[2] = B0*Rm*iR * (sigmaI * r_qRm * costheta);
	magdata.Babs = B0*Rm*iR * sqr;
    magdata.J    = this->jacobian;

	magdata.gradB[0] =-DBr*costheta*cosphi + DB0*sintheta*cosphi;
	magdata.gradB[1] =-DBr*costheta*sinphi + DB0*sintheta*sinphi;
	magdata.gradB[2] = DBr*sintheta        + DB0*costheta;

	magdata.curlB[0] =-ddr_rB0 * sinphi;
	magdata.curlB[1] = ddr_rB0 * cosphi;
	magdata.curlB[2] = 0;

    // Bx
    magdata.J[0][0] = B0*Rm*( sigmaI*Z/(q*Rm*R*R)*(1.0-2.0*cosphi*cosphi)  +  2.0*sigmaB/(R*R)*sinphi*cosphi  +  sigmaI*rdq_dr/(q*q*Rm*R)*cosphi*cosphi*sintheta*costheta );
    magdata.J[0][1] = B0*Rm*( 2.0*sigmaB/(R*R)*sinphi*sinphi  -  sigmaB/(R*R)  -  2.0*sigmaI*Z/(q*Rm*R*R)*sinphi*cosphi  +  sigmaI*rdq_dr/(q*q*Rm*R)*sinphi*cosphi*sintheta*costheta );
    magdata.J[0][2] = sigmaI*B0/(q*R)*cosphi * ( 1 - rdq_dr / q * sintheta*sintheta );
    // By
    magdata.J[1][0] = B0*Rm*( sigmaB/(R*R)  -  2.0*sigmaB/(R*R)*cosphi*cosphi  -  2.0*sigmaI*Z/(q*Rm*R*R)*sinphi*cosphi  +  sigmaI*rdq_dr/(q*q*Rm*R)*sinphi*cosphi*sintheta*costheta );
    magdata.J[1][1] = B0*Rm*( sigmaI*Z/(q*Rm*R*R)*(1.0-2.0*sinphi*sinphi)  -  2.0*sigmaB/(R*R)*sinphi*cosphi  +  sigmaI*rdq_dr/(q*q*Rm*R)*sinphi*sinphi*sintheta*costheta );
    magdata.J[1][2] = sigmaI*B0/(q*R)*sinphi*( 1.0 - rdq_dr/q*sintheta*sintheta );
    // Bz
    magdata.J[2][0] = sigmaI*B0/(q*R)*cosphi * ( rdq_dr/q*costheta*costheta  -  Rm/R );
    magdata.J[2][1] = sigmaI*B0/(q*R)*sinphi * ( rdq_dr/q*costheta*costheta  -  Rm/R );
    magdata.J[2][2] = -sigmaI*B0/(q*q*R) * rdq_dr * sintheta*costheta;

	return magdata;
}

/**
 * Evaluate magnetic flux in the given point.
 *
 * x, y, z: Cartesian coordinates of the point to evaluate
 *          the magnetic flux in.
 *
 * RETURNS the value of the magnetic flux in the
 * given point.
 */
slibreal_t MagneticFieldAnalytical2D::EvalFlux(
    slibreal_t x, slibreal_t y, slibreal_t z
) {
    slibreal_t R  = sqrt(x*x + y*y);
    slibreal_t Z  = z-this->maxis_z;
    slibreal_t r2 = Z*Z + (Rm - R)*(Rm - R);
	slibreal_t q  = __GetSafetyFactor(sqrt(r2));

    slibreal_t psi = sigmaI * M_PI*B0*r2 / q;

    return psi;
}

/**
 * Returns the derivatives of the magnetic flux with respect
 * to R and Z respectively, in the given point.
 *
 * x, y, z: Cartesian coordinates of the point to evaluate
 *          the magnetic flux derivatives in.
 *
 * RETURNS the R and Z derivatives of the magnetic flux
 * in the given point as a 2-vector. The vector components are
 *
 *    [0] = d(psi) / dR
 *    [1] = d(psi) / dZ
 */
struct flux_diff *MagneticFieldAnalytical2D::EvalFluxDerivatives(
    slibreal_t x, slibreal_t y, slibreal_t z
) {
    slibreal_t R     = sqrt(x*x + y*y);
    slibreal_t Z     = z-this->maxis_z;
    slibreal_t r2    = Z*Z + (Rm - R)*(Rm - R);
    //slibreal_t d     = sqrt(Rm*Rm - r2);
    //slibreal_t p     = 2.0*M_PI*B0;
	slibreal_t r     = sqrt(r2);
	slibreal_t q     = __GetSafetyFactor(r);
	//slibreal_t rdqdr = __GetrDqDr(r);

	// Toroidal flux
    /*flux_retval.psi     = p*Rm*(1.0 - sqrt(1.0 - r2/(Rm*Rm)));
    flux_retval.dpsi_dR = p * (R - Rm) / d;
    flux_retval.dpsi_dZ = p * z / d;*/

	flux_retval.psi     = sigmaI * M_PI*B0*r2 / q;
	flux_retval.dpsi_dR = sigmaI * M_PI*B0*(R-Rm)/q * 2.0;//- rdqdr / q);
	flux_retval.dpsi_dZ = sigmaI * M_PI*B0*Z/q      * 2.0;//- rdqdr / q);

    return &flux_retval;
}

/**
 * Creates a numeric magnetic field from this analytical
 * magnetic field.
 *
 * nr: Number of radial points in B-fields.
 * nz: Number of vertical points in B-fields.
 */
MagneticFieldNumeric2D *MagneticFieldAnalytical2D::ToNumeric2D(const unsigned int nr, const unsigned int nz) {
    slibreal_t *R, *Z, *Br, *Bphi, *Bz, *Psi,
        *rwall, *zwall, *temp,
        Rmin, Rmax, Zmin, Zmax;
    unsigned int i, j, nwall;

    Rmax = GetRMajor()   + GetRMinor();
    Rmin = GetRMajor()   - GetRMinor();
    Zmax = this->maxis_z + GetRMinor();
    Zmin = this->maxis_z - GetRMinor();

    R    = new slibreal_t[nr];
    Z    = new slibreal_t[nz];
    Br   = new slibreal_t[nr*nz];
    Bphi = new slibreal_t[nr*nz];
    Bz   = new slibreal_t[nr*nz];

    Psi  = new slibreal_t[nr*nz];

    for (i = 0; i < nr; i++)
        R[i] = Rmin + (Rmax-Rmin) * ((slibreal_t)i) / ((slibreal_t)(nr-1.0));
    for (i = 0; i < nz; i++)
        Z[i] = Zmin + (Zmax-Zmin) * ((slibreal_t)i) / ((slibreal_t)(nz-1.0));

    for (i = 0; i < nr; i++) {
        for (j = 0; j < nz; j++) {
            Vector<3> B = Eval(R[i], 0.0, Z[j]);

            Br  [i + j*nr] = B[0];
            Bphi[i + j*nr] = B[1];
            Bz  [i + j*nr] = B[2];
        }
    }

    nwall = GetNDomain();
    rwall = new slibreal_t[nwall];
    zwall = new slibreal_t[nwall];

    temp = GetRDomain();
    for (i = 0; i < nwall; i++)
        rwall[i] = temp[i];
    temp = GetZDomain();
    for (i = 0; i < nwall; i++)
        zwall[i] = temp[i];

    return new MagneticFieldNumeric2D(
        GetName(), GetDescription(), this->COCOS,
        R, Z, nr, nz, Br, Bphi, Bz, Psi,
        GetMagneticAxisR(), GetMagneticAxisZ(),
        nullptr, nullptr, 0,
        rwall, zwall, nwall
    );
}

