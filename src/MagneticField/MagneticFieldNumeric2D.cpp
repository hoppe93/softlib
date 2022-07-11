/**
 * Implementation of a 2D numeric magnetic field.
 *
 */

#include <cmath>
#include <string>
using namespace std;

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>
#include <softlib/SFile.h>
#include <softlib/SFileException.h>
#include <softlib/SOFTLibException.h>

/**
 * Minimal contsructor. Intended for use with
 * derived classes that load different kinds of files.
 */
MagneticFieldNumeric2D::MagneticFieldNumeric2D() : MagneticField2D() {
	BaseInit();
}

/**
 * Construct from file
 * 
 * filename: Name of file containing magnetic field data.
 * ftype: Reader to use for reading file.
 */
MagneticFieldNumeric2D::MagneticFieldNumeric2D(const string& filename) : MagneticField2D() {
	Load(filename);
	InitInterpolation();

	BaseInit();

	// Verify that COCOS number makes sense
	this->VerifyCocos(this->cocos);
}
MagneticFieldNumeric2D::MagneticFieldNumeric2D(const string& filename, enum sfile_type ftype) : MagneticField2D() {
	Load(filename, ftype);
	InitInterpolation();

	BaseInit();
}
MagneticFieldNumeric2D::~MagneticFieldNumeric2D() {
	if (interpval != nullptr)
		delete [] interpval;
	if (jacobian != nullptr) {
		delete [] jacobian[0];
		delete [] jacobian;
	}

	if (R != nullptr) delete [] R;
	if (Z != nullptr) delete [] Z;

	if (Br != nullptr) delete [] Br;
	if (Bphi != nullptr) delete [] Bphi;
	if (Bz != nullptr) delete [] Bz;

	if (nwall > 0)
		delete [] rwall, delete [] zwall;
	if (nsep > 0)
		delete [] rsep, delete [] zsep;
}

/**
 * Do basic initialization of this class.
 */
void MagneticFieldNumeric2D::BaseInit() {
	interpval = new slibreal_t[3];

    // The 'jacobian' is a 4-by-3 matrix here.
    // The first 3 rows (3-by-3 matrix) is the
    // actual jacobian, while the last row
    // (1-by-3 vector) is the magnetic field vector.
	jacobian  = new slibreal_t*[4];
	jacobian[0] = new slibreal_t[4*3];
	jacobian[1] = jacobian[0] + 3;
	jacobian[2] = jacobian[1] + 3;
	jacobian[3] = jacobian[2] + 3;
}

/**
 * Construct from data *
 * name: Name of magnetic field.
 * description: Description of magnetic field.
 * cocos: COCOS number for magnetic field.
 * R, Z: R and Z coordinates of magnetic field components. Simple vectors.
 * Br: Radial component of magnetic field.
 * Bphi: Toroidal component of magnetic field.
 * Bz: Vertical component of magnetic field.
 * rsep, zsep: R and Z points of separatrix data. Can be nullptr if rwall/zwall is not.
 * rwall, zwall: R and Z points of wall data. Can be nullptr if rsep/zsep is not.
 */
MagneticFieldNumeric2D::MagneticFieldNumeric2D(
	const string& name, const string& description, const int cocos,
	slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
	slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz,
    slibreal_t *Psi, slibreal_t raxis, slibreal_t zaxis,
	slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
	slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
) {
	interpval = new slibreal_t[3];
	jacobian  = new slibreal_t*[4];
	jacobian[0] = new slibreal_t[4*3];
	jacobian[1] = jacobian[0] + 3;
	jacobian[2] = jacobian[1] + 3;
	jacobian[3] = jacobian[2] + 3;

	Init(name, description, cocos, R, Z, nr, nz, Br, Bphi, Bz, Psi, raxis, zaxis, rsep, zsep, nsep, rwall, zwall, nwall);
}

/**
 * Deep-clones this magnetic field object.
 * Copies all arrays and passes them to a new
 * MagneticFieldNumeric2D object.
 */
MagneticFieldNumeric2D *MagneticFieldNumeric2D::Clone() {
    unsigned int i, nB;
    slibreal_t *nR, *nZ, *nBr, *nBphi, *nBz, *nPsi=nullptr,
        *nrsep, *nzsep, *nrwall, *nzwall;

    nB     = this->nr*this->nz;
    nR     = new slibreal_t[this->nr];
    nZ     = new slibreal_t[this->nz];
    nBr    = new slibreal_t[nB];
    nBphi  = new slibreal_t[nB];
    nBz    = new slibreal_t[nB];
    nrsep  = new slibreal_t[this->nsep];
    nzsep  = new slibreal_t[this->nsep];
    nrwall = new slibreal_t[this->nwall];
    nzwall = new slibreal_t[this->nwall];

    for (i = 0; i < this->nr; i++)
        nR[i] = this->R[i];
    for (i = 0; i < this->nz; i++)
        nZ[i] = this->Z[i];
    for (i = 0; i < nB; i++)
        nBr[i] = this->Br[i];
    for (i = 0; i < nB; i++)
        nBphi[i] = this->Bphi[i];
    for (i = 0; i < nB; i++)
        nBz[i] = this->Bz[i];
    for (i = 0; i < this->nsep; i++)
        nrsep[i] = this->rsep[i];
    for (i = 0; i < this->nsep; i++)
        nzsep[i] = this->zsep[i];
    for (i = 0; i < this->nwall; i++)
        nrwall[i] = this->rwall[i];
    for (i = 0; i < this->nwall; i++)
        nzwall[i] = this->zwall[i];

    if (this->Psi != nullptr) {
        nPsi   = new slibreal_t[nB];

        for (i = 0; i < nB; i++)
            nPsi[i] = this->Psi[i];
    }

    return new MagneticFieldNumeric2D(
        this->name, this->description, this->cocos,
        nR, nZ, this->nr, this->nz,
        nBr, nBphi, nBz, nPsi,
        this->magnetic_axis[0], this->magnetic_axis[1],
        nrsep, nzsep, this->nsep,
        nrwall, nzwall, this->nwall
    );
}

/**
 * Locates the maximum allowed radius at which
 * particles can be dropped in the plasma. This
 * is done by finding the maximum radius of the
 * domain at the vertical level of the magnetic
 * axis. If available, the separatrix is chosen
 * as the domain, even if it's not used as _the_
 * domain. Otherwise, the wall is used.
 */
slibreal_t MagneticFieldNumeric2D::FindMaxRadius() {
	if (this->nsep > 0)
		return MagneticField2D::FindMaxRadius(nsep, rsep, zsep);
	else
		return MagneticField2D::FindMaxRadius(nwall, rwall, zwall);
}

/**
 * Locates the minimum allowed radius at which
 * particles can be dropped in the plasma. This
 * is done by finding the minimum radius of the
 * domain at the vertical level of the magnetic
 * axis. If available, the separatrix is chosen
 * as the domain, even if it's not used as _the_
 * domain. Otherwise, the wall is used.
 */
slibreal_t MagneticFieldNumeric2D::FindMinRadius() {
	if (this->nsep > 0)
		return MagneticField2D::FindMinRadius(nsep, rsep, zsep);
	else
		return MagneticField2D::FindMinRadius(nwall, rwall, zwall);
}

/**
 * Initialize this MagneticFieldNumeric2D object,
 * setting all properties in one call.
 *
 * name: Name of magnetic field.
 * description: Description of magnetic field.
 * cocos: COCOS number for magnetic field.
 * R, Z: R and Z coordinates of the magnetic field components. Simple vectors.
 * Br: Radial component of the magnetic field.
 * Bphi: Toroidal component of the magnetic field.
 * Bz: Vertical component of the magnetic field.
 * raxis, zaxis: R and Z coordinates of the magnetic axis.
 * rsep, zsep: R and Z points of separatrix data. Can be nullptr if rwall/zwall is not.
 * rwall, zwall: R and Z points of wall data. Can be nullptr if rsep/zsep is not.
 */
void MagneticFieldNumeric2D::Init(
	const string& name, const string& description, const int cocos,
	slibreal_t *R, slibreal_t *Z, unsigned int nr, unsigned int nz,
	slibreal_t *Br, slibreal_t *Bphi, slibreal_t *Bz,
    slibreal_t *Psi, slibreal_t raxis, slibreal_t zaxis,
	slibreal_t *rsep, slibreal_t *zsep, unsigned int nsep,
	slibreal_t *rwall, slibreal_t *zwall, unsigned int nwall
) {
	if ((rsep == nullptr || zsep == nullptr) &&
		(rwall == nullptr || zwall == nullptr))
		throw SOFTLibException("Wall and separatrix cannot both be nullptr.");

	this->name = name;
	this->description = description;
	this->cocos = cocos;
	this->R = R;
	this->Z = Z;
	this->nr = nr;
	this->nz = nz;
	this->Br = Br;
	this->Bphi = Bphi;
	this->Bz = Bz;
    this->Psi = Psi;
    this->magnetic_axis[0] = raxis;
    this->magnetic_axis[1] = zaxis;
	this->rsep = rsep;
	this->zsep = zsep;
	this->nsep = nsep;
	this->rwall = rwall;
	this->zwall = zwall;
	this->nwall = nwall;

    this->hasFluxCoordinates = (Psi != nullptr);

    if (rwall != nullptr && zwall != nullptr)
        SetDomain(this->rwall, this->zwall, this->nwall);
    else
        SetDomain(this->rsep, this->zsep, this->nsep);

    // Always set 'separatrix', even if it is 'null'.
    // (This quantity is only used as metadata)
    SetSeparatrix(this->rsep, this->zsep, this->nsep);

	InitInterpolation();
}

/**
 * Initialize the interpolation library, used for
 * interpolating in the numeric magnetic field.
 */
void MagneticFieldNumeric2D::InitInterpolation() {
#ifdef INTERP_SPLINTER
	/* Use SPLINTER for interpolation */
#	error "SPLINTER support has not been implemented yet."
#else
	/* Use GSL for interpolation */
	ra = gsl_interp_accel_alloc();
	za = gsl_interp_accel_alloc();

	sBr   = gsl_spline2d_alloc(gsl_interp2d_bicubic, nr, nz);
	sBphi = gsl_spline2d_alloc(gsl_interp2d_bicubic, nr, nz);
	sBz   = gsl_spline2d_alloc(gsl_interp2d_bicubic, nr, nz);

	gsl_spline2d_init(sBr, R, Z, Br, nr, nz);
	gsl_spline2d_init(sBphi, R, Z, Bphi, nr, nz);
	gsl_spline2d_init(sBz, R, Z, Bz, nr, nz);

    if (HasMagneticFlux()) {
        sPsi = gsl_spline2d_alloc(gsl_interp2d_bicubic, nr, nz);
        gsl_spline2d_init(sPsi, R, Z, Psi, nr, nz);
    }

	rmin = R[0];
	rmax = R[nr-1];
	zmin = Z[0];
	zmax = Z[nz-1];
#endif
}

/**
 * Evaluate the magnetic field vector in the given point.
 *
 * xyz/x,y,z: Cartesian coordinates of the point at which
 *            to evaluate the magnetic field.
 */
slibreal_t *MagneticFieldNumeric2D::Eval(slibreal_t *xyz) {
	return Eval(xyz[0], xyz[1], xyz[2]);
}
slibreal_t *MagneticFieldNumeric2D::Eval(slibreal_t x, slibreal_t y, slibreal_t z) {
	slibreal_t r = hypot(x, y), br, bphi, bz, sin0, cos0;

	/* Make sure requested point is within domain */
		 if (r < rmin) r = rmin;
	else if (r > rmax) r = rmax;

		 if (z < zmin) z = zmin;
	else if (z > zmax) z = zmax;

    br   = (slibreal_t)gsl_spline2d_eval(sBr,   r, z, ra, za);
    bphi = (slibreal_t)gsl_spline2d_eval(sBphi, r, z, ra, za);
    bz   = (slibreal_t)gsl_spline2d_eval(sBz,   r, z, ra, za);

    if (r != 0) { sin0 = y/r, cos0 = x/r; }
    else { sin0 = 0, cos0 = 1; }
    interpval[0] = br * cos0 - bphi * sin0;
    interpval[1] = br * sin0 + bphi * cos0;
    interpval[2] = bz;

	return interpval;
}

/**
 * Calculate the magnetic field jacobian in the
 * given point. The returned object is a 4-by-3
 * matrix. The first 3 rows contain the Jacobian
 * (in order (r, theta, phi)), while the 4th row
 * contains [0] = Bx, [1] = By, [2] = Bz.
 *
 * xyz/x,y,z: Cartesian coordinates of the point at which
 *            to evaluate the magnetic field.
 */
slibreal_t **MagneticFieldNumeric2D::__EvalJacobian(slibreal_t x, slibreal_t y, slibreal_t z) {
    slibreal_t r = hypot(x, y),
        sin0, cos0, sin20, cos20, sc0,
        dsin0_dx, dcos0_dx, dsin0_dy, dcos0_dy,
        dBr_dr, dBr_dz, dB0_dr, dB0_dz, dBz_dr, dBz_dz,
        br, b0, bz;

         if (r < rmin) r = rmin;
    else if (r > rmax) r = rmax;
         if (z < zmin) z = zmin;
    else if (z > zmax) z = zmax;

    /* dBr/dr */
    dBr_dr = gsl_spline2d_eval_deriv_x(sBr,   r, z, ra, za);
    /* dBr/dz */
    dBr_dz = gsl_spline2d_eval_deriv_y(sBr,   r, z, ra, za);
    /* dB0/dr */
    dB0_dr = gsl_spline2d_eval_deriv_x(sBphi, r, z, ra, za);
    /* dB0/dz */
    dB0_dz = gsl_spline2d_eval_deriv_y(sBphi, r, z, ra, za);
    /* dBz/dr */
    dBz_dr = gsl_spline2d_eval_deriv_x(sBz,   r, z, ra, za);
    /* dBz/dz */
    dBz_dz = gsl_spline2d_eval_deriv_y(sBz,   r, z, ra, za);

    br = gsl_spline2d_eval(sBr,   r, z, ra, za);
    b0 = gsl_spline2d_eval(sBphi, r, z, ra, za);
    bz = gsl_spline2d_eval(sBz,   r, z, ra, za);

    if (r != 0) {
        sin0 = y/r, cos0 = x/r;
        sin20 = sin0*sin0, cos20 = cos0*cos0;
        sc0 = sin0*cos0;

        dsin0_dx = -sc0/r,  dcos0_dx = sin20/r;
        dsin0_dy = cos20/r, dcos0_dy = -sc0/r;
    } else {
        sin0 = sin20 = sc0 = 0;
        cos0 = cos20 = 1;

        dsin0_dx = 0, dcos0_dx = 0;
        dsin0_dy = 0, dcos0_dy = 0;
    }

    interpval[0] = br * cos0 - b0 * sin0;
    interpval[1] = br * sin0 + b0 * cos0;
    interpval[2] = bz;

    /* Bx */
    jacobian[0][0] = cos20*dBr_dr + br*dcos0_dx - sc0  *dB0_dr - b0*dsin0_dx;
    jacobian[0][1] = sc0  *dBr_dr + br*dcos0_dy - sin20*dB0_dr - b0*dsin0_dy;
    jacobian[0][2] = cos0 *dBr_dz - sin0*dB0_dz;
    /* By */
    jacobian[1][0] = sc0  *dBr_dr + br*dsin0_dx + cos20*dB0_dr + b0*dcos0_dx;
    jacobian[1][1] = sin20*dBr_dr + br*dsin0_dy + sc0  *dB0_dr + b0*dcos0_dy;
    jacobian[1][2] = sin0 *dBr_dz + cos0*dB0_dz;
    /* Bz */
    jacobian[2][0] = cos0*dBz_dr;
    jacobian[2][1] = sin0*dBz_dr;
    jacobian[2][2] = dBz_dz;

    return jacobian;
}

struct magnetic_field_data& MagneticFieldNumeric2D::EvalDerivatives(slibreal_t x, slibreal_t y, slibreal_t z) {
	slibreal_t **J = __EvalJacobian(x, y, z);

	slibreal_t
		dBx_dx=J[0][0],  dBx_dy=J[0][1],  dBx_dz=J[0][2],
		dBy_dx=J[1][0],  dBy_dy=J[1][1],  dBy_dz=J[1][2],
		dBz_dx=J[2][0],  dBz_dy=J[2][1],  dBz_dz=J[2][2],
		Bx=interpval[0], By=interpval[1], Bz=interpval[2],
		aBx=fabs(Bx),    aBy=fabs(By),    aBz=fabs(Bz),
		r1, r2, m, Babs;

	if (aBx > aBy && aBx > aBz)
		r1 = aBy / aBx, r2 = aBz / aBx, m = aBx;
	else if (aBy > aBz)
		r1 = aBx / aBy, r2 = aBz / aBy, m = aBy;
	else
		r1 = aBx / aBz, r2 = aBy / aBz, m = aBz;
	
	Babs = m * sqrt(1.0 + r1*r1 + r2*r2);
	magdata.Babs = Babs;
	magdata.B[0] = Bx;
	magdata.B[1] = By;
	magdata.B[2] = Bz;
    magdata.J    = J;

	magdata.gradB[0] = (Bx*dBx_dx + By*dBy_dx + Bz*dBz_dx) / Babs;
	magdata.gradB[1] = (Bx*dBx_dy + By*dBy_dy + Bz*dBz_dy) / Babs;
	magdata.gradB[2] = (Bx*dBx_dz + By*dBy_dz + Bz*dBz_dz) / Babs;

    magdata.curlB[0] = dBz_dy - dBy_dz;
    magdata.curlB[1] = dBx_dz - dBz_dx;
    magdata.curlB[2] = dBy_dx - dBx_dy;

	return magdata;
}

/**
 * Evaluate the magnetic flux in the given point.
 * 
 * x, y, z: Cartesian coordinates of the point to
 *          evaluate the magnetic flux in.
 */
slibreal_t MagneticFieldNumeric2D::EvalFlux(
    slibreal_t x, slibreal_t y, slibreal_t z
) {
    if (!HasMagneticFlux())
        throw SOFTLibException("Trying to evaluate magnetic flux coordinates which have not been defined.");

    slibreal_t r = hypot(x, y);
    return (slibreal_t)gsl_spline2d_eval(sPsi, r, z, ra, za);
}

/**
 * Evaluate the R and Z derivatives of the magnetic
 * flux in the given point.
 *
 * x, y, z: Cartesian coordinates of the point to
 *          evaluate the R and Z derivatives of the
 *          magnetic flux in.
 *
 * RETURNS a 2-vector containing the derivatives of
 * the magnetic flux. The contents are as follows:
 *
 *    [0] = d(psi) / dR
 *    [1] = d(psi) / dZ
 */
struct flux_diff *MagneticFieldNumeric2D::EvalFluxDerivatives(
    slibreal_t x, slibreal_t y, slibreal_t z
) {
    if (!HasMagneticFlux())
        throw SOFTLibException("Trying to evaluate magnetic flux coordinates which have not been defined.");

    slibreal_t r = hypot(x, y);

    flux_data.psi     = (slibreal_t)gsl_spline2d_eval(sPsi, r, z, ra, za);
    flux_data.dpsi_dR = (slibreal_t)gsl_spline2d_eval_deriv_x(sPsi, r, z, ra, za);
    flux_data.dpsi_dZ = (slibreal_t)gsl_spline2d_eval_deriv_y(sPsi, r, z, ra, za);

    return &flux_data;
}

/**
 * Loads a 2D numeric magnetic field from the
 * file named 'filename'.
 *
 * filename: Name of file to load magnetic field from.
 * (optional) ftype: SFile-interface to use when loading file.
 */
void MagneticFieldNumeric2D::Load(const string& filename) {
	Load(filename, SFile::TypeOfFile(filename));
}
void MagneticFieldNumeric2D::Load(const string& filename, enum sfile_type ftype) {
	string name, desc;
	double *_maxis, *_Rgrid, *_Zgrid, **_Br, **_Bphi, **_Bz, **_Psi=nullptr,
		**_temp, *_rwall=nullptr, *_zwall=nullptr, *_rsep=nullptr, *_zsep=nullptr,
        *_verBr=nullptr, *_verBphi=nullptr, *_verBz=nullptr;
	sfilesize_t fs[2], i, j, nverBr=0, nverBphi=0, nverBz=0;

	SFile *sf = SFile::Create(filename, SFILE_MODE_READ, ftype);

	/* Load all variables (to native types) */
	name = sf->GetString("name");
	desc = sf->GetString("desc");

	_maxis = sf->GetList("maxis", fs);
	if (fs[0] != 2)
		throw SFileException(filename+": The magnetic axis must be given as a point with two coordinates.");
	
	_Rgrid = sf->GetList("r", fs);
	this->nr = fs[0];
	_Zgrid = sf->GetList("z", fs);
	this->nz = fs[0];

	_Br = sf->GetDoubles("Br", fs);
	if (fs[0] != (sfilesize_t)this->nz)
		_Br = Transpose(_Br, fs[0], fs[1]);
	_Bphi = sf->GetDoubles("Bphi", fs);
	if (fs[0] != (sfilesize_t)this->nz)
		_Bphi = Transpose(_Bphi, fs[0] ,fs[1]);
	_Bz = sf->GetDoubles("Bz", fs);
	if (fs[0] != (sfilesize_t)this->nz)
		_Bz = Transpose(_Bz, fs[0], fs[1]);

    /* Load magnetic flux */
    try {
        _Psi = sf->GetDoubles("Psi", fs);
        if (fs[0] != (sfilesize_t)this->nz)
            _Psi = Transpose(_Psi, fs[0], fs[1]);
    } catch (SFileException& ex) {}

	/* Load COCOS number */
	try {
		this->cocos = sf->GetInt("cocos");
	} catch (SFileException& ex) {
		// Assume COCOS number = 1
		this->cocos = 1;
	}

    /* Load verification vectors (if present) */
    try {
        _verBr = sf->GetList("verBr", fs);
        nverBr = fs[0];
	} catch (SFileException& ex) {}
    try {
        _verBphi = sf->GetList("verBphi", fs);
        nverBphi = fs[0];
    } catch (SFileException& ex) {}
    try {
        _verBz = sf->GetList("verBz", fs);
        nverBz = fs[0];
    } catch (SFileException& ex) {}

    /* Verify magnetic fields (check along R dimension) */
    if (nverBr > this->nr)
        throw SFileException(filename+": Verification vector 'verBr' has more than nr (= %u) elements.", this->nr);
    for (i = 0; i < nverBr; i++) {
        if (_Br[0][i] != _verBr[i])
            throw SFileException(filename+": Magnetic field component 'Br' did not match its verification vector.");
    }

    if (nverBphi > this->nr)
        throw SFileException(filename+": Verification vector 'verBphi' has more than nr (= %u) elements.", this->nr);
    for (i = 0; i < nverBphi; i++) {
        if (_Bphi[0][i] != _verBphi[i])
            throw SFileException(filename+": Magnetic field component 'Bphi' did not match its verification vector.");
    }

    if (nverBz > this->nr)
        throw SFileException(filename+": Verification vector 'verBz' has more than nr (= %u) elements.", this->nr);
    for (i = 0; i < nverBz; i++) {
        if (_Bz[0][i] != _verBz[i])
            throw SFileException(filename+": Magnetic field component 'Bz' did not match its verification vector.");
    }

    delete [] _verBr;
    delete [] _verBphi;
    delete [] _verBz;

	try {
		_temp = sf->GetDoubles("wall", fs);
        if (fs[0] == 2) {
            _rwall = _temp[0];
            _zwall = _temp[1];
            this->nwall = fs[1];
            delete [] _temp;
        } else {
            this->nwall = fs[0];
            _rwall = new slibreal_t[this->nwall];
            _zwall = new slibreal_t[this->nwall];

            for (i = 0; i < this->nwall; i++) {
                _rwall[i] = _temp[i][0];
                _zwall[i] = _temp[i][1];
            }

            delete [] _temp[0];
            delete [] _temp;
        }
	} catch (SFileException& ex) {
		this->nwall = 0;
		this->rwall = nullptr;
		this->zwall = nullptr;
	}

	try {
		_temp = sf->GetDoubles("separatrix", fs);
        if (fs[0] == 2) {
            _rsep = _temp[0];
            _zsep = _temp[1];
            this->nsep = fs[1];
            delete [] _temp;
        } else {
            this->nsep = fs[0];
            _rsep = new slibreal_t[this->nsep];
            _zsep = new slibreal_t[this->nsep];

            for (i = 0; i < this->nsep; i++) {
                _rsep[i] = _temp[i][0];
                _zsep[i] = _temp[i][1];
            }

            delete [] _temp[0];
            delete [] _temp;
        }
	} catch (SFileException& ex) {
		this->nsep = 0;
		this->rsep = nullptr;
		this->zsep = nullptr;
	}

	sf->Close();
	
	/* Convert doubles to 'slibreal_t' */
	this->name = name;
	this->description = desc;
	this->magnetic_axis[0] = (slibreal_t)_maxis[0];
	this->magnetic_axis[1] = (slibreal_t)_maxis[1];
	this->R = new slibreal_t[this->nr];
	this->Z = new slibreal_t[this->nz];

	for (i = 0; i < this->nr; i++)
		this->R[i] = (slibreal_t)_Rgrid[i];
	for (i = 0; i < this->nz; i++)
		this->Z[i] = (slibreal_t)_Zgrid[i];
	
	this->Br = new slibreal_t[this->nr*this->nz];
	this->Bphi = new slibreal_t[this->nr*this->nz];
	this->Bz = new slibreal_t[this->nr*this->nz];

	for (i = 0; i < this->nr; i++)
		for (j = 0; j < this->nz; j++)
			this->Br[i + j*this->nr] = (slibreal_t)_Br[j][i];

	for (i = 0; i < this->nr; i++)
		for (j = 0; j < this->nz; j++)
			this->Bphi[i + j*this->nr] = (slibreal_t)_Bphi[j][i];

	for (i = 0; i < this->nr; i++)
		for (j = 0; j < this->nz; j++)
		    this->Bz[i + j*this->nr] = (slibreal_t)_Bz[j][i];

    if (_Psi != nullptr) {
        this->Psi = new slibreal_t[this->nr*this->nz];

        for (i = 0; i < this->nr; i++)
            for (j = 0; j < this->nz; j++)
                this->Psi[i + j*this->nr] = (slibreal_t)_Psi[j][i];
        
        this->hasFluxCoordinates = true;
    } else
        this->Psi = nullptr;
	
	delete [] _maxis;
	delete [] _Rgrid, delete [] _Zgrid;
	delete [] _Br[0], delete [] _Bphi[0], delete [] _Bz[0];
	delete [] _Br, delete [] _Bphi, delete [] _Bz;

	if (_rwall == nullptr && _zsep == nullptr)
		throw SOFTLibException(filename+": At least one of 'wall' and 'separatrix' must be provided in the magnetic equilibrium file.");

	if (_rwall != nullptr) {
		this->rwall = new slibreal_t[this->nwall];
		this->zwall = new slibreal_t[this->nwall];

		for (i = 0; i < this->nwall; i++) {
			this->rwall[i] = (slibreal_t)_rwall[i];
			this->zwall[i] = (slibreal_t)_zwall[i];
		}
		
		delete [] _rwall;
	}
	if (_rsep != nullptr) {
		this->rsep  = new slibreal_t[this->nsep];
		this->zsep  = new slibreal_t[this->nsep];

		for (i = 0; i < this->nsep; i++) {
			this->rsep[i] = (slibreal_t)_rsep[i];
			this->zsep[i] = (slibreal_t)_zsep[i];
		}

		delete [] _rsep;
	}

    if (rwall != nullptr && zwall != nullptr)
        SetDomain(this->rwall, this->zwall, this->nwall);
    else
        SetDomain(this->rsep, this->zsep, this->nsep);

    // Always set 'separatrix', even if it is 'null'.
    // (This quantity is only used as metadata)
    SetSeparatrix(this->rsep, this->zsep, this->nsep);
}

/**
 * Return the COCOS parameter \sigma_{B_p}, indicating how the
 * poloidal flux is defined relative to the poloidal field.
 */
slibreal_t MagneticFieldNumeric2D::GetCocosSigmaBp() {
	switch (this->cocos) {
		case 1:  case 2:  case 5:  case 6:
		case 11: case 12: case 15: case 16:
			return 1.0;
		
		case 3:  case 4:  case 7:  case 8:
		case 13: case 14: case 17: case 18:
			return -1.0;

		// Invalid COCOS number
		default: return 1.0;
	}
}

/**
 * Return the COCOS parameter e_{B_p}, indicating whether
 * the poloidal flux is normalized by 2*pi (e_{B_p}=1) or
 * not (e_{B_p}=0).
 */
slibreal_t MagneticFieldNumeric2D::GetCocosEBp() {
	if (this->cocos >= 10)
		return 1;
	else
		return 0;
}

/**
 * Transpose the given matrix of size rows-by-cols.
 * The returned matrix will be of size cols-by-rows.
 *
 * a: Matrix to transpose.
 * rows: Number of rows in the input matrix.
 * cols: Number of cols in the input matrix.
 */
double **MagneticFieldNumeric2D::Transpose(double **a, sfilesize_t rows, sfilesize_t cols) {
	sfilesize_t i, j;
	double **r = new double*[cols],
			*p = new double[rows*cols];
	
	for (i = 0; i < cols; i++) {
		r[i] = p + i*rows;
		for (j = 0; j < rows; j++) {
			r[i][j] = a[j][i];
		}
	}

	delete [] a[0];
	delete [] a;

	return r;
}

/**
 * Save this numeric magnetic field to the output
 * file with the given name.
 *
 * filename: Name of output file to save magnetic field to.
 */
void MagneticFieldNumeric2D::Save(const string& filename) {
	Save(filename, SFile::TypeOfFile(filename));
}
void MagneticFieldNumeric2D::Save(const string& filename, enum sfile_type ftype) {
	SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE, ftype);
	sfilesize_t dims[2] = {this->nz, this->nr};

	sf->WriteMultiArray("Bphi", this->Bphi, 2, dims);
	sf->WriteMultiArray("Br",   this->Br,   2, dims);
	sf->WriteMultiArray("Bz",   this->Bz,   2, dims);

	int32_t cc = this->cocos;
	sf->WriteInt32List("cocos", &cc, 1);

	sf->WriteString("desc", this->description);
	sf->WriteString("name", this->name);

	sf->WriteList("maxis", this->magnetic_axis, 2);

	sf->WriteMultiArray("Psi", this->Psi, 2, dims);

	sf->WriteList("r", this->R, this->nr);
	sf->WriteList("z", this->Z, this->nz);

	sf->WriteList("verBphi", this->Bphi, this->nr);
	sf->WriteList("verBr",   this->Br,   this->nr);
	sf->WriteList("verBz",   this->Bz,   this->nr);

	if (this->nsep > 0) {
		double **separatrix = new double*[2];
		separatrix[0] = new double[2*this->nsep];
		separatrix[1] = separatrix[0] + this->nsep;

		for (unsigned int i = 0; i < this->nsep; i++) {
			separatrix[0][i] = rsep[i];
			separatrix[1][i] = zsep[i];
		}

		sf->WriteArray("separatrix", separatrix, 2, this->nsep);

		delete [] separatrix[0];
		delete [] separatrix;
	}

	if (this->nwall > 0) {
		double **wall = new double*[2];
		wall[0] = new double[2*this->nwall];
		wall[1] = wall[0] + this->nwall;

		for (unsigned int i = 0; i < this->nwall; i++) {
			wall[0][i] = rwall[i];
			wall[1][i] = zwall[i];
		}

		sf->WriteArray("wall", wall, 2, this->nwall);

		delete [] wall[0];
		delete [] wall;
	}

	sf->Close();
}

