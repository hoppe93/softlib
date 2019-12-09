/**
 * Wrapper for loading a LUKE magnetic equilibrium file.
 *
 * LUKE magnetic fields are contained in MATLAB files prefixed
 * with 'EQUIL_' and with extension '.mat'. The file contains
 * one struct named 'equil', which in turn contains a number
 * of fields.
 *
 * The fields of most interest to this implementation are
 *   
 *   equil/psi_apRp     -- Poloidal flux
 *   equil/theta		-- Poloidal angle
 *   equil/ptx
 *   equil/pty          -- X- and Y coordinates of poloidal
 *                         grid, on which the magnetic field
 *                         is given. (NOTE: relative to magnetic
 *                         axis location, so that (x,y) = (0,0)
 *                         is where the magnetic axis is)
 *   equil/ptBx         -- Radial component of magnetic field.
 *   equil/ptBy         -- Vertical component of magnetic field.
 *   equil/ptBPHI       -- Toroidal component of magnetic field.
 *                         (NOTE: the magnetic field components are
 *                         given on the grid defined by 'psi_apRp'
 *                         and 'theta', the poloidal flux and angle.)
 *   equil/Rp           -- Magnetic axis radial position (in global
 *                         coordinates)
 *   equil/Zp           -- Magnetic axis vertical position (in global
 *                         coordinates)
 *   
 * For mapping from R-Z coordinates to flux coordinates,
 * we also need
 *
 *   equil/eqdsk/x
 *   equil/eqdsk/y      -- X and Y coordinates of the grid on
 *                         which the poloidal flux is given.
 *   equil/eqdsk/xypsin -- Normalized poloidal flux, given on the
 *                         grid defined by 'x' and 'y'.
 *
 * To get the magnetic field in the format we need, we first
 * need to create a regular R-Z grid, then interpolate in the
 * xypsin array to get the corresponding poloidal flux value.
 */

#include <string>
#include <softlib/MagneticField/MagneticFieldLUKE.h>
#include <softlib/SFileException.h>

using namespace std;

MagneticFieldLUKE::MagneticFieldLUKE(const string& name) {
	Load(name);
	InitInterpolation();
}
MagneticFieldLUKE::MagneticFieldLUKE(const string& name, enum sfile_type type) {
	Load(name, type);
	InitInterpolation();
}
MagneticFieldLUKE::MagneticFieldLUKE(const string& name, const string& wallname) {
	Load(name, wallname);
	InitInterpolation();
}
MagneticFieldLUKE::MagneticFieldLUKE(const string& name, enum sfile_type type, const string& wallname) {
	Load(name, type, wallname);
	InitInterpolation();
}
MagneticFieldLUKE::MagneticFieldLUKE(const string& name, enum sfile_type type, const string& wallname, enum sfile_type walltype) {
	Load(name, type, wallname, walltype);
	InitInterpolation();
}

/**
 * Loads the LUKE magnetic equilibrium located in the file
 * named 'name' assuming it is of the type 'type'. Optionally,
 * a file containing the wall of the tokamak can be specified
 * as 'wallname'. If 'wallname' is an empty string, no wall
 * is loaded.
 *
 * name:     Name of file containing equilibrium data.
 * type:     Type of file 'name'.
 * wallname: Name of file containing tokamak wall to load.
 * walltype: Type of file 'wallname'.
 */
void MagneticFieldLUKE::Load(const string& name) {
	Load(name, SFile::TypeOfFile(name));
}
void MagneticFieldLUKE::Load(const string& name, const string& wallname) {
	Load(name, SFile::TypeOfFile(name), wallname, SFile::TypeOfFile(wallname));
}
void MagneticFieldLUKE::Load(const string& name, enum sfile_type type) {
	Load(name, type, "", SFILE_TYPE_UNKNOWN);
}
void MagneticFieldLUKE::Load(const string& name, enum sfile_type type, const string& wallname) {
	Load(name, type, wallname, SFile::TypeOfFile(wallname));
}
void MagneticFieldLUKE::Load(const string& name, enum sfile_type type, const string& wallname, enum sfile_type walltype) {
	sfilesize_t fs[2];

	SFile *sf = SFile::Create(name, SFILE_MODE_READ, type);
	
	this->name = sf->GetString("equil/id");
	this->description = this->name;

	double _maxis[2];
	_maxis[0] = sf->GetScalar("equil/Rp");
	_maxis[1] = sf->GetScalar("equil/Zp");
	this->magnetic_axis[0] = (slibreal_t)_maxis[0];
	this->magnetic_axis[1] = (slibreal_t)_maxis[1];

	double *_xgpsi    = sf->GetList("equil/eqdsk/x", fs);
	sfilesize_t nxgpsi = fs[0];

	double *_ygpsi = sf->GetList("equil/eqdsk/y", fs);
	sfilesize_t nygpsi = fs[0];

	double **_Psi 	= sf->GetDoubles("equil/eqdsk/xypsin", fs);
	if (nygpsi != fs[0] || nxgpsi != fs[1])
		throw SFileException(
			"%s: The normalized flux grid has the wrong dimensions: %zu x %zu. Expected %zu x %zu.", 
			name.c_str(), fs[0], fs[1], nygpsi, nxgpsi
		);
	
	// Poloidal flux on LCFS
	double psia = sf->GetScalar("equil/eqdsk/psia");
	
	// Build poloidal angle grid
	double **_theta = new double*[nygpsi];
	_theta[0] = new double[nxgpsi*nygpsi];
	for (sfilesize_t i = 1; i < nygpsi; i++)
		_theta[i] = _theta[i-1] + nxgpsi;
	
	for (sfilesize_t i = 0; i < nxgpsi; i++) {
		for (sfilesize_t j = 0; j < nygpsi; j++) {
			_theta[j][i] = atan2(_ygpsi[j], _xgpsi[i]);

			if (_theta[j][i] < 0)
				_theta[j][i] = 2.0*M_PI + _theta[j][i];
		}
	}

	// Set R/Z grid
	this->nr = nxgpsi;
	this->nz = nygpsi;
	this->R = new slibreal_t[this->nr];
	this->Z = new slibreal_t[this->nz];

	for (sfilesize_t i = 0; i < this->nr; i++)
		this->R[i] = _xgpsi[i] + _maxis[0];
	for (sfilesize_t i = 0; i < this->nz; i++)
		this->Z[i] = _ygpsi[i] + _maxis[1];
	
	// Load magnetic fields
	// Psi grid (on which MF is given)
	double *_psig = sf->GetList("equil/psi_apRp", fs);
	sfilesize_t npsig = fs[0];

	// GSL requires "strictly increasing" x/y-vectors,
	// so we need to take the absolute value of every element.
	// This is only used for interpolating from psi/theta to R/Z,
	// and does not enter anywhere else.
	// While we're at it, we also normalize to the value on the
	// LCFS, so that we get the normalized flux.
	for (unsigned int i = 0; i < npsig; i++)
		_psig[i] = fabs(_psig[i] / _psig[npsig-1]);
	
	// Theta grid (on which MF is given)
	double *_tsig = sf->GetList("equil/theta", fs);
	sfilesize_t ntsig = fs[0];

	double **_Br = sf->GetDoubles("equil/ptBx", fs);
	if (npsig != fs[1] || ntsig != fs[0])
		throw SFileException(
			"%s: Invalid dimensions of magnetic field variable 'ptBx': %zu x %zu. Expected %zu x %zu.",
			name.c_str(), fs[0], fs[1], ntsig, npsig
		);
	
	double **_Bz = sf->GetDoubles("equil/ptBy", fs);
	if (npsig != fs[1] || ntsig != fs[0])
		throw SFileException(
			"%s: Invalid dimensions of magnetic field variable 'ptBy': %zu x %zu. Expected %zu x %zu.",
			name.c_str(), fs[0], fs[1], ntsig, npsig
		);
	
	double **_Bt = sf->GetDoubles("equil/ptBPHI", fs);
	if (npsig != fs[1] || ntsig != fs[0])
		throw SFileException(
			"%s: Invalid dimensions of magnetic field variable 'ptBPHI': %zu x %zu. Expected %zu x %zu.",
			name.c_str(), fs[0], fs[1], ntsig, npsig
		);
	
	// Setup interpolation
	// Bicubic interpolation seems to become unstable for these types of
	// magnetic fields, and so we are required to use bilinear interpolation instead.
	/*gsl_interp2d *wspaceR = gsl_interp2d_alloc(gsl_interp2d_bicubic, npsig, ntsig);
	gsl_interp2d *wspaceZ = gsl_interp2d_alloc(gsl_interp2d_bicubic, npsig, ntsig);
	gsl_interp2d *wspaceT = gsl_interp2d_alloc(gsl_interp2d_bicubic, npsig, ntsig);*/

	gsl_interp2d *wspaceR = gsl_interp2d_alloc(gsl_interp2d_bilinear, npsig, ntsig);
	gsl_interp2d *wspaceZ = gsl_interp2d_alloc(gsl_interp2d_bilinear, npsig, ntsig);
	gsl_interp2d *wspaceT = gsl_interp2d_alloc(gsl_interp2d_bilinear, npsig, ntsig);

	gsl_interp2d_init(wspaceR, _psig, _tsig, _Br[0], npsig, ntsig);
	gsl_interp2d_init(wspaceZ, _psig, _tsig, _Bz[0], npsig, ntsig);
	gsl_interp2d_init(wspaceT, _psig, _tsig, _Bt[0], npsig, ntsig);

	gsl_interp_accel
		*raccR = gsl_interp_accel_alloc(),
		*raccZ = gsl_interp_accel_alloc(),
		*raccT = gsl_interp_accel_alloc(),
		*zaccR = gsl_interp_accel_alloc(),
		*zaccZ = gsl_interp_accel_alloc(),
		*zaccT = gsl_interp_accel_alloc();

	// Interpolate magnetic fields onto R/Z grid
	this->Br   = new slibreal_t[this->nr*this->nz];
	this->Bphi = new slibreal_t[this->nr*this->nz];
	this->Bz   = new slibreal_t[this->nr*this->nz];
	this->Psi  = new slibreal_t[this->nr*this->nz];

	for (unsigned int j = 0; j < this->nz; j++) {
		for (unsigned int i = 0; i < this->nr; i++) { 
			slibreal_t psival          = _Psi[j][i];
			this->Psi[i + j*this->nr]  = psival * psia;

			this->Br[i + j*this->nr]   = gsl_interp2d_eval_extrap(wspaceR, _psig, _tsig, _Br[0], fabs(psival), _theta[j][i], raccR, zaccR);
			this->Bphi[i + j*this->nr] = gsl_interp2d_eval_extrap(wspaceR, _psig, _tsig, _Bt[0], fabs(psival), _theta[j][i], raccT, zaccT);
			this->Bz[i + j*this->nr]   = gsl_interp2d_eval_extrap(wspaceR, _psig, _tsig, _Bz[0], fabs(psival), _theta[j][i], raccZ, zaccZ);
		}
	}

	this->hasFluxCoordinates = true;

	// Get separatrix/limiter
	double **_ptx = sf->GetDoubles("equil/ptx", fs);
	if (fs[0] != ntsig || fs[1] != npsig)
		throw SFileException(
			"%s: Invalid dimensions of variable 'ptx': %zu x %zu. Expected %zu x %zu.",
			name.c_str(), fs[0], fs[1], ntsig, npsig
		);

	double **_pty = sf->GetDoubles("equil/pty", fs);
	if (fs[0] != ntsig || fs[1] != npsig)
		throw SFileException(
			"%s: Invalid dimensions of variable 'pty': %zu x %zu. Expected %zu x %zu.",
			name.c_str(), fs[0], fs[1], ntsig, npsig
		);
	
	this->nsep = ntsig;
	this->rsep = new slibreal_t[ntsig];
	this->zsep = new slibreal_t[ntsig];

	const unsigned int sepi = npsig-1;
	for (unsigned int i = 0; i < ntsig; i++) {
		this->rsep[i] = (slibreal_t)_ptx[i][sepi] + _maxis[0];
		this->zsep[i] = (slibreal_t)_pty[i][sepi] + _maxis[1];
	}

	gsl_interp2d_free(wspaceT);
	gsl_interp2d_free(wspaceZ);
	gsl_interp2d_free(wspaceR);

	gsl_interp_accel_free(raccR);
	gsl_interp_accel_free(raccT);
	gsl_interp_accel_free(raccZ);

	gsl_interp_accel_free(zaccR);
	gsl_interp_accel_free(zaccT);
	gsl_interp_accel_free(zaccZ);

	sf->Close();

	delete [] _xgpsi;
	delete [] _ygpsi;
	delete [] _Psi[0], delete [] _Psi;
	delete [] _theta[0], delete [] _theta;
	delete [] _psig, delete [] _tsig;

	delete [] _Br[0], delete [] _Br;
	delete [] _Bt[0], delete [] _Bt;
	delete [] _Bz[0], delete [] _Bz;

	delete [] _pty, delete [] _ptx;

	LoadWall(wallname, walltype);
}

/**
 * Loads tokamak wall from the file named 'wallname'. The
 * type of the file can also be specified, and if not, it
 * is assumed from the file name extension.
 *
 * wallname: Name of file containing wall to load.
 * walltype: Type of file containing wall to load.
 */
void MagneticFieldLUKE::LoadWall(const string& wallname) {
	LoadWall(wallname, SFile::TypeOfFile(wallname));
}
void MagneticFieldLUKE::LoadWall(const string& wallname, enum sfile_type walltype) {
	// Load wall (if given)
	this->nwall = 0;
	this->rwall = nullptr;
	this->zwall = nullptr;

	if (!wallname.empty()) {
		SFile *sf = SFile::Create(wallname, SFILE_MODE_READ, walltype);
		sfilesize_t fs[2];
		
		double **_wall = sf->GetDoubles("wall", fs);

		if (fs[0] == 2) {
			this->nwall = fs[1];

			this->rwall = new slibreal_t[this->nwall];
			this->zwall = new slibreal_t[this->nwall];

			for (unsigned int i = 0; i < this->nwall; i++) {
				this->rwall[i] = _wall[0][i];
				this->zwall[i] = _wall[1][i];
			}
		} else if (fs[1] == 2) {
			this->nwall = fs[0];

			this->rwall = new slibreal_t[this->nwall];
			this->zwall = new slibreal_t[this->nwall];

			for (unsigned int i = 0; i < this->nwall; i++) {
				this->rwall[i] = _wall[i][0];
				this->zwall[i] = _wall[i][1];
			}
		} else
			throw SFileException("%s: Invalid dimensions of wall: %zu x %zu. Expected 2 x many.",
				wallname.c_str(), fs[0], fs[1]
			);

		delete [] _wall[0];
		delete [] _wall;
	}

	if (this->rwall != nullptr && this->zwall != nullptr)
		SetDomain(this->rwall, this->zwall, this->nwall);
	else
		SetDomain(this->rsep, this->zsep, this->nsep);

    SetSeparatrix(this->rsep, this->zsep, this->nsep);
}

