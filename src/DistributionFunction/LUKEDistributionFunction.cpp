/**
 * Implementation of the LUKE Distribution Function
 * subclass of the Distribution Function class.
 */

#include <cstdio>
#include <cstring>
#include <string>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/LUKEDistributionFunction.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>

using namespace std;

LUKEDistributionFunction::LUKEDistributionFunction(const string &fname, MagneticField2D *magfield, bool logarithmic, int interptype) {
	this->magfield = magfield;

    Load(fname, logarithmic, interptype);
}

/**
 * Loads the LUKE distribution located in the file 'fname'
 * into this object. Initializes the interpolation objects
 * of this distribution function.
 *
 * fname:       Name of file containing LUKE distribution
 *              function to load.
 * logarithmic: If true, stores the logarithm of f and interpolates
 *              in that value instead. The returned value is f though
 *              (exp(log(f))).
 * interptype:  Type of interpolation to do (NumericDistributionFunction::INTERPOLATION_???).
 */
void LUKEDistributionFunction::Load(const string& fname, bool logarithmic, int interptype) {
    struct lukedf_data *dat;
    dat = LUKEDistributionFunction::__Load(fname, this->magfield);

    if (logarithmic) {
        this->InitializeLog(
            dat->nr, dat->np, dat->nxi,
            dat->r, dat->p, dat->xi, dat->f,
            interptype
        );
    } else {
        this->Initialize(
            dat->nr, dat->np, dat->nxi,
            dat->r, dat->p, dat->xi, dat->f,
            interptype
        );
    }

    delete dat;
}

/**
 * Loads the LUKE distribution function located in file 'fname'
 * and initializes a new momentum-space distribution function
 * with the loaded function at one radial location.
 *
 * fname:       Name of file containing LUKE distribution
 *              function to load.
 * mf:          Magnetic field in which the distribution
 *              function lives.
 * logarithmic: If true, stores the logarithm of f instead of
 *              f itself and interpolates in the logarithm.
 *              The value returned by 'Eval()' later is f though.
 * radindex:    Index of the radial point corresponding to
 *              the momentum-space distribution function to
 *              extract (default = 0).
 * interptype:  Type of interpolation to do (NumericMomentumSpaceDistributionFunction::INTERPOLATION_???).
 *
 * RETURNS a newly allocated NumericMomentumSpaceDistributionFunction
 * object initialized with the specified momentum-space
 * distribution function.
 */
NumericMomentumSpaceDistributionFunction
*LUKEDistributionFunction::LoadMomentumSpace(
    const string& fname, MagneticField2D *mf,
	bool logarithmic, unsigned int radindex, int interptype
) {
    NumericMomentumSpaceDistributionFunction *msdf
        = new NumericMomentumSpaceDistributionFunction();
    struct lukedf_data *dat;
    
    dat = __Load(fname, mf);

    if (radindex >= dat->nr)
        throw SOFTLibException("Invalid radial index selected for the LUKE distribution function: %u > %u (number of radial points).", radindex, dat->nr-1);

    // Copy selected part of distribution function
    slibreal_t *tf = new slibreal_t[dat->np*dat->nxi];
    memcpy(tf, dat->f+radindex*(dat->np*dat->nxi), sizeof(slibreal_t)*dat->np*dat->nxi);
    
    if (logarithmic)
        msdf->InitializeLog(dat->np, dat->nxi, dat->p, dat->xi, tf, interptype);
    else
        msdf->Initialize(dat->np, dat->nxi, dat->p, dat->xi, tf, interptype);

    delete [] dat->r;
    delete [] dat->f;
    // p and xi are still used!
    delete dat;

    return msdf;
}

/**
 * Loads the LUKE distribution located in the file 'fname'
 * into a LUKE distribution function data structure that
 * is returned. This is a static method.
 *
 * fname: Name of file containing LUKE distribution
 *    function to load.
 *
 * RETURNS a data structure representing the LUKE data file.
 */
struct LUKEDistributionFunction::lukedf_data *LUKEDistributionFunction::__Load(const string &fname, MagneticField2D *mf) {
    string punits;
    sfilesize_t fsize[2];
    unsigned int i, nnr, nnp, nnxi, dn;
    struct lukedf_data *dat;

    SFile *sf = SFile::Create(fname, SFILE_MODE_READ);

    // r grid
    double **tr = sf->GetDoubles("xrhoG", fsize);
    if (fsize[0] == 1) nnr = fsize[1];
    else if (fsize[1] == 1) nnr = fsize[0];
    else throw SOFTLibException("Invalid size of LUKE distribution function 'xrhoG' vector: %llu x %llu.", fsize[0], fsize[1]);

	slibreal_t
		minr = mf->GetMagneticAxisR(),
		maxr = mf->GetMaxRadius();
	
	// Transform radial grid from normalized -> major radius
	for (i = 0; i < nnr; i++)
		tr[0][i] = minr + tr[0][i] * (maxr - minr);

    // p grid
    double **tp = sf->GetDoubles("pn", fsize);
    if (fsize[0] == 1) nnp = fsize[1];
    else if (fsize[1] == 1) nnp = fsize[0];
    else throw SOFTLibException("Invalid size of LUKE distribution function 'pn' vector: %llu x %llu.", fsize[0], fsize[1]);

	double **tbetath = sf->GetDoubles("betath_ref", fsize);
	if (fsize[0] != 1 || fsize[1] != 1)
		throw SOFTLibException("Invalid size of LUKE distribution function 'betath_ref' variable: %llu x %llu.", fsize[0], fsize[1]);
	
	// Normalize momentum
	for (i = 0; i < nnp; i++)
		tp[0][i] *= **tbetath;

    // xi grid
    double **txi = sf->GetDoubles("mhu", fsize);
    if (fsize[0] == 1) nnxi = fsize[1];
    else if (fsize[1] == 1) nnxi = fsize[0];
    else throw SOFTLibException("Invalid size of LUKE distribution function 'mhu' vector: %llu x %llu.", fsize[0], fsize[1]);

    // Distribution function
	sfilesize_t ndims, fsize3[3];
    double *tf = sf->GetMultiArray_linear("f", 3, ndims, fsize3);

	if (tf == nullptr || ndims != 3)
		throw SOFTLibException("Invalid number of dimensions in LUKE distribution function 'f': %llu.", ndims);
	if (fsize3[0] != nnr || fsize3[1] != nnxi || fsize3[2] != nnp)
		throw SOFTLibException(
			"Invalid size of LUKE distribution function 'f': %llu x %llu x %llu. Expected %llu x %llu x %llu.",
			fsize3[0], fsize3[1], fsize3[2], nnr, nnxi, nnp
		);

    sf->Close();

    //this->Initialize(nnp, nnxi, tp[0], txi[0], tf[0]);
    dat     = new struct lukedf_data;
    dat->nr  = nnr;
    dat->np  = nnp;
    dat->nxi = nnxi;

    // Convert 'double' to 'slibreal_t' if necessary
    if (std::is_same<slibreal_t,double>::value) {
        // Not necessary, just copy
        dat->r  = tr[0];
        dat->p  = tp[0];
        dat->xi = txi[0];
        dat->f  = tf;
    } else {
        // Necessary, copy
        dat->r = new slibreal_t[nnr];
        for (i = 0; i < nnr; i++)
            dat->r[i] = (slibreal_t)tr[0][i];

        dat->p = new slibreal_t[nnp];
        for (i = 0; i < nnp; i++)
            dat->p[i] = (slibreal_t)tp[0][i];

        dat->xi = new slibreal_t[nnxi];
        for (i = 0; i < nnxi; i++)
            dat->xi[i] = (slibreal_t)txi[0][i];

        dn = nnp*nnxi*nnr;
        dat->f    = new slibreal_t[dn];
        for (i = 0; i < dn; i++) {
			dat->f[i] = (slibreal_t)tf[i];
        }

        delete [] tr[0];
        delete [] tr;
        delete [] tp[0];
        delete [] tp;
        delete [] txi[0];
        delete [] txi;
        delete [] tf;
    }

    return dat;
}

