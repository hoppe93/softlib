/**
 * Implementation of the SOFT Distribution Function
 * subclass of the Distribution Function class.
 */

#include <cstdio>
#include <cstring>
#include <string>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/SOFTDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>


using namespace std;

// Type string (used to identify distribution function files)
const string SOFTDistributionFunction::MAGIC = "distribution/soft";

SOFTDistributionFunction::SOFTDistributionFunction(const string &fname, MagneticField2D *mf, bool logarithmic, int interptype) {
    Load(fname, mf, logarithmic, interptype);
}

/**
 * Loads the SOFT distribution located in the file 'fname'
 * into this object. Initializes the interpolation objects
 * of this distribution function.
 *
 * fname:       Name of file containing SOFT distribution
 *              function to load.
 * logarithmic: If true, stores the logarithm of f and interpolates
 *              in that value instead. The returned value is f though
 *              (exp(log(f))).
 * interptype:  Type of interpolation to do (NumericDistributionFunction::INTERPOLATION_???).
 */
void SOFTDistributionFunction::Load(const string& fname, MagneticField2D *mf, bool logarithmic, int interptype) {
    struct softdf_data *dat;
    dat = SOFTDistributionFunction::__Load(fname, mf);

    if (logarithmic) {
        for (unsigned int i = 0; i < dat->nr; i++) {
            struct softmdf_data *mdf = dat->fr+i;

            this->InsertMomentumSpaceDistributionLog(
                mdf->np, mdf->nxi, dat->r[i],
                mdf->p, mdf->xi, mdf->f, interptype
            );
        }
    } else {
        for (unsigned int i = 0; i < dat->nr; i++) {
            struct softmdf_data *mdf = dat->fr+i;

            this->InsertMomentumSpaceDistribution(
                mdf->np, mdf->nxi, dat->r[i],
                mdf->p, mdf->xi, mdf->f, interptype
            );
        }
    }

    delete [] dat->r;
    delete dat;
}

/**
 * Loads the SOFT distribution function located in file 'fname'
 * and initializes a new momentum-space distribution function
 * with the loaded function at one radial location.
 *
 * fname:       Name of file containing SOFT distribution
 *              function to load.
 * mf:          Magnetic field in which the distribution
 *              function lives. Not used, but here for
 *              compatibility with other numerical
 *              distribution functions.
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
*SOFTDistributionFunction::LoadMomentumSpace(
    const string& fname, MagneticField2D *mf,
	bool logarithmic, unsigned int radindex,
    int interptype
) {
    NumericMomentumSpaceDistributionFunction *msdf
        = new NumericMomentumSpaceDistributionFunction();
    struct softdf_data *dat;
    
    dat = __Load(fname, mf);

    if (radindex >= dat->nr)
        throw SOFTLibException("Invalid radial index selected for the SOFT distribution function: %u > %u (number of radial points).", radindex, dat->nr-1);

    // Copy selected part of distribution function
    struct softmdf_data *mdf = dat->fr+radindex;
    
    if (logarithmic)
        msdf->InitializeLog(mdf->np, mdf->nxi, mdf->p, mdf->xi, mdf->f, interptype);
    else
        msdf->Initialize(mdf->np, mdf->nxi, mdf->p, mdf->xi, mdf->f, interptype);

    delete [] dat->r;
    for (unsigned int i = 0; i < dat->nr; i++) {
        if (i == radindex)
            continue;
        
        mdf = dat->fr+i;

        delete [] mdf->p;
        delete [] mdf->xi;
        delete [] mdf->f;
    }

    delete [] dat->fr;

    // p and xi are still used!
    delete dat;

    return msdf;
}

/**
 * Loads the SOFT distribution located in the file 'fname'
 * into a SOFT distribution function data structure that
 * is returned. This is a static method.
 *
 * fname: Name of file containing SOFT distribution
 *    function to load.
 *
 * RETURNS a data structure representing the SOFT data file.
 */
struct SOFTDistributionFunction::softdf_data *SOFTDistributionFunction::__Load(const string &fname, MagneticField2D *mf) {
    sfilesize_t fsize[2];
    unsigned int nnr;
    struct softdf_data *dat;

    SFile *sf = SFile::Create(fname, SFILE_MODE_READ);

    // Check if type is specified
    if (sf->HasVariable("type")) {
        string type = sf->GetString("type");
        if (type != SOFTDistributionFunction::MAGIC)
            throw SOFTLibException("The given file is not a SOFT distribution function file.");
    }

    // r grid
    double **tr = sf->GetDoubles("r", fsize);
    if (fsize[0] == 1) nnr = fsize[1];
    else if (fsize[1] == 1) nnr = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function 'r' vector: %llu x %llu.", fsize[0], fsize[1]);

    dat     = new struct softdf_data;
    dat->nr = nnr;
    dat->fr = new struct softmdf_data[nnr];

    // Load individual momentum space distribution functions
    for (unsigned int i = 0; i < nnr; i++) {
        string groupname = "r"+to_string(i)+"/";
        __LoadRadius(sf, groupname, dat->fr+i);
    }

    sf->Close();

    // Convert 'double' to 'slibreal_t' if necessary
    if (std::is_same<slibreal_t,double>::value) {
        // Not necessary, just copy
        dat->r  = tr[0];
    } else {
        // Necessary, copy
        dat->r = new slibreal_t[nnr];
        for (unsigned int i = 0; i < nnr; i++)
            dat->r[i] = (slibreal_t)tr[0][i];

        delete [] tr[0];
        delete [] tr;
    }

    // Shift radial grid (convert from normalized minor radius -> major radius)
    for (unsigned int i = 0; i < nnr; i++) {
        dat->r[i] += mf->GetMagneticAxisR();

        if (i > 0 && dat->r[i] <= dat->r[i-1])
            throw SOFTLibException("The radial grid in the SOFT distribution function is not strictly increasing.");
    }

    return dat;
}

/**
 * Load a momentum space distribution function at a
 * single radius.
 *
 * f:         SFile handle to read file with.
 * groupname: Name of group (in the file) from which the
 *            the distribution function should be read.
 *
 * RETURNS a data structure representing the SOFT data file.
 */
void SOFTDistributionFunction::__LoadRadius(
    SFile *sf, const string groupname, struct softmdf_data *dat
) {
    sfilesize_t fsize[2];
    unsigned int nnp, nnxi, dn;

    // p grid
    double **tp = sf->GetDoubles(groupname+"p", fsize);
    if (fsize[0] == 1) nnp = fsize[1];
    else if (fsize[1] == 1) nnp = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function '%sp' vector: %llu x %llu.", groupname.c_str(), fsize[0], fsize[1]);

    // xi grid
    double **txi = sf->GetDoubles(groupname+"xi", fsize);
    if (fsize[0] == 1) nnxi = fsize[1];
    else if (fsize[1] == 1) nnxi = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function '%sxi' vector: %llu x %llu.", groupname.c_str(), fsize[0], fsize[1]);

    // Distribution function
    double **ttf = sf->GetDoubles(groupname+"f", fsize);
	double *tf = ttf[0];
    delete ttf;

    if (fsize[0] != nnxi || fsize[1] != nnp)
        throw SOFTLibException("Invalid size of SOFT distribution function '%sf' vector: %llu x %llu. Must have size nxi-by-np (%llu x %llu).", groupname.c_str(), fsize[0], fsize[1], nnxi, nnp);

    dat->np = nnp;
    dat->nxi = nnxi;
    
    if (std::is_same<slibreal_t,double>::value) {
        // Not necessary, just copy
        dat->p  = tp[0];
        dat->xi = txi[0];
        dat->f  = tf;
    } else {
        // Necessary, copy
        dat->p = new slibreal_t[nnp];
        for (unsigned int i = 0; i < nnp; i++)
            dat->p[i] = (slibreal_t)tp[0][i];

        dat->xi = new slibreal_t[nnxi];
        for (unsigned int i = 0; i < nnxi; i++)
            dat->xi[i] = (slibreal_t)txi[0][i];

        dn = nnxi*nnp;
        dat->f = new slibreal_t[dn];
        for (unsigned int i = 0; i < dn; i++) {
			dat->f[i] = (slibreal_t)tf[i];
        }

        delete [] tp[0];
        delete [] tp;
        delete [] txi[0];
        delete [] txi;
        delete [] tf;
    }
}

