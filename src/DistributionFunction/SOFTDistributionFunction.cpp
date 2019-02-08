/**
 * Implementation of the SOFT Distribution Function
 * subclass of the Distribution Function class.
 */

#include <cstdio>
#include <cstring>
#include <string>

#include <softlib/constants.h>
#include <softlib/DistributionFunction/SOFTDistributionFunction.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>

using namespace std;

SOFTDistributionFunction::SOFTDistributionFunction(const string &fname, bool logarithmic, int interptype) {
    Load(fname, logarithmic, interptype);
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
void SOFTDistributionFunction::Load(const string& fname, bool logarithmic, int interptype) {
    struct softdf_data *dat;
    dat = SOFTDistributionFunction::__Load(fname);

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
 * Loads the SOFT distribution function located in file 'fname'
 * and initializes a new momentum-space distribution function
 * with the loaded function at one radial location.
 *
 * fname:       Name of file containing SOFT distribution
 *              function to load.
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
    const string& fname, bool logarithmic, unsigned int radindex,
    int interptype
) {
    NumericMomentumSpaceDistributionFunction *msdf
        = new NumericMomentumSpaceDistributionFunction();
    struct softdf_data *dat;
    
    dat = __Load(fname);

    if (radindex >= dat->nr)
        throw SOFTLibException("Invalid radial index selected for the SOFT distribution function: %u > %u (number of radial points).", radindex, dat->nr-1);

    // Copy selected part of distribution function
    slibreal_t *tf = new slibreal_t[dat->np*dat->nxi];
    memcpy(tf, dat->f[radindex], sizeof(slibreal_t)*dat->np*dat->nxi);
    
    if (logarithmic)
        msdf->InitializeLog(dat->np, dat->nxi, dat->p, dat->xi, tf, interptype);
    else
        msdf->Initialize(dat->np, dat->nxi, dat->p, dat->xi, tf, interptype);

    delete [] dat->r;
    delete [] dat->f[0];
    delete [] dat->f;
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
struct softdf_data *SOFTDistributionFunction::__Load(const string &fname) {
    SFile *sf;
    double **tf, **tr, **tp, **txi;
    string punits;
    sfilesize_t fsize[2];
    slibreal_t normf;
    unsigned int i, j, nnr, nnp, nnxi, dn;
    struct softdf_data *dat;

    sf = SFile::Create(fname, SFILE_MODE_READ);

    // r grid
    tr = sf->GetDoubles("r", fsize);
    if (fsize[0] == 1) nnr = fsize[1];
    else if (fsize[1] == 1) nnr = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function 'r' vector: %llu x %llu.", fsize[0], fsize[1]);

    // p grid
    tp = sf->GetDoubles("p", fsize);
    if (fsize[0] == 1) nnp = fsize[1];
    else if (fsize[1] == 1) nnp = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function 'p' vector: %llu x %llu.", fsize[0], fsize[1]);

    // xi grid
    txi = sf->GetDoubles("xi", fsize);
    if (fsize[0] == 1) nnxi = fsize[1];
    else if (fsize[1] == 1) nnxi = fsize[0];
    else throw SOFTLibException("Invalid size of SOFT distribution function 'xi' vector: %llu x %llu.", fsize[0], fsize[1]);

    // Distribution function
    tf = sf->GetDoubles("f", fsize);
    if (fsize[0] != nnr || fsize[1] != (nnp*nnxi))
        throw SOFTLibException("Invalid size of SOFT distribution function 'f' vector: %llu x %llu.", fsize[0], fsize[1]);

    // f(p,xi0)
    if (!Verify("fp0", tf, nnp, 1, sf))
        throw SOFTLibException("Verification of distribution function failed in 'p' dimension.");
    if (!Verify("fxi0", tf, nnxi, nnp, sf))
        throw SOFTLibException("Verification of distribution function failed in 'xi' dimension.");
    if (!Verify("fr0", tf, nnr, nnp*nnxi, sf))
        throw SOFTLibException("Verification of distribution function failed in 'r' dimension.");

    punits = sf->GetString("punits");

    sf->Close();

    // Normalize if necessary
    if (punits == "ev") {
        normf = ELECTRON_MASS_EV;
    } else if (punits == "normalized") {
        normf = 1.0;
    } else if (punits == "si") {
        normf = ELECTRON_MASS * LIGHTSPEED;
    } else
        throw SOFTLibException("Unrecognized momentum units of SOFT distribution function: '%s'.", punits.c_str());

    if (normf != 1.0)
        for (i = 0; i < nnp; i++)
            tp[0][i] /= normf;
    
    //this->Initialize(nnp, nnxi, tp[0], txi[0], tf[0]);
    dat     = new struct softdf_data;
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

        dn = nnp*nnxi;
        dat->f    = new slibreal_t*[nnr];
        dat->f[0] = new slibreal_t[nnr*dn];
        for (i = 0; i < nnr; i++) {
            if (i > 0) dat->f[i] = dat->f[i-1] + dn;

            for (j = 0; j < dn; j++) {
                dat->f[i][j] = (slibreal_t)tf[i][j];
            }
        }

        delete [] tr[0];
        delete [] tr;
        delete [] tp[0];
        delete [] tp;
        delete [] txi[0];
        delete [] txi;
        delete [] tf[0];
        delete [] tf;
    }

    return dat;
}

/**
 * Load the verification vector named 'var' from the file
 * represented by 'sf' and make sure that the vector matrix
 * 'v' matches the verification.
 *
 * The matrix v is verified by checking that
 *
 *   |v[0][step*0] - var[0]| < eps
 *   |v[0][step*1] - var[1]| < eps
 *   ...
 *   |v[0][step*n(var)] - var[n(var)]| < eps
 *
 * var:  Name of verification vector to load. The structure of
 *       the vector in file is expected to be 1-by-(<=n) or
 *       (<=n)-by-1.
 * v:    Matrix to verify (i.e. the distribution function).
 * n:    Number of elements in dimension to verify.
 * step: Step to take to get to the "next" element in v.
 * sf:   File to load verification vector from.
 *
 * RETURNS true if the matrix 'v' was successfully verified
 * against 'var'. If 'var' does not exist, also returns without
 * any further checking. Returns 'false' in all other cases.
 */
bool SOFTDistributionFunction::Verify(const string& var, double **v, sfilesize_t n, sfilesize_t step, SFile *sf) {
    double **f0 = nullptr, err;
    sfilesize_t nf0, fsize[2], i;

    try {
        f0 = sf->GetDoubles(var, fsize);
    } catch (SOFTLibException& e) {
        return true;
    }

    if (fsize[0] == 1 && fsize[1] <= n)
        nf0 = fsize[1];
    else if (fsize[1] == 1 && fsize[0] <= n)
        nf0 = fsize[0];
    else
        throw SOFTLibException("Invalid size of SOFT distribution function verification vector '%s': %llu x %llu.", var.c_str(), fsize[0], fsize[1]);

    for (i = 0; i < nf0; i++) {
        err = fabs(f0[0][i] - v[0][i*step]);

        if (f0[0][i] == 0.0) {
            if (err > DBL_EPSILON)
                return false;
        } else if (err/fabs(f0[0][i]) > DBL_EPSILON)
            return false;
    }

    return true;
}

