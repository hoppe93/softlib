/**
 * Interface for loading GO+CODE distribution functions.
 */

#include <string>
#include "softlib/DistributionFunction/CODEDistributionFunction.h"
#include "softlib/DistributionFunction/GOCODEDistributionFunction.h"
#include "softlib/MagneticField/MagneticField2D.h"

using namespace std;


// Type string (used to identify distribution function files)
const string GOCODEDistributionFunction::MAGIC = "distribution/gocode";

/**
 * Constructor.
 */
GOCODEDistributionFunction::GOCODEDistributionFunction(
    const string& fname, MagneticField2D *mf, int time, int interptype
) {
    this->Load(fname, mf, time, interptype);
}

/**
 * Destructor.
 */
GOCODEDistributionFunction::~GOCODEDistributionFunction() {
    delete [] this->code;
}

/**
 * Evaluate this distribution function.
 */
slibreal_t GOCODEDistributionFunction::Eval(
    const slibreal_t r, const slibreal_t p,
    const slibreal_t xi, const slibreal_t drift_shift
) {
    slibreal_t rho = r-drift_shift, f0, f1;

    if (!this->allowExtrapolation) {
        if (rho < this->rhomin || rho > this->rhomax)
            throw SOFTLibException("GOCODEDistributionFunction: No data avilable for particle position. rho = %e, (shifted by %e from r = %e).", rho, drift_shift, r);
    }

    unsigned int ir = (unsigned int)__FindNearestR(rho);
    slibreal_t r0, r1;

    if (ir == 0) {
        r0 = 0;
        r1 = this->rho[0];
        slibreal_t r2 = this->rho[1], f2;

        f1 = this->code[0]->Eval(p, xi);
        f2 = this->code[1]->Eval(p, xi);
        f0 = (r2*f1 - r1*f2) / (r2-r1);
    } else if (ir >= this->nr-1) {
        r0 = this->rho[this->nr-1];
        r1 = r0 + (2*r0 - this->rho[nr-2]);

        f0 = this->code[this->nr-1]->Eval(p, xi);
        f1 = 2*f0 - this->code[this->nr-2]->Eval(p, xi);
    } else {
        r0 = this->rho[ir];
        r1 = this->rho[ir+1];

        f0 = this->code[ir]->Eval(p, xi);
        f1 = this->code[ir+1]->Eval(p, xi);
    }

    return ((r1-rho)*f0 + (rho-r0)*f1)/(r1-r0);
}

/**
 * Finds the index of the point in the radial grid
 * that is closest to (and less than or equal to) r.
 *
 * r: Radial point for which to find the nearest preceeding
 *    point in the radial grid.
 * 
 * RETURNS the index of the point nearest preceeding r
 * in 'this->rho'.
 */
int GOCODEDistributionFunction::__FindNearestR(const slibreal_t r) {
    bool ascnd;
    int il, im, iu, nnr = (int)nr;

    if (nr == 1) return 0;
    ascnd = (this->rho[1] > this->rho[0]);
    
    il = 0;
    iu = nr-1;
    while (iu-il > 1) {
        im = (iu+il) >> 1;
        if ((r >= this->rho[im]) == ascnd)
            il = im;
        else
            iu = im;
    }

    if (il > nnr-1)
        return nr-1;
    else
        return il;
}

/**
 * Initialize this GO+CODE distribution function from a set
 * of radii and CODE distribution functions.
 */
void GOCODEDistributionFunction::Initialize(unsigned int nr, slibreal_t *rho, CODEDistributionFunction **cdf) {
    this->nr   = nr;
    this->rho  = rho;
    this->code = cdf;
}

/**
 * Load the GO+CODE distribution function stored in
 * the file with the given name.
 *
 * fname:       Name of file containing GO+CODE distribution function
 *              to load.
 * time:        Time index of distribution function to load.
 * logarithmic: If true, stores the logarithm of f and interpolates
 *              in that value instead. The returned value is f though
 *              (exp(log(f))).
 * interptype:  Type of interpolation to do (NumericDistributionFunction::INTERPOLATION_???).
 */
void GOCODEDistributionFunction::Load(const string& fname, MagneticField2D *mf, int time, int interptype) {
    SFile *sf = SFile::Create(fname, SFILE_MODE_READ);

    sfilesize_t fsize[2];
    unsigned int nnr;

    // Check if type is specified
    if (sf->HasVariable("type")) {
        string type = sf->GetString("type");
        if (type != GOCODEDistributionFunction::MAGIC)
            throw SOFTLibException("The given file is not a SOFT distribution function file.");
    }

    // ///////////
    // r grid
    // ///////////
    double **tr = sf->GetDoubles("r", fsize);
    if (fsize[0] == 1) nnr = fsize[1];
    else if (fsize[1] == 1) nnr = fsize[0];
    else throw SOFTLibException("Invalid size of GO+CODE distribution function 'x_CODE' vector: %llu x %llu.", fsize[0], fsize[1]);

    this->nr = nnr;
    this->rho = new slibreal_t[nnr];
    for (unsigned int i = 0; i < nnr; i++)
        this->rho[i] = tr[0][i];

    delete [] tr[0];
    delete [] tr;

    // Shift radii by magnetic field major radius
    // and verify that the array is ordered
    for (unsigned int i = 0; i < nnr; i++) {
        this->rho[i] += mf->GetMagneticAxisR();
        
        if (i > 0 && this->rho[i] <= this->rho[i-1])
            throw SOFTLibException("The radial grid in the GO+CODE distribution function is not strictly increasing.");
    }

    this->rhomin = this->rho[0];
    this->rhomax = this->rho[this->nr-1];

    // ///////////
    // nt
    // ///////////
    int ntimes = (int)sf->GetScalar("nt");
    
    // Convert negative time idnex to
    // index relative to end. This way,
    //   -1  =>  end
    //   -2  =>  end-1
    //   ...
    if (time < 0) {
        time = ntimes+time;
    }

    if (time >= ntimes || time < 0)
        throw SOFTLibException("Invalid time index select for GO+CODE distribution function: %d.", time);

    // ////////////////////////////
    // CODE distribution functions
    // ////////////////////////////
    this->code = new CODEDistributionFunction*[nnr];
    
    string tname = "t"+to_string(time)+"/";
    for (unsigned int i = 0; i < nnr; i++) {
        string groupname = tname+"r"+to_string(i)+"/";

        this->code[i] = new CODEDistributionFunction(sf, 0, interptype, groupname);
    }
}

/**
 * Clone this distribution function.
 */
GOCODEDistributionFunction *GOCODEDistributionFunction::MinClone() {
    GOCODEDistributionFunction *gcdf = new GOCODEDistributionFunction();

    CODEDistributionFunction **cdf = new CODEDistributionFunction*[this->nr];
    slibreal_t *newr = new slibreal_t[this->nr];

    for (unsigned int i = 0; i < this->nr; i++) {
        cdf[i]  = this->code[i]->MinClone();
        newr[i] = this->rho[i];
    }

    gcdf->Initialize(this->nr, newr, cdf);

    return gcdf;
}

