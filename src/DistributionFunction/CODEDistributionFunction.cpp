/**
 * Implementation of the CODE distribution function
 * handler.
 */

#include <cmath>
#include <cstdio>
#include <string>

#include <softlib/config.h>
#include <softlib/constants.h>
#include <softlib/DistributionFunction/CODEDistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>

#ifdef INTERP_SPLINTER
#   error "SPLINTER interpolation has not been implemented in the CODEDistributionFunction class yet."
#else
#   include <gsl/gsl_interp.h>
#endif

using namespace std;

/**
 * CONSTRUCTOR
 *
 * fname: Name of file containing CODE
 *    distribution function.
 * time:  Timeslice to extract from distribution.
 */
CODEDistributionFunction::CODEDistributionFunction() { }
CODEDistributionFunction::CODEDistributionFunction(const string& fname, int time, int interptype, const string& path) {
    Load(fname, time, interptype, path);
}
CODEDistributionFunction::CODEDistributionFunction(SFile *sf, int time, int interptype, const string& path) {
    Load(sf, time, interptype, path);
}

/**
 * DESTRUCTOR
 */
CODEDistributionFunction::~CODEDistributionFunction() {
    if (this->interp_accel != nullptr)
        gsl_interp_accel_free(this->interp_accel);
    if (this->interp != nullptr) {
        for (unsigned int i = 0; i < this->nleg; i++) {
            gsl_interp_free(this->interp[i]);
        }
        delete [] this->interp;
    }
}

/**
 * Load a CODE distribution from the
 * designated file. This routine assumes
 * that the contents of the CODE output struct
 * have been saved in the file, and not the actual
 * output struct itself.
 *
 * fname:      Name of file containing CODE
 *             distribution function.
 * time:       Time index of distribution slice to load.
 * interptype: Interpolation method to use (CODEDistributionFunction::INTERPOLATION_???).
 * path:       Path in input file to load distribution from.
 */
void CODEDistributionFunction::Load(const string& fname, int time, int interptype, const string& path) {
    SFile *sf;
    sf = SFile::Create(fname, SFILE_MODE_READ);

    this->Load(sf, time, interptype, path);
    
    sf->Close();
}
void CODEDistributionFunction::Load(SFile *sf, int time, int interptype, const string& path) {
    double **tf, **ty, **tdelta, **tNxi, **tnref, **tTref;
    slibreal_t *ff, *pp;
    sfilesize_t fsize[2];
    int i, nnp, Nleg, ntimes;

    string modpath = "";
    if (!path.empty()) {
        modpath = path;
        if (path.back() != '/')
            modpath += '/';
    }

    ty = sf->GetDoubles(modpath+"y", fsize);
    if (fsize[0] == 1) nnp = fsize[1];
    else if (fsize[1] == 1) nnp = fsize[0];
    else throw SOFTLibException("Invalid size of CODE distribution function momentum grid: %llu x %llu.", fsize[0], fsize[1]);

    tdelta = sf->GetDoubles(modpath+"delta", fsize);
    if (fsize[0] != 1 || fsize[1] != 1)
        throw SOFTLibException("Invalid size of CODE distribution 'delta' variable: %llu x %llu.", fsize[0], fsize[1]);

    // Optional reference value
    if (sf->HasVariable(modpath+"nref")) {
        tnref = sf->GetDoubles(modpath+"nref", fsize);
        if (fsize[0] != 1 || fsize[1] != 1)
            throw SOFTLibException("Invalid size of CODE distribution 'nref' variable: %llu x %llu.", fsize[0], fsize[1]);
    } else {
        // This somewhat weird code is necessary due to that
        // 'tnref' must be a 2D array (with only one element)
        tnref = new double*[1];
        tnref[0] = new double[1];
        tnref[0][0] = 1;
    }

    // Optional reference value
    if (sf->HasVariable(modpath+"Tref")) {
        tTref = sf->GetDoubles(modpath+"Tref", fsize);
        if (fsize[0] != 1 || fsize[1] != 1)
            throw SOFTLibException("Invalid size of CODE distribution 'Tref' variable: %llu x %llu.", fsize[0], fsize[1]);
    } else {
        // This somewhat weird code is necessary due to that
        // 'tTref' must be a 2D array (with only one element)
        tTref = new double*[1];
        tTref[0] = new double[1];
        tTref[0][0] = 1;
    }

    tNxi = sf->GetDoubles(modpath+"Nxi", fsize);
    if (fsize[0] != 1 || fsize[1] != 1)
        throw SOFTLibException("Invalid size of CODE distribution 'Nxi' variable: %llu x %llu.", fsize[0], fsize[1]);

    Nleg = (int)tNxi[0][0];

    // Convert from y to p
    pp = new slibreal_t[nnp];
    for (i = 0; i < nnp; i++)
        pp[i] = (slibreal_t)(ty[0][i] * tdelta[0][0]);

    tf = sf->GetDoubles(modpath+"f", fsize);
    if (fsize[1] != (sfilesize_t)nnp*Nleg)
        throw SOFTLibException("Invalid size of CODE distribution 'f' variable: %llu x %llu.", fsize[0], fsize[1]);
    ntimes = fsize[0];

    /* Convert negative time index to
     * index relative to end. This way
     *   -1  =>  end
     *   -2  =>  end-1
     *   ... */
    if (time < 0) {
        time = ntimes+time;
        this->timestep = time;
    }

    if (time >= ntimes || time < 0)
        throw SOFTLibException("Invalid time index selected for CODE distribution function: %d.", time);

    // Copy to f
    ff = new slibreal_t[nnp*Nleg];
    double norm_div = 2*M_PI*ELECTRON_MASS*tTref[0][0];
    double norm     = tnref[0][0] / (sqrt(norm_div) * norm_div);
    for (i = 0; i < nnp*Nleg; i++) {
        ff[i] = (slibreal_t)(tf[time][i] * norm);
    }

    delete [] tf[0];
    delete [] tf;
    delete [] ty[0];
    delete [] ty;
    delete [] tdelta[0];
    delete [] tdelta;
    delete [] tnref[0];
    delete [] tnref;
    delete [] tTref[0];
    delete [] tTref;
    delete [] tNxi[0];
    delete [] tNxi;

    Initialize(nnp, Nleg, pp, ff, interptype);
}

/**
 * Evaluates the Nleg first Legendre
 * polynomials in x=xi.
 *
 * legmodes: Array to store Legendre modes in.
 * xi:       Point to evaluate the polynomials in.
 * Nleg      Number of Legendre polynomials to evaluate.
 */
void CODEDistributionFunction::ComputeLegendrePolynomials(
    slibreal_t *legmodes, slibreal_t xi, int Nleg
) {
    int i;
    slibreal_t n;

    legmodes[0] = 1.0;
    if (Nleg >= 1)
        legmodes[1] = xi;

    for (i = 0; i < Nleg-2; i++) {
        n = (slibreal_t)(i+1);
        legmodes[i+2] = (2.0*n+1.0)/(n+1.0) * xi * legmodes[i+1] - n/(n+1.0) * legmodes[i];
    }
}

/**
 * Destroy this object by deallocating all dynamically
 * allocated properties.
 */
void CODEDistributionFunction::Destroy() {
    if (this->f != nullptr) {
        delete [] this->f;
        this->f = nullptr;
    }
    if (this->p != nullptr) {
        delete [] this->p;
        this->p = nullptr;
    }
    if (this->legendre != nullptr) {
        delete [] this->legendre;
        this->legendre = nullptr;
    }

    if (this->interp_accel != nullptr) {
        gsl_interp_accel_free(this->interp_accel);
        this->interp_accel = nullptr;
    }
    if (this->interp != nullptr) {
        for (unsigned int i = 0; i < this->nleg; i++) {
            gsl_interp_free(this->interp[i]);
        }
        delete [] this->interp;
        this->interp = nullptr;
    }
}

/**
 * Evaluate the CODE distribution function in the given
 * point of momentum-space.
 * 
 * p:  Particle momentum.
 * xi: Particle pitch.
 */
slibreal_t CODEDistributionFunction::Eval(const slibreal_t p, const slibreal_t xi) {
    unsigned int i;
    slibreal_t s = 0.0, fval;
    slibreal_t sgn[2] = {1.0,-1.0};

    if (p > this->GetPMax() && p < this->GetPMin())
        return 0;

    if (xi != this->xileg) {
        ComputeLegendrePolynomials(this->legendre, xi, this->nleg);
        this->xileg = xi;
    }

    for (i = 0; i < this->nleg; i++) {
        fval = gsl_interp_eval(this->interp[i], this->p, this->f+(i*this->np), p, this->interp_accel);
        s += fval * this->legendre[i] * sgn[i%2];
    }

    return s;
}

/**
 * Get maximum momentum on grid.
 */
slibreal_t CODEDistributionFunction::GetPMax() {
    return ((p[np-1] > p[0]) ? p[np-1] : p[0]);
}
/**
 * Get minimum momentum on grid.
 */
slibreal_t CODEDistributionFunction::GetPMin() {
    return ((p[np-1] > p[0]) ? p[0] : p[np-1]);
}

/**
 * Initialize this distribution from p and f loaded
 * from file.
 */
void CODEDistributionFunction::Initialize(const unsigned int np, const unsigned int nleg, slibreal_t *p, slibreal_t *f, int interptype) {
    this->np = np;
    this->nleg = nleg;
    this->p = p;
    this->f = f;
    this->xileg = nan("");
    
    // Allocate space for legendre modes
    this->legendre = new slibreal_t[nleg];

    InitInterpolation(interptype);
}

/**
 * Initialize interpolation objects. Should
 * be called after f, p, np and nleg have been
 * initialized.
 *
 * interptype: Interpolation method to use (CODEDistributionFunction::INTERPOLATION_???).
 */
void CODEDistributionFunction::InitInterpolation(int interptype) {
    unsigned int i;

    if (this->f == nullptr || this->p == nullptr)
        throw SOFTLibException("CODE distribution function: f and p have not been initialized yet.");

    this->interptype = interptype;

    this->interp_accel = gsl_interp_accel_alloc();
    this->interp = new gsl_interp*[this->np];
    for (i = 0; i < this->nleg; i++) {
        switch (interptype) {
            case INTERPOLATION_LINEAR:
                this->interp[i] = gsl_interp_alloc(gsl_interp_linear, this->np);
                break;
            case INTERPOLATION_POLYNOMIAL:
                this->interp[i] = gsl_interp_alloc(gsl_interp_polynomial, this->np);
                break;
            case INTERPOLATION_CSPLINE:
                this->interp[i] = gsl_interp_alloc(gsl_interp_cspline, this->np);
                break;
            case INTERPOLATION_CSPLINE_PERIODIC:
                this->interp[i] = gsl_interp_alloc(gsl_interp_cspline_periodic, this->np);
                break;
            case INTERPOLATION_AKIMA:
                this->interp[i] = gsl_interp_alloc(gsl_interp_akima, this->np);
                break;
            case INTERPOLATION_AKIMA_PERIODIC:
                this->interp[i] = gsl_interp_alloc(gsl_interp_akima_periodic, this->np);
                break;
            case INTERPOLATION_STEFFEN:
                this->interp[i] = gsl_interp_alloc(gsl_interp_steffen, this->np);
                break;
            default:
                throw SOFTLibException("CODE distribution function: Unrecognized interpolation method requested: %d.", interptype);
        }

        gsl_interp_init(this->interp[i], this->p, this->f+(i*this->np), this->np);
    }
}

/**
 * Clone this distribution function.
 */
CODEDistributionFunction *CODEDistributionFunction::MinClone() {
    CODEDistributionFunction *cdf
        = new CODEDistributionFunction();

    cdf->Initialize(this->np, this->nleg, this->p, this->f, this->interptype);
    return cdf;
}

/**
 * Convert the CODE distribution to a numeric distribution
 * function gridded in momentum-space (in difference to
 * Legendre-mode decomposed momentum-space).
 *
 * (optional) uniformTheta: If true, generates a xi grid that is
 *    uniform in theta rather than uniform in xi = cos(theta). (default = false)
 * (optional) nnxi: Number of points in cosine of pitch angle. (default = 200)
 * (optional) interptype: Interpolation method to use for the new
 *    numeric momentum-space distribution function
 *    (NumericMomentumSpaceDistributionFunction::INTERPOLATION_???).
 */
NumericMomentumSpaceDistributionFunction *CODEDistributionFunction::ToMomentumSpace(
    bool logarithmic, bool uniformTheta, unsigned int nnxi, int interptype
) {
    slibreal_t *txi, *ff, sgn[2] = {1.0,-1.0};
    unsigned int elindex, i, j, k;
    NumericMomentumSpaceDistributionFunction *msdf =
        new NumericMomentumSpaceDistributionFunction();

    // Generate xi
    txi = new slibreal_t[nnxi];
    for (i = 0; i < nnxi; i++) {
        if (uniformTheta)
            txi[i] = cos(M_PI*((slibreal_t)(nnxi-i-1))/((slibreal_t)nnxi-1.0));
        else
            txi[i] = -1.0 + 2*((slibreal_t)i)/((slibreal_t)nnxi-1.0);
    }

    // Generate f
    ff = new slibreal_t[nnxi*np];
    for (i = 0; i < nnxi; i++) {
        ComputeLegendrePolynomials(this->legendre, txi[i], this->nleg);

        for (j = 0; j < np; j++) {
            elindex = i*np + j;
            ff[elindex] = 0.0;
            for (k = 0; k < this->nleg; k++) {
                ff[elindex] += this->f[k*np+j] * this->legendre[k] * sgn[k%2];
            }

            if (logarithmic) {
                if (ff[elindex] < 0)
                    ff[elindex] = 0.0;
            }
        }
    }

    if (logarithmic)
        msdf->InitializeLog(np, nnxi, p, txi, ff, interptype);
    else
        msdf->Initialize(np, nnxi, p, txi, ff, interptype);

    return msdf;
}

