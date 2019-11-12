/**
 * Implementation of the 'NumericDistributionFunction' class,
 * which represents a numerical 3-D (1D+2P) distribution
 * function.
 */

#include <cmath>
#include <limits>

#include <softlib/config.h>
#include <softlib/DistributionFunction/NumericDistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>
#include <softlib/SOFTLibException.h>

/**
 * Evaluate the distribution function in the given
 * point, taking orbit drift shift into account (if any).
 *
 * r:           Radial location of the particle.
 * p:           Magnitude of momentum of the particle.
 * xi:          Pitch of the particle.
 * drift_shift: Radial drift shift of the drift surface.
 */
slibreal_t NumericDistributionFunction::Eval(
    const slibreal_t r, const slibreal_t p,
    const slibreal_t xi, const slibreal_t drift_shift
) {
    slibreal_t rho = r-drift_shift, f0, f1;
    unsigned int ir;

	if (!this->allowExtrapolation) {
		if (rho < this->rmin || rho > this->rmax)
			throw SOFTLibException("Numeric distribution function: No data available for particle position. rho = %e, (shifted by %e from r = %e).", rho, drift_shift, r);
	}

    ir = (unsigned int)__FindNearestR(rho);
    slibreal_t r0, r1;

    if (ir == 0) {
        r0 = 0;
        r1 = this->r[0];
        slibreal_t r2 = this->r[1], f2;

        f1 = this->msdf[0].Eval(p, xi);
        f2 = this->msdf[1].Eval(p, xi);
        f0 = (r2*f1 - r1*f2) / (r2-r1);
    } else if (ir >= this->nr-1) {
        r0 = this->r[this->nr-1];
        r1 = r0 + (2*r0 - this->r[nr-2]);

        f0 = this->msdf[this->nr-1].Eval(p, xi);
        f1 = 2*f0 - this->msdf[this->nr-2].Eval(p, xi);
    } else {
        r0 = this->r[ir];
        r1 = this->r[ir+1];

        f0 = this->msdf[ir].Eval(p, xi);
        f1 = this->msdf[ir+1].Eval(p, xi);
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
 * in 'this->r'.
 */
int NumericDistributionFunction::__FindNearestR(const slibreal_t r) {
    bool ascnd;
    int il, im, iu, nnr = (int)nr;

    if (nr == 1) return 0;
    ascnd = (this->r[1] > this->r[0]);
    
    il = 0;
    iu = nr-1;
    while (iu-il > 1) {
        im = (iu+il) >> 1;
        if ((r >= this->r[im]) == ascnd)
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
 * Initializes the numeric distribution function
 * with the given components.
 *
 * nr:         Number of grid points in r.
 * np:         Number of grid points in p.
 * nxi:        Number of grid points in xi.
 * r:          Radial location (vector of size nr).
 * p:          Momentum of particle (vector of size np).
 * xi:         Pitch of particle (vector of size nxi).
 * f:          Distribution function. Linear array of size
 *             nr*nxi*np. The order of elements in the array
 *             should be p-first, xi-second, r-third
 *               f(p0,xi0,r0)  f(p1,xi0,r0)  f(p2,xi0,r0) ... f(pN,xi0,r0)  f(p0,xi1,r0) f(p1,xi1,r0) ... f(pN,xiN,r0) f(p0,xi0,r1) ... f(pN,xiN,rN)
 *             where N = np-1.
 * interptype: interpolation method to use
 *    (NumericDistributionFunction::INTERPOLATION_???).
 */
void NumericDistributionFunction::Initialize(
    const unsigned int nr, const unsigned int np, const unsigned int nxi,
    slibreal_t *r, slibreal_t *p, slibreal_t *xi,
    slibreal_t *f, int interptype
) {
    bool ascnd = true;
    unsigned int i;
    this->msdf = new NumericMomentumSpaceDistributionFunction[nr];
    this->r = r;
    this->p = p;
    this->xi = xi;
    this->f = f;
    this->nr = nr;
    this->np = np;
    this->nxi = nxi;
    this->interptype = interptype;

    for (i = 0; i < nr; i++) {
        if (interptype == INTERPOLATION_LINEAR)
            this->msdf[i].Initialize(np, nxi, p, xi, f+i*(np*nxi), NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR);
        else
            this->msdf[i].Initialize(np, nxi, p, xi, f+i*(np*nxi), NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC);
    }

    // Ensure that the radial grid is monotonically increasing or decreasing
    if (nr > 2)
        ascnd = (r[1] >= r[0]);

    for (i = 1; i < nr; i++) {
        if ((ascnd && r[i-1] >= r[i]) || (!ascnd && r[i-1] <= r[i]))
            throw SOFTLibException("Numerical distribution function: The radial grid must be monotonically increasing or decreasing.");
    }

    this->rmin = r[0];
    this->rmax = r[nr-1];
    if (this->rmin > this->rmax) {
        slibreal_t t = this->rmin;
        this->rmin = this->rmax;
        this->rmax = t;
    }
}

/**
 * Same as 'Initialize()', but first logarithmizes
 * f. Unless explicitly requested, f is substituted
 * with its corresponding logarithm.
 *
 * nr:         Number of grid points in r.
 * np:         Number of grid points in p.
 * nxi:        Number of grid points in xi.
 * r:          Radial location (vector of size nr).
 * p:          Momentum of particle (vector of size np).
 * xi:         Pitch of particle (vector of size nxi).
 * f:          Distribution function. Linear array of size
 *             nr*nxi*np. The order of elements in the array
 *             should be
 *               f(p0,xi0,r0)  f(p1,xi0,r0)  f(p2,xi0,r0) ... f(pN,xi0,r0)  f(p0,xi1,r0) f(p1,xi1,r0) ... f(pN,xiN,r0) f(p0,xi0,r1) ... f(pN,xiN,rN)
 *             where N = np-1.
 * interptype: interpolation method to use
 *               (NumericDistributionFunction::INTERPOLATION_???).
 * alloc:      If true, allocates a new array to store
 *             the logarithm of f in. Otherwise, the matrix
 *             f given here is overwritten (default).
 */
void NumericDistributionFunction::InitializeLog(
    const unsigned int nr, const unsigned int np, const unsigned int nxi,
    slibreal_t *r, slibreal_t *p, slibreal_t *xi,
    slibreal_t *f, int interptype, bool alloc
) {
	slibreal_t *tf;

    if (alloc) {
        tf = new slibreal_t[nr*np*nxi];
        for (unsigned int i = 1; i < nr; i++)
            tf[i] = tf[i-1] + np*nxi;
    } else tf = f;

    for (unsigned int i = 0; i < nr*np*nxi; i++) {
		if (f[i] <= 0)
			tf[i] = -std::numeric_limits<slibreal_t>::infinity();
		else
			tf[i] = log(f[i]);
    }

    this->logarithmic = true;
    this->Initialize(nr, np, nxi, r, p, xi, tf, interptype);

    for (unsigned int i = 0; i < nr; i++)
        this->msdf[i].SetLogarithmic(true);
}

/**
 * Clone this distribution function while
 * keeping a minimal memory footprint.
 */
NumericDistributionFunction *NumericDistributionFunction::MinClone() {
    NumericDistributionFunction *df = new NumericDistributionFunction();

    // No need to check if object is logarithmic, as the
    // stored data already is logarithmized in that case.
    df->Initialize(this->nr, this->np, this->nxi, this->r, this->p, this->xi, this->f, this->interptype);
    return df;
}

