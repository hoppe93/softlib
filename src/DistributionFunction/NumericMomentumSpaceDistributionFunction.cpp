/**
 * Implementation of a distribution function.
 */

#include <cmath>
#include <cstdio>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <softlib/DistributionFunction/MomentumSpaceDistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Evaluate this distribution function in the
 * given phase-space point.
 *
 * r:   Radial location to evaluate f at.
 * p:   Magnitude of momentum to evaluate f at.
 * xi:  Cosing of pitch angle to evaluate f at.
 * drift_shift: Shift of drift surface (radial coordinate)
 *    due to orbit drifts.
 */
slibreal_t NumericMomentumSpaceDistributionFunction::Eval(const slibreal_t p, const slibreal_t xi) {
    slibreal_t fval;
    if (p < this->pmin || p > this->pmax) return 0.0;
    if (xi < this->ximin || xi > this->ximax) return 0.0;

    if (this->flipPitchSign)
        fval = gsl_spline2d_eval(fspline, p, -xi, pa, xia);
    else
        fval = gsl_spline2d_eval(fspline, p, xi, pa, xia);

    if (this->logarithmic) {
        if (isnan(fval))
            return 0.0;
        else
            return (slibreal_t)exp(fval);
    } else
        return (slibreal_t)fval;
}

/**
 * Initialize the momentum-space distribution
 * function.
 *
 * np:     Number of momentum points (resolution in p)
 * nxi:    Number of pitch points (resolution in xi)
 * p:      Momentum grid (p in units of mc)
 * xi:     Pitch grid (cosine of pitch angle)
 * f:      Momentum-space distribution function
 * interp: Interpolation method to use
 *         (NumericMomentumSpaceDistributionFunction::INTERPOLATION_???)
 */
void NumericMomentumSpaceDistributionFunction::Initialize(
    const unsigned int np, const unsigned int nxi,
    slibreal_t *p, slibreal_t *xi, slibreal_t *f,
    int interp
) {
    this->np = np;
    this->nxi = nxi;
    this->p = p;
    this->xi = xi;
    this->f = f;

    this->pmin = p[0];
    this->pmax = p[np-1];
    this->ximin = xi[0];
    this->ximax = xi[nxi-1];

    if (this->pmin > this->pmax)
        throw SOFTLibException("The p-grid must be strictly increasing.");

    if (this->ximin > this->ximax)
        throw SOFTLibException("The xi-grid must be strictly increasing.");

    // Initialize interpolation
    pa = gsl_interp_accel_alloc();
    xia = gsl_interp_accel_alloc();

    if (interp == NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR) {
        fspline = gsl_spline2d_alloc(gsl_interp2d_bilinear, np, nxi);
        this->interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR;
    } else {
        fspline = gsl_spline2d_alloc(gsl_interp2d_bicubic, np, nxi);
        this->interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC;
    }

    gsl_spline2d_init(fspline, this->p, this->xi, this->f, np, nxi);
}

/**
 * Initialize the momentum-space distribution function
 * with the logarithm of the given f.
 *
 * np:     Number of momentum points (resolution in p)
 * nxi:    Number of pitch points (resolution in xi)
 * p:      Momentum grid (p in units of mc)
 * xi:     Pitch grid (cosine of pitch angle)
 * f:      Momentum-space distribution function. This
 *    argument will be logarithmizead (base e)
 *    before interpolation.
 * interp: Interpolation method to use
 *         (NumericMomentumSpaceDistributionFunction::INTERPOLATION_???)
 */
void NumericMomentumSpaceDistributionFunction::InitializeLog(
    const unsigned int np, const unsigned int nxi,
    slibreal_t *p, slibreal_t *xi, slibreal_t *f,
    int interp
) {
    unsigned int i;
    slibreal_t *logf = new slibreal_t[np*nxi];

    for (i = 0; i < np*nxi; i++) {
        if (f[i] < 0.0)
            throw SOFTLibException("Distribution function cannot be negative.");
        logf[i] = log(f[i]);
    }

    this->logarithmic = true;
    this->Initialize(np, nxi, p, xi, logf, interp);
}

/**
 * Clone this distribution function while
 * keeping a minimal memory footprint.
 */
NumericMomentumSpaceDistributionFunction *NumericMomentumSpaceDistributionFunction::MinClone() {
    NumericMomentumSpaceDistributionFunction *df = new NumericMomentumSpaceDistributionFunction();

    df->Initialize(this->np, this->nxi, this->p, this->xi, this->f, this->interptype);
    df->FlipPitchSign(this->flipPitchSign);
    return df;
}

