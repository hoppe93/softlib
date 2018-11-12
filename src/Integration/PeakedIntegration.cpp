/**
 * Implementation of the 'peaked integration scheme'.
 */

#include <cmath>
#include <functional>
#include <softlib/config.h>
#include <softlib/Integration/PeakedIntegration.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructor.
 */
PeakedIntegration::PeakedIntegration() {}

/**
 * Evaluate the integral.
 *
 * f:       Function to integrate.
 * a, b:    Integration limits.
 * falloff: Characteristic falloff of f.
 * cutoff:  Fraction of maximum at which to cut off
 *          the remaning integration intervals.
 * maximum: (optional) Location of of the maximum
 *          of f. If omitted, the maximum is
 *          located automatically.
 */
slibreal_t PeakedIntegration::Evaluate(
    function<slibreal_t(slibreal_t)> &f,
    const slibreal_t a, const slibreal_t b,
    const slibreal_t falloff, const slibreal_t cutoff
) {
    slibreal_t m = FindMaximum(f);
    slibreal_t fmax = f(m);

    return Evaluate(f, a, b, falloff, cutoff, m, fmax);
}
slibreal_t PeakedIntegration::Evaluate(
    function<slibreal_t(slibreal_t)> &f,
    const slibreal_t a, const slibreal_t b,
    const slibreal_t falloff, const slibreal_t cutoff,
    const slibreal_t maximum, const slibreal_t fmax
) {
    slibreal_t s = fmax;
    slibreal_t
        cf = s*cutoff,
        fpval, fmval, delta;

    if (falloff == 0)
        throw SOFTLibException("PeakedIntegral: The fall-off must be non-zero.");
    
    delta = falloff;
    fpval = f(maximum+delta);
    fmval = f(maximum-delta);

    while (max(fpval, fmval) > cf) {
        if (fpval > fmax)
            return Evaluate(f, a, b, falloff, cutoff, maximum+delta, fpval);
        else if (fmval > fmax)
            return Evaluate(f, a, b, falloff, cutoff, maximum-delta, fmval);

        s += fpval + fmval;
        delta += falloff;

        if (maximum+delta < b)
            fpval = f(maximum+delta);
        else
            fpval = 0;

        if (maximum-delta > a)
            fmval = f(maximum-delta);
        else
            fmval = 0;
    }

    s += 0.5*(fpval + fmval);

    return (falloff*s);
}

/**
 * Locate the maximum of the given
 * unimodal function.
 *
 * f:   Unimodal function to find maximum of.
 * tol: Desired tolerance of the maximum.
 *
 * RETURNS the location of the maximum of f.
 */
slibreal_t PeakedIntegration::FindMaximum(
    function<slibreal_t(slibreal_t)> &f,
    const slibreal_t tol
) {
    const slibreal_t invphi = 0.5*(sqrt(5.0)-1.0), invphi2 = invphi*invphi;
    unsigned int i, nmax;
    slibreal_t a, b, c, d, h, yc, yd;

    // To leading order, we expect
    // a maximum at pi/2
    a = 0.0;
    b = 2.0*M_PI;
    h = b-a;

    c  = a + invphi2 * h;
    d  = a + invphi  * h;
    yc = f(c);
    yd = f(d);

    nmax = (unsigned int)(log(tol/h) / log(invphi));

    for (i = 0; i < nmax; i++) {
        if (yc > yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi*h;
            c = a + invphi2*h;
            yc = f(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi*h;
            d = a + invphi*h;
            yd = f(d);
        }
    }

    if (yc < yd)
        return 0.5*(a+d);
    else
        return 0.5*(b+c);
}

