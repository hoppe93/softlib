/**
 * Implements general methods for testing the distribution function
 * classes of SOFTLib.
 */

#include <cmath>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include "Test_DistributionFunction.h"

/**
 * Evaluates the given distribution function
 * in the given points and compares the result
 * to the given table of values.
 *
 * df:   SOFTLib distribution function to test.
 * r:    1-D list of radial points at which f is given.
 * p:    1-D list of momentum points at which f is given.
 * xi:   1-D list of pitch points at which f is given.
 * f:    Expected value of the distribution function
 *       in the given points.
 * ntab: Number of points in table, i.e. number of
 *       elements in each of r, p, xi and f.
 * tol:  The maximum allowed relative error for the test
 *       to be considered a success.
 *
 * RETURNS the number of points that had a relative error
 * larger than 'tol'.
 */
unsigned int Test_DistributionFunction::CompareToTable(
    DistributionFunction *df, const slibreal_t *r,
    const slibreal_t *p, const slibreal_t *xi,
    const slibreal_t *f, const unsigned int ntab,
    const slibreal_t tol
) {
    slibreal_t tf, err;
    unsigned int i, failed=0;
    for (i = 0; i < ntab; i++) {
        tf = df->Eval(r[i], p[i], xi[i]);

        if (f[i] != 0.0) err = fabs((tf-f[i])/f[i]);
        else err = fabs(tf-f[i]);

        if (isnan(err) || err > tol) {
            //printf("%10e  %10e  %+12e  %+12e  %4.2f%%\n", p[i], xi[i], f[i], tf, err*100.0);
            failed++;
        }
    }

    return failed;
}
