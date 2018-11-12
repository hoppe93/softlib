/**
 * Implementation of unit tests for CODE distribution function
 * in SOFTLib.
 */

#include <cmath>
#include <ctime>

#include <softlib/config.h>
#include <softlib/DistributionFunction/CODEDistributionFunction.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>

#include "Test_CODEDistributionFunction.h"

/**
 * Create a numeric distribution function from the
 * given CODE distribution and compare the two in
 * random points.
 */
bool Test_CODEDistributionFunction::CreateNumericAndTest(CODEDistributionFunction *cdf) {
    unsigned int failed, i;
    slibreal_t *r, *p, *xi, *f, pmin, pmax;
    NumericMomentumSpaceDistributionFunction *msdf
        = cdf->ToMomentumSpace(false, false, 200, NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR);

    pmin = cdf->GetPMin();
    pmax = cdf->GetPMax();

    srand(time(NULL));

    // Generate random points
    r  = new slibreal_t[NPOINTS];
    p  = new slibreal_t[NPOINTS];
    xi = new slibreal_t[NPOINTS];
    f  = new slibreal_t[NPOINTS];
    for (i = 0; i < NPOINTS; i++) {
        r[i]  = 0.0;
        p[i]  = pmin + (pmax-pmin) * ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
        // Numeric distributions are often under-resolved at xi < 0,
        // so we ignore those points...
        xi[i] = ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);

        f[i] = cdf->Eval(p[i], xi[i]);
    }

    failed = this->CompareToTable(msdf, r, p, xi, f, NPOINTS, TOLERANCE2);

    delete [] r;
    delete [] p;
    delete [] xi;
    delete [] f;

    if (failed > NPOINTS*ALLOWED_FAILURES) {
        this->PrintError(
            "Re-gridded function: Too many failures. %u of %u (%.2f%%) test points had an error larger than %.2f%%.",
            failed, NPOINTS, ((double)failed)/((double)NPOINTS)*100.0, TOLERANCE2*100.0
        );
        return false;
    } else if (failed != 0) {
        this->PrintWarning(
            "Re-gridded function: %u of %u (%.2f%%) test points had an error larger than %.2f%%.",
            failed, NPOINTS, ((double)failed)/((double)NPOINTS)*100.0, TOLERANCE2*100.0
        );
    }

    this->PrintOK("Conversion of CODE distribution to numeric momentum-space distribution function passed.");
    return true;
}

/**
 * Load the 'correct values' values file and compare
 * to the values of the given CODE distribution function.
 *
 * fname: Name of file containing correct values.
 * cdf:   CODE distribution function object.
 */
bool Test_CODEDistributionFunction::LoadAndTestCorrect(
    const string& fname, CODEDistributionFunction *cdf
) {
    double **tp, **txi, **tf, *tr;
    unsigned int failed;
    sfilesize_t fsize[2], npoints, i;

    SFile *sf = SFile::Create(fname, SFILE_MODE_READ);

    tp = sf->GetDoubles("p", fsize);
    npoints = fsize[1];

    txi = sf->GetDoubles("xi", fsize);
    if (fsize[1] != npoints)
        throw SOFTLibException("Invalid number of points in 'correct data' xi variable. Expected %llu, got %llu.", npoints, fsize[0]);

    tf = sf->GetDoubles("f", fsize);
    if (fsize[1] != npoints)
        throw SOFTLibException("Invalid number of points in 'correct data' f variable. Expected %llu, got %llu.", npoints, fsize[0]);

    sf->Close();

    tr = new slibreal_t[npoints];
    for (i = 0; i < npoints; i++)
        tr[i] = 0.0;

    failed = this->CompareToTable(cdf, tr, tp[0], txi[0], tf[0], (unsigned int)npoints, TOLERANCE1);
    
    delete [] tr;

    if (failed > npoints*ALLOWED_FAILURES) {
        this->PrintError(
            "Comparison to correct: Too many failures. %u of %u (%.2f%%) test points had an error larger than %.2f%%.",
            failed, npoints, ((double)failed)/((double)npoints)*100.0, TOLERANCE1*100.0
        );
        return false;
    } else if (failed != 0) {
        this->PrintWarning(
            "Comparison to correct: %u of %u (%.2f%%) test points had an error larger than %.2f%%.",
            failed, npoints, ((double)failed)/((double)npoints)*100.0, TOLERANCE1*100.0
        );
    }

    this->PrintOK("Direct evaluation of CODE distribution function passed.");
    return true;
}

/**
 * Run the tests.
 */
bool Test_CODEDistributionFunction::Run() {
    bool success = true;
    CODEDistributionFunction *cdf
        = new CODEDistributionFunction(this->distfilename);

    // Test 1 (compare to correct values)
    success = LoadAndTestCorrect(this->corrfilename, cdf);

    // Test 2 (create numeric momentum-space distribution and compare)
    success = CreateNumericAndTest(cdf);

    cdf->Destroy();

    return success;
}

