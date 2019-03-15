#ifndef _TEST_NUMERIC_DISTRIBUTION_FUNCTION_H
#define _TEST_NUMERIC_DISTRIBUTION_FUNCTION_H

#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/NumericMomentumSpaceDistributionFunction.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Test_DistributionFunction.h"
#include "runtest.h"

class Test_NumericDistributionFunction : public Test_DistributionFunction {
    private:
        const unsigned int
            NPOINTS1        = 300,
            NPOINTS2        = 300;
        const slibreal_t
            EHAT            = 4.0,
            LOGLAMBDA       = 17.0,
            ZEFF            = 2.0,

            RMIN            = 1.6,
            RMAX            = 2.2,
            PMIN            = 0.1,
            PMAX            = 100.0,
            TOLERANCE1      = 0.2,
            TOLERANCE2      = 0.05,
            ALLOWED_FAILURES= 0.1;

    public:
        Test_NumericDistributionFunction(const string& name) : Test_DistributionFunction(name) {}
        RadialDistributionFunction *GeneratePhaseSpaceDistribution(
            slibreal_t, slibreal_t, slibreal_t,
            slibreal_t, slibreal_t
        );
        unsigned int RandomPoints(
            DistributionFunction*, DistributionFunction*,
            const unsigned int, const slibreal_t,
            const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t
        );
		template<class T> bool RunInternal(const string&);
};

/******************************************
 * Implementation of templated functions. *
 ******************************************/
/**
 * Run the tests.
 */
template<class T>
bool Test_NumericDistributionFunction::RunInternal(const string& inputfilename) {
    unsigned int failed;
    slibreal_t perc;
    bool success = true;

	// Create dummy magnetic field (we just need the
	// radial interval)
	MagneticFieldAnalytical2D *mfa =
		new MagneticFieldAnalytical2D(5.0, 1.6, 0.6, MFATFS_CW, MFASF_CONSTANT, 1.0, 1.0);

    // Test full phase-space distribution
    RadialDistributionFunction *adf1 = GeneratePhaseSpaceDistribution(
		EHAT, LOGLAMBDA, ZEFF, RMIN, RMAX
	);
    T *df1
        = new T(inputfilename, mfa);

	// Generate points on a slightly larger grid to
	// test extrapolation
    failed = RandomPoints(
        adf1, df1, NPOINTS1, TOLERANCE1, RMIN-0.05, RMAX-0.05, PMIN, PMAX
    );

    if (failed > 0) {
        perc = ((slibreal_t)failed) / ((slibreal_t)NPOINTS1) * 100.0;
        this->PrintWarning(
            "--> Numeric distribution: %u of %u (%.2f%%) tests had an error larger than %.2f%%.",
            failed, NPOINTS1, perc, TOLERANCE1*100.0
		);
    }
    
    if (failed > NPOINTS1 * ALLOWED_FAILURES) {
        this->PrintError("--> Too many failures for numeric distribution function.");
        success = false;
    }

    // Test conversion to momentum-space distribution
    AnalyticalAvalanche *adf2 = new AnalyticalAvalanche();
    adf2->Initialize(EHAT, LOGLAMBDA, ZEFF);

    NumericMomentumSpaceDistributionFunction *df2
        = T::LoadMomentumSpace(inputfilename, mfa);

    failed = RandomPoints(
        adf2, df2, NPOINTS2, TOLERANCE2, 0.0, 0.0, PMIN, PMAX
    );

    if (failed > 0) {
        perc = ((slibreal_t)failed) / ((slibreal_t)NPOINTS2) * 100.0;
        this->PrintWarning(
            "--> Numeric distribution: %u of %u (%.2f%%) tests had an error larger than %.2f%%.",
            failed, NPOINTS2, perc, TOLERANCE2*100.0);
    }
    
    if (failed > NPOINTS2 * ALLOWED_FAILURES) {
        this->PrintError("--> Too many failures for numeric distribution function.");
        success = false;
    }

    return success;
}

#endif/*_TEST_NUMERIC_DISTRIBUTION_FUNCTION_H*/
