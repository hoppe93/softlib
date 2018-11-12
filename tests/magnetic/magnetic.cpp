
#include <cmath>
#include <iostream>

#include <softlib/config.h>
#include <softlib/SOFTLibException.h>

#include "magfield_points.h"
#include "magnetic.h"
#include "runtest.h"

/**
 * Constructor
 */
Test_MagneticField::Test_MagneticField(const string& name) : UnitTest(name) {
	this->threshold = MAGNETIC_FIELD_TEST_GUARANTEED_PRECISION;
}
Test_MagneticField::Test_MagneticField(const string& name, const slibreal_t threshold) : UnitTest(name) {
	this->threshold = threshold;
}
/**
 * Run through the test points, evaluate them with softlib
 * and compare them to the expected values.
 *
 * points: The list of test points to run and compare to.
 * mf:     The magnetic field object to test.
 * B0:     The magnetic field strength on-axis (for normalization of result).
 * qprof:  Name of the safety factor used.
 */
bool Test_MagneticField::ComparePoints(const slibreal_t points[MAGNETIC_FIELD_TEST_NPOINTS][12], MagneticField2D *mf, slibreal_t B0, const string &qprof) {
	unsigned int i;
	slibreal_t x[3], B[3], *mB, d[3];
	for (i = 0; i < MAGNETIC_FIELD_TEST_NPOINTS; i++) {
		x[0] = points[i][0];
		x[1] = points[i][1];
		x[2] = points[i][2];

		B[0] = points[i][3];
		B[1] = points[i][4];
		B[2] = points[i][5];

		mB = mf->Eval(x);
		d[0] = fabs(B[0]-mB[0])/B0;
		d[1] = fabs(B[1]-mB[1])/B0;
		d[2] = fabs(B[2]-mB[2])/B0;

		if (d[0] > threshold || d[1] > threshold || d[2] > threshold) {
			throw SOFTLibException(
				"Computed magnetic field with '%s' q is above threshold. Max Delta = %e.",
                qprof.c_str(), max(d[0], max(d[1], d[2]))
			);
		}
	}

	return true;
}

/**
 * Run through the test points, evaluate the gradients
 * and curls of the magnetic field in those points, and
 * compare to the expected values.
 *
 * points: The list of test points to run and compare to.
 * mf:     The magnetic field object to test.
 * B0:     The magnetic field strength on-axis (for normalization of result).
 * qprof:  Name of the safety factor used.
 */
bool Test_MagneticField::CompareDerivatives(const slibreal_t points[MAGNETIC_FIELD_TEST_NPOINTS][12], MagneticField2D *mf, slibreal_t B0, const string &qprof) {
	unsigned int i;
	slibreal_t x[3],
		B[3], gradB[3], curlB[3],
		dB[3], dGrad[3], dCurl[3];
	struct magnetic_field_data mfd;

	for (i = 0; i < MAGNETIC_FIELD_TEST_NPOINTS; i++) {
		x[0] = points[i][0];
		x[1] = points[i][1];
		x[2] = points[i][2];

		B[0] = points[i][3];
		B[1] = points[i][4];
		B[2] = points[i][5];

		gradB[0] = points[i][6];
		gradB[1] = points[i][7];
		gradB[2] = points[i][8];

		curlB[0] = points[i][9];
		curlB[1] = points[i][10];
		curlB[2] = points[i][11];

		mfd = mf->EvalDerivatives(x);

		dB[0]    = fabs(mfd.B[0]-B[0])/B0;
		dB[1]    = fabs(mfd.B[1]-B[1])/B0;
		dB[2]    = fabs(mfd.B[2]-B[2])/B0;

		dGrad[0] = fabs(mfd.gradB[0]-gradB[0])/B0;
		dGrad[1] = fabs(mfd.gradB[1]-gradB[1])/B0;
		dGrad[2] = fabs(mfd.gradB[2]-gradB[2])/B0;

		dCurl[0] = fabs(mfd.curlB[0]-curlB[0])/B0;
		dCurl[1] = fabs(mfd.curlB[1]-curlB[1])/B0;
		dCurl[2] = fabs(mfd.curlB[2]-curlB[2])/B0;

		if (dB[0] > threshold || dB[1] > threshold || dB[2] > threshold) {
			char buffer[32];
			snprintf(buffer, 31, "%e", max(dB[0], max(dB[1], dB[2])));

			throw SOFTLibException(
				"Computed magnetic field with "+qprof+" q is above threshold. "+
				"max d = "+buffer
			);
		}
		if (dGrad[0] > threshold || dGrad[1] > threshold || dGrad[2] > threshold) {
			char buffer[32];
			snprintf(buffer, 31, "%e", max(dGrad[0], max(dGrad[1], dGrad[2])));

			throw SOFTLibException(
				"Computed gradient of magnetic field with "+qprof+" q is above threshold. "+
				"max d = "+buffer
			);
		}
		if (dCurl[0] > threshold || dCurl[1] > threshold || dCurl[2] > threshold) {
			char buffer[32];
			snprintf(buffer, 31, "%e", max(dCurl[0], max(dCurl[1], dCurl[2])));

			throw SOFTLibException(
				"Computed curl of magnetic field with "+qprof+" q is above threshold. "+
				"max d = "+buffer
			);
		}
	}

	return true;
}

/**
 * Set the threshold value, above which
 * the compared values should be considered
 * as deviating.
 *
 * thr: Threshold value.
 */
void Test_MagneticField::SetThreshold(const slibreal_t thr) {
	this->threshold = thr;
}
