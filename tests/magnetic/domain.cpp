/**
 * Test of the domain handling logic of 'MagneticField2D'.
 */

#include <cmath>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>
#include "domain.h"

/* Wall section drawn by Mathias */
const unsigned int NWALL=34;
slibreal_t test_domain[2][NWALL] = {
	{2.3066,2.3066,2.3953,2.3953,2.7294,2.7294,2.7844,3.1734,3.1734,2.9323,2.8774,3.1564,3.2791,3.296,3.389,3.4186,3.537,3.5666,3.4186,3.4059,3.5539,3.7104,3.685,3.6216,3.5877,3.4482,3.1311,3.1311,3.2833,3.2114,2.9323,2.9027,2.3066,2.3066},
    {1.4974,-0.16582,-0.27806,-0.81888,-0.8699,-1.0689,-1.1862,-1.1862,-1.0893,-1.074,-0.75255,-0.57908,-0.55357,-0.58929,-0.58418,-0.6301,-0.58929,-0.24745,-0.22194,-0.09949,0.084184,0.42602,0.59949,0.8852,1.0842,1.1811,1.2372,1.3699,1.375,1.574,1.6097,1.6862,1.6862,1.4974}
};
const unsigned int NINPOINTS=24;
slibreal_t test_domain_inside[2][NINPOINTS] = {
    {2.4503,2.4249,2.4799,2.6364,2.8224,2.8351,2.8647,2.7844,3.4355,3.4778,3.334,3.3044,3.592,2.4968,3.5412,3.055,3.0677,3.1818,3.186,2.9239,2.5899,2.7505,2.9831,2.7928},
    {1.4821,0.34439,-0.32908,-0.74745,-0.83929,-1.0995,-0.84439,-0.65561,-0.47194,-0.34439,-0.26786,0.017857,0.47194,0.47704,0.99745,1.2168,1.4413,1.4158,1.4872,1.5434,1.6046,0.875,0.25765,-0.21684}
};
const unsigned int NOUTPOINTS=18;
slibreal_t test_domain_outside[2][NOUTPOINTS] = {
    {2.2008,2.6321,2.5941,2.852,3.0846,3.055,3.63,3.4609,3.7273,3.4017,3.3044,3.0254,3.0169,2.5856,2.2431,2.5137,2.2093,2.9323},
    {-0.66071,-0.6301,-1.2168,-1.0842,-0.875,-0.42602,-0.55357,0.038265,0.70663,1.0842,1.2628,1.3342,1.7321,1.6097,1.574,0.72194,0.22194,0.45663}
};
const slibreal_t Test_Domain_magnetic_axis[2] = {2.8,0.0};    // Arbitrary numbers, inside domain

Test_Domain::Test_Domain(const string& name) : UnitTest(name) {}

MagneticField2D *Test_Domain::GenerateMF() {
	unsigned int i, np = 5;
	slibreal_t R[np], Z[np],
		Br[np*np], Bphi[np*np], Bz[np*np];

	for (i = 0; i < np; i++)
		R[i] = Z[i] = (slibreal_t)i;
	for (i = 0; i < np*np; i++)
		Br[i] = Bphi[i] = Bz[i] = 0.0;

	return new MagneticFieldNumeric2D(
		"domain", "domain", R, Z, np, np, Br, Bphi, Bz, nullptr,
        Test_Domain_magnetic_axis[0], Test_Domain_magnetic_axis[1],
		NULL, NULL, 0, test_domain[0], test_domain[1], NWALL
	);
}

/**
 * Run the test.
 * RETURNS true on success, false on failure.
 * On failure, an error message is also printed.
 */
bool Test_Domain::Run() {
	unsigned int i;
	slibreal_t x, y, s, c, x2, y2;

	MagneticField2D *mf = GenerateMF();
	srand(time(NULL));

	/* Test 'inside points' */
	for (i = 0; i < NINPOINTS-1; i++) {
		x = test_domain_inside[0][i];
		y = test_domain_inside[1][i];
		s = ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
		c = sqrt(1.0 - s*s);

		/* Test with next 'inside' point */
		x2 = test_domain_inside[0][i+1];
		y2 = test_domain_inside[1][i+1];

		/* Test with another inside point */
		if (mf->CrossesDomain(x*s, x*c, y, x2, 0.0, y2)) {
			this->PrintError("[Domain]: The two points (%u) and (%u), both located inside, are marked as lying outside domain.", i, i+1);
			return false;
		}
	}

	/* Test 'outside' points */
	for (i = 0; i < NOUTPOINTS; i++) {
		x = test_domain_outside[0][i];
		y = test_domain_outside[1][i];
		s = ((slibreal_t)rand()) / ((slibreal_t)RAND_MAX);
		c = sqrt(1.0 - s*s);

		/* Test with next 'outside' point */
		x2 = test_domain_outside[0][i+1];
		y2 = test_domain_outside[1][i+1];

		/* Test with another inside point */
		if (!mf->CrossesDomain(x*s, x*c, y, x2, 0.0, y2)) {
			this->PrintError("[Domain]: The two points (%u) and (%u), located inside and outside, are marked as lying inside domain.", i, i+1);
			return false;
		}
	}

	return true;
}

