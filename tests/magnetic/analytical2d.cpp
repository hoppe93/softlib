
#include <cmath>
#include <iostream>
#include <limits>
#include <runtest.h>

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>
#include <softlib/SOFTLibException.h>

#include "analytical2d.h"
#include "magfield_points.h"

using namespace std;

#define THRESHOLD MAGNETIC_FIELD_TEST_GUARANTEED_PRECISION

/**
 * Constructor
 */
Test_MagneticFieldAnalytical2D::Test_MagneticFieldAnalytical2D(const string& name) : Test_MagneticField(name) {}

/**
 * Run the test.
 * RETURNS true on success, false on failure.
 * On failure, an error message is also printed.
 */
bool Test_MagneticFieldAnalytical2D::Run() {
    bool success = true;
	slibreal_t
		B0       = magnetic_field_test_data_B0,
		Rm       = magnetic_field_test_data_Rm,
		rminor   = magnetic_field_test_data_rminor,
		qa_const = magnetic_field_test_data_qa_const,
		qa_lin   = magnetic_field_test_data_qa_lin,
		qa_quad  = magnetic_field_test_data_qa_quad,
		qa_exp   = magnetic_field_test_data_qa_exp;

	MagneticFieldAnalytical2D *mf_const, *mf_lin, *mf_quad, *mf_exp;

	mf_const = new MagneticFieldAnalytical2D(B0, Rm, rminor, MFAFS_CW, MFAFS_CCW, MFASF_CONSTANT, qa_const, 0.0);
	mf_lin   = new MagneticFieldAnalytical2D(B0, Rm, rminor, MFAFS_CW, MFAFS_CCW, MFASF_LINEAR, qa_lin, 1.0);
	mf_quad  = new MagneticFieldAnalytical2D(B0, Rm, rminor, MFAFS_CW, MFAFS_CCW, MFASF_QUADRATIC, qa_quad, 1.0);
	mf_exp   = new MagneticFieldAnalytical2D(B0, Rm, rminor, MFAFS_CW, MFAFS_CCW, MFASF_EXPONENTIAL, qa_exp, 0.0);

	// Test list of values
	try {
		success &= ComparePoints(magnetic_field_test_data_const, mf_const, B0, "constant", true);
		success &= ComparePoints(magnetic_field_test_data_linear, mf_lin, B0, "linear", true);
		success &= ComparePoints(magnetic_field_test_data_quadratic, mf_quad, B0, "quadratic", true);
		success &= ComparePoints(magnetic_field_test_data_exponential, mf_exp, B0, "exponential", true);
	} catch (SOFTLibException &ex) {
		this->PrintError("[MagneticFieldAnalytical2D]: "+ex.whats());
		return false;
	}

    if (success)
        this->PrintOK("Magnetic field evaluation");
    else
        return false;

	// Test derivatives
	try {
		success &= CompareDerivatives(magnetic_field_test_data_const, mf_const, B0, "constant");
		success &= CompareDerivatives(magnetic_field_test_data_linear, mf_lin, B0, "linear");
		success &= CompareDerivatives(magnetic_field_test_data_quadratic, mf_quad, B0, "quadratic");
		success &= CompareDerivatives(magnetic_field_test_data_exponential, mf_exp, B0, "exponential");
	} catch (SOFTLibException &ex) {
		this->PrintError("[MagneticFieldAnalytical2D]: "+ex.whats());
		return false;
	}

    if (success)
        this->PrintOK("Magnetic field derivative evaluation");
    else
        return false;

    // Test conversion to numeric magnetic field
    try {
        success &= TestConversion(mf_const, "constant");
        success &= TestConversion(mf_lin, "linear");
        success &= TestConversion(mf_quad, "quadratic");
        success &= TestConversion(mf_exp, "exponential");
    } catch (SOFTLibException& ex) {
        this->PrintError("[MagneticFieldAnalytical2D]: "+ex.whats());
        return false;
    }

    if (success)
        this->PrintOK("Magnetic field analytical to numeric conversion.");
    else
        return false;

    // Test Jacobian 1
    success &= TestJacobian1(mf_const);
    success &= TestJacobian1(mf_lin);
    success &= TestJacobian1(mf_quad);
    success &= TestJacobian1(mf_exp);

    if (success)
        this->PrintOK("Magnetic field Jacobian evaluation #1.");
    else
        return false;

    // Test Jacobian 2 (more thorough)
    success &= TestJacobian2(mf_const);
    success &= TestJacobian2(mf_lin);
    success &= TestJacobian2(mf_quad);
    success &= TestJacobian2(mf_exp);

    if (success)
        this->PrintOK("Magnetic field Jacobian evaluation #2.");
    else
        return false;

    delete mf_const;
    delete mf_lin;
    delete mf_quad;
    delete mf_exp;

	return success;
}

/**
 * Check if the analytical magnetic field can be converted into a
 * numeric magnetic field without errors.
 *
 * mf:   Analytical magnetic field to convert.
 * name: Name of magnetic field.
 */
bool Test_MagneticFieldAnalytical2D::TestConversion(MagneticFieldAnalytical2D *mf, const string& name) {
    unsigned int i;
    slibreal_t Delta, Rmin, Rmax, Zmin, Zmax, r, p, z;

    Rmax = mf->GetRMajor() - mf->GetRMinor();
    Rmin = mf->GetRMajor() + mf->GetRMinor();
    Zmax = mf->GetRMinor();
    Zmin =-mf->GetRMinor();

    MagneticFieldNumeric2D *mfn = mf->ToNumeric2D();
    for (i = 0; i < NRANDOM; i++) {
        r = Rand(Rmin, Rmax);
        p = Rand(0.0, 2.0*M_PI);
        z = Rand(Zmin, Zmax);

        Vector<3> Ba =  mf->Eval(r*cos(p), r*sin(p), z);
        Vector<3> Bn = mfn->Eval(r*cos(p), r*sin(p), z);

        Delta = fabs((Ba[0]-Bn[0])/Ba[0]);
        if (Delta > CONVTOL) {
            this->PrintError("%s: Analytical to numerical-converted magnetic field x-component does not evaluate to the correct thing. Delta = %e.", name.c_str(), Delta);
            delete mfn;
            return false;
        }

        Delta = fabs((Ba[1]-Bn[1])/Ba[1]);
        if (Delta > CONVTOL) {
            this->PrintError("%s: Analytical to numerical-converted magnetic field y-component does not evaluate to the correct thing. Delta = %e.", name.c_str(), Delta);
            delete mfn;
            return false;
        }

        Delta = fabs((Ba[2]-Bn[2])/Ba[2]);
        if (Delta > CONVTOL) {
            this->PrintError("%s: Analytical to numerical-converted magnetic field z-component does not evaluate to the correct thing. Delta = %e.", name.c_str(), Delta);
            delete mfn;
            return false;
        }
    }

    return true;
}

/**
 * Test whether the magnetic field Jacobian is being evaluated
 * correctly by comparing it to the divergence and curl of B.
 * By Gauss' law, the divergence must be zero, and since the
 * curl of B is evaluated and tested separately in the analytical
 * field, it can be compared to.
 *
 * mf: Magnetic field object to test.
 */
bool Test_MagneticFieldAnalytical2D::TestJacobian1(MagneticFieldAnalytical2D *mf) {
    unsigned int i;
    slibreal_t r, p, z, Rmin, Rmax, Zmin, Zmax, Delta;

    Rmax = mf->GetRMajor() - mf->GetRMinor();
    Rmin = mf->GetRMajor() + mf->GetRMinor();
    Zmax = mf->GetRMinor();
    Zmin =-mf->GetRMinor();

    for (i = 0; i < NRANDOM; i++) {
        r = Rand(Rmin, Rmax);
        p = Rand(0.0, 2.0*M_PI);
        z = Rand(Zmin, Zmax);

        struct magnetic_field_data mfd = mf->EvalDerivatives(r*cos(p), r*sin(p), z);
        
        // Divergence of B (should vanish)
        Delta = fabs(mfd.J[0][0] + mfd.J[1][1] + mfd.J[2][2]);
        if (Delta > JAC1TOL) {
            this->PrintError("%u: Jacobian trace does not vanish. Delta = %e.", i, Delta);
            return false;
        }

        // (curl B)_x
        if (mfd.curlB[0] == 0.0)
            Delta = fabs((mfd.J[2][1]-mfd.J[1][2]) - mfd.curlB[0]);
        else
            Delta = fabs(((mfd.J[2][1]-mfd.J[1][2]) - mfd.curlB[0]) / mfd.curlB[0]);

        if (Delta > JAC1TOL) {
            this->PrintError("%u: x-component of curl B calculated from Jacobian is wrong. Delta = %e.", i, Delta);
            return false;
        }

        // (curl B)_y
        if (mfd.curlB[1] == 0.0)
            Delta = fabs((mfd.J[0][2]-mfd.J[2][0]) - mfd.curlB[1]);
        else
            Delta = fabs(((mfd.J[0][2]-mfd.J[2][0]) - mfd.curlB[1]) / mfd.curlB[1]);

        if (Delta > JAC1TOL) {
            this->PrintError("%u: y-component of curl B calculated from Jacobian is wrong. Delta = %e.", i, Delta);
            return false;
        }

        // (curl B)_z
        if (mfd.curlB[2] == 0.0)
            Delta = fabs((mfd.J[1][0]-mfd.J[0][1]) - mfd.curlB[2]);
        else
            Delta = fabs(((mfd.J[1][0]-mfd.J[0][1]) - mfd.curlB[2]) / mfd.curlB[2]);

        if (Delta > JAC1TOL) {
            this->PrintError("%u: z-component of curl B calculated from Jacobian is wrong. Delta = %e.", i, Delta);
            return false;
        }
    }

    return true;
}

/**
 * Test whether the magnetic field Jacobian is being evaluated
 * correctly by generating a numerical magnetic field from which
 * the jacobian can be evaluated numerically. The numeric magnetic
 * field only depends on the correct evaluation of B, not its
 * derivatives, and so this test can detect algebraic errors in
 * the analytical magnetic field.
 *
 * mf: Magnetic field object to test.
 */
bool Test_MagneticFieldAnalytical2D::TestJacobian2(MagneticFieldAnalytical2D *mf) {
    slibreal_t r, p, z, Delta,
        Rmin, Rmax, Zmin, Zmax;
    MagneticFieldNumeric2D *mfn;
    unsigned int i, j, k;
    
    mfn = mf->ToNumeric2D(NR, NZ);

    Rmax = mf->GetRMajor() + mf->GetRMinor();
    Rmin = mf->GetRMajor() - mf->GetRMinor();
    Zmax = mf->GetRMinor();
    Zmin =-mf->GetRMinor();

    /* Do the testing */
    struct magnetic_field_data a, n;
    for (i = 0; i < NRANDOM; i++) {
        r = Rand(Rmin, Rmax);
        p = Rand(0.0, 2.0*M_PI);
        z = Rand(Zmin, Zmax);

        a =  mf->EvalDerivatives(r*cos(p), r*sin(p), z);
        n = mfn->EvalDerivatives(r*cos(p), r*sin(p), z);

        Delta = fabs((a.Babs-n.Babs)/a.Babs);
        if (Delta > CONVTOL) {
            this->PrintError(
                "Magnetic field absolute value did not match numeric equivalent. Delta = %e.",
                Delta
            );

            delete mfn;
            return false;
        }

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                if (a.J[j][k] == 0.0)
                    Delta = fabs(a.J[j][k] - n.J[j][k]);
                else
                    Delta = fabs((a.J[j][k] - n.J[j][k]) / n.J[j][k]);

                if (Delta > JAC2TOL) {
                    this->PrintError(
                        "%u: Jacobian (%u, %u) component did not match numeric equivalent. Delta = %e.",
                        i, j, k, Delta
                    );

                    delete mfn;
                    return false;
                }
            }
        }
    }

    delete mfn;

    return true;
}

