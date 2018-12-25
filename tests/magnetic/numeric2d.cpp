/**
 * Test of the numeric 2D magnetic field.
 */ 
#include <iostream>
#include <runtest.h>

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/SFileException.h>
#include <softlib/SOFTLibException.h>
#include <cmath>
#include "numeric2d.h"
#include "magfield_points.h"

using namespace std;

#define NR 99
#define NZ 95
#define NWALL 100
#define TESTFILE "../tests/magnetic/circular-benchmark.mat"

#define THRESHOLD 5e-7
#define PI 3.14159265359

/**
 * Constructor
 */
Test_MagneticFieldNumeric2D::Test_MagneticFieldNumeric2D(const string& name) : Test_MagneticField(name, THRESHOLD) {}

MagneticFieldNumeric2D *Test_MagneticFieldNumeric2D::GenerateMF(
	MagneticFieldAnalytical2D *mfa, slibreal_t Rm, slibreal_t rminor,
	unsigned int nr, unsigned int nz, unsigned int nwall
) {
	unsigned int i, j;
	slibreal_t *Br, *Bphi, *Bz, *R, *Z, *B, *rwall, *zwall, x;

	/* Generate the grid */
	R = new slibreal_t[nr];
	Z = new slibreal_t[nz];

	for (i = 0; i < nr; i++)
		R[i] = Rm + rminor*((((slibreal_t)(2*i))/(nr-1)) - 1);
	for (i = 0; i < nz; i++)
		Z[i] = rminor*((((slibreal_t)(2*i))/(nz-1)) - 1);

	/* Generate the magnetic data */
	Br  = new slibreal_t[nr*nz];
	Bphi= new slibreal_t[nr*nz];
	Bz  = new slibreal_t[nr*nz];

	for (i = 0; i < nr; i++) {
		for (j = 0; j < nz; j++) {
			B = mfa->Eval(R[i], 0.0, Z[j]);

			/* When phi = 0,
			 * 	 Br = Bx, Bphi = By, Bz =Bz */
			/*Br[i*nz + j]   = B[0];
			Bphi[i*nz + j] = B[1];
			Bz[i*nz + j]   = B[2];*/
            Br[i + j*nr]   = B[0];
            Bphi[i + j*nr] = B[1];
            Bz[i + j*nr]   = B[2];
		}
	}
	
	/* Generate the wall */
	rwall = new slibreal_t[nwall];
	zwall = new slibreal_t[nwall];
	for (i = 0; i < nwall; i++) {
		x = 2.0*PI*(((slibreal_t)i)/(nwall-1));
		rwall[i] = Rm + rminor*cos(x);
		zwall[i] =      rminor*sin(x);
	}

	return new MagneticFieldNumeric2D(
		"numeric", "numeric", R, Z, nr, nz, Br, Bphi, Bz, nullptr,
        Rm, 0.0, NULL, NULL, 0, rwall, zwall, nwall
	);
}

bool Test_MagneticFieldNumeric2D::CheckNumericDerivatives(
    MagneticFieldAnalytical2D *mfa,
    MagneticFieldNumeric2D *mfn
) {
    slibreal_t maxerrJ[3][3], minerrJ[3][3],
        Rmin, Rmax, Zmin, Zmax,
        r, z, p, Delta, maxerr=0;
    unsigned int i, j, k, maxj=0, maxk=0;

    Rmin = mfa->GetRMajor() - mfa->GetRMinor();
    Rmax = mfa->GetRMajor() + mfa->GetRMinor();
    Zmin =-mfa->GetRMinor();
    Zmax = mfa->GetRMinor();

    for (i = 0; i < NPOINTS; i++) {
        /*
        r = Rand(Rmin, Rmax);
        p = Rand(0.0, 2.0*M_PI);
        z = Rand(Zmin, Zmax);
        */
        r = Rmin + (Rmax-Rmin)*((slibreal_t)i) / ((slibreal_t)(NPOINTS-1));
        p =        (2.0*PI)   *((slibreal_t)i) / ((slibreal_t)(NPOINTS-1));
        z = Zmin + (Zmax-Zmin)*((slibreal_t)i) / ((slibreal_t)(NPOINTS-1));

        struct magnetic_field_data ma = mfa->EvalDerivatives(r*cos(p), r*sin(p), z);
        struct magnetic_field_data mn = mfn->EvalDerivatives(r*cos(p), r*sin(p), z);

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                if (ma.J[j][k] != 0)
                    Delta = fabs((ma.J[j][k]-mn.J[j][k])/ma.J[j][k]);
                else
                    Delta = fabs(ma.J[j][k]-mn.J[j][k]);

                if (i == 0) {
                    maxerrJ[j][k] = Delta;
                    minerrJ[j][k] = Delta;
                    maxerr = Delta;
                } else if (Delta > maxerrJ[j][k])
                    maxerrJ[j][k] = Delta;
                else if (Delta < minerrJ[j][k])
                    minerrJ[j][k] = Delta;

                if (maxerrJ[j][k] > maxerr) {
                    maxerr = maxerrJ[j][k];
                    maxj = j;
                    maxk = k;
                }
            }
        }
    }

    if (maxerr > 3e-3) {
        this->PrintError("Error in derivatives is larger than expected (%u, %u). Delta = %e.", maxj, maxk, maxerr);
        return false;
    }

    /*
    bool printJ = true;
    if (printJ) {
        printf("Maximum errors in J:\n");
        printf("  %10.3e  %10.3e  %10.3e\n", maxerrJ[0][0], maxerrJ[0][1], maxerrJ[0][2]);
        printf("  %10.3e  %10.3e  %10.3e\n", maxerrJ[1][0], maxerrJ[1][1], maxerrJ[1][2]);
        printf("  %10.3e  %10.3e  %10.3e\n\n", maxerrJ[2][0], maxerrJ[2][1], maxerrJ[2][2]);

        printf("Minimum errors in J:\n");
        printf("  %10.3e  %10.3e  %10.3e\n", minerrJ[0][0], minerrJ[0][1], minerrJ[0][2]);
        printf("  %10.3e  %10.3e  %10.3e\n", minerrJ[1][0], minerrJ[1][1], minerrJ[1][2]);
        printf("  %10.3e  %10.3e  %10.3e\n\n", minerrJ[2][0], minerrJ[2][1], minerrJ[2][2]);
    } else
        printf("Maximum error in J: %.3e (%u,%u).\n", maxerr, maxj, maxk);
    */
    
    return true;
}

/**
 * Run the test.
 * RETURNS true on success, false on failure.
 * On failure, an error message is also printed.
 */
bool Test_MagneticFieldNumeric2D::Run() {
	slibreal_t
		B0       = magnetic_field_test_data_B0,
		Rm       = magnetic_field_test_data_Rm,
		rminor   = magnetic_field_test_data_rminor,
		qa_const = magnetic_field_test_data_qa_const;

	MagneticFieldNumeric2D *mf, *mff;
	MagneticFieldAnalytical2D *mfa;

	// Construct the analytical magnetic field object
	mfa = new MagneticFieldAnalytical2D(B0, Rm, rminor, MFATFS_CW, MFASF_CONSTANT, qa_const, 0.0);
	
	// Generate the numeric magnetic field
	mf = GenerateMF(mfa, Rm, rminor, NR, NZ, NWALL);

	try {
		ComparePoints(magnetic_field_test_data_const, mf, B0, "constant");
	} catch (SOFTLibException &ex) {
		this->PrintError("[MagneticFieldNumeric2D]: "+ex.whats());
		return false;
	}

	this->PrintOK("Magnetic field evaluation");

	// Load magnetic field from file
	try {
		mff = new MagneticFieldNumeric2D(TESTFILE);
		try {
			ComparePoints(magnetic_field_test_data_const, mff, B0, "constant");
		} catch (SOFTLibException &ex) {
			this->PrintError("[MagneticFieldNumeric2D]: "+ex.whats());
			return false;
		}
	} catch (SFileException& ex) {
		this->PrintWarning("[MagneticFieldNumeric2D]: Unable to open '" TESTFILE "'. Skipping test... "+ex.whats());
	}

    this->PrintOK("Magnetic field evaluation -- compared to table");

    // Check numeric magnetic field derivative scaling
    MagneticFieldNumeric2D *mfn = mfa->ToNumeric2D(200, 200);
    if (!CheckNumericDerivatives(mfa, mfn))
        return false;
    else
        this->PrintOK("Magnetic field derivative evaluation");

	return true;
}

