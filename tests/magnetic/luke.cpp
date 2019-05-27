/**
 * Test of the 'LUKE' magnetic field format.
 */

#include <runtest.h>
#include "luke.h"
#include "magfield_points.h"

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldLUKE.h>
#include <softlib/SFileException.h>

#define THRESHOLD 			1e-3
#define TESTFILE_LUKE 		"test_equil_luke.mat"
#define TESTFILE_LUKE_WALL	"test_equil_luke_wall.mat"
#define NR					99
#define NZ					95
#define NTHETA				93
#define NWALL				100

Test_MagneticFieldLUKE::Test_MagneticFieldLUKE(const string& name) : Test_MagneticField(name, THRESHOLD) {}

void Test_MagneticFieldLUKE::GenerateLUKEMF(
	MagneticFieldAnalytical2D *mfa, slibreal_t Rm, slibreal_t rminor,
	unsigned int nr, unsigned int nz,
	unsigned int ntheta, unsigned int nwall
) {
	const double magaxis_roffs = Rm;
	const double magaxis_zoffs = 0;

	// Construct R/Z grid
	double *R = new double[nr];
	double *Z = new double[nz];

	for (unsigned int i = 0; i < nr; i++)
		R[i] = rminor*((((double)(2*i))/(nr-1)) - 1);
	for (unsigned int i = 0; i < nz; i++)
		Z[i] = rminor*((((double)(2*i))/(nz-1)) - 1);
	
	// Generate magnetic flux and poloidal angle
	double **Psi = new double*[nz];
	Psi[0] = new double[nr*nz];
	for (unsigned int i = 1; i < nz; i++)
		Psi[i] = Psi[i-1] + nr;

	// Poloidal flux on LCFS
	double psia = mfa->EvalFlux(Rm + rminor, 0.0, 0.0);

	for (unsigned int j = 0; j < nz; j++) {
		for (unsigned int i = 0; i < nr; i++) {
			Psi[j][i] = mfa->EvalFlux(R[i] + Rm, 0.0, Z[j]) / psia;
		}
	}

	// XXX Assume circular flux surfaces
	unsigned int
		npsi = nr/2 + nr%2,
		nrbase = nr/2,
		nzbase = nz/2;
	double *pgrid = new double[npsi];
	double *tgrid = new double[ntheta];

	for (unsigned int i = 0; i < npsi; i++)
		pgrid[i] = Psi[nzbase][nrbase + i] * psia;
	
	for (unsigned int i = 0; i < ntheta; i++)
		tgrid[i] = 2.0*M_PI * ((double)i)/((double)(ntheta-1));
	
	double **ptx = new double*[ntheta];
	double **pty = new double*[ntheta];
	ptx[0] = new double[ntheta*npsi];
	pty[0] = new double[ntheta*npsi];

	for (unsigned int i = 1; i < ntheta; i++) {
		ptx[i] = ptx[i-1] + npsi;
		pty[i] = pty[i-1] + npsi;
	}

	for (unsigned int i = 0; i < ntheta; i++) {
		double ct = cos(tgrid[i]);
		double st = sin(tgrid[i]);
		for (unsigned int j = 0; j < npsi; j++) {
			ptx[i][j] = R[nrbase + j] * ct;
			pty[i][j] = R[nrbase + j] * st;
		}
	}

	// XXX assumption done

	// Evaluate magnetic field on psi/theta grid
	double **Br = new double*[ntheta];
	double **Bt = new double*[ntheta];
	double **Bz = new double*[ntheta];

	Br[0] = new double[ntheta*npsi];
	Bt[0] = new double[ntheta*npsi];
	Bz[0] = new double[ntheta*npsi];
	for (unsigned int i = 1; i < ntheta; i++) {
		Br[i] = Br[i-1] + npsi;
		Bt[i] = Bt[i-1] + npsi;
		Bz[i] = Bz[i-1] + npsi;
	}

	for (unsigned int i = 0; i < ntheta; i++) {
		for (unsigned int j = 0; j < npsi; j++) {
			slibreal_t *B = mfa->Eval(ptx[i][j] + Rm, 0.0, pty[i][j]);
			Br[i][j] = B[0];
			Bt[i][j] = B[1];
			Bz[i][j] = B[2];
		}
	}

	// Generate wall
	double **wall = new double*[nwall];
	wall[0] = new double[2*nwall];
	wall[1] = wall[0] + nwall;

	for (unsigned int i = 0; i < nwall; i++) {
		double t = 2.0*M_PI * ((slibreal_t)i) / ((slibreal_t)(nwall-1));
		wall[0][i] = Rm + rminor*cos(t);
		wall[1][i] = rminor*sin(t);
	}

	// Write magnetic field
	SFile *sf = SFile::Create(TESTFILE_LUKE, SFILE_MODE_WRITE);

	sf->CreateStruct("equil");
	sf->WriteString("equil/id", "EQUIL_LUKE_TEST");
	sf->WriteList("equil/psi_apRp", pgrid, npsi);
	sf->WriteList("equil/theta", tgrid, ntheta);
	sf->WriteArray("equil/ptx", ptx, ntheta, npsi);
	sf->WriteArray("equil/pty", pty, ntheta, npsi);
	sf->WriteArray("equil/ptBx", Br, ntheta, npsi);
	sf->WriteArray("equil/ptBy", Bz, ntheta, npsi);
	sf->WriteArray("equil/ptBPHI", Bt, ntheta, npsi);

	sf->WriteScalar("equil/Rp", magaxis_roffs);
	sf->WriteScalar("equil/Zp", magaxis_zoffs);

	sf->CreateStruct("equil/eqdsk");
	sf->WriteList("equil/eqdsk/x", R, nr);
	sf->WriteList("equil/eqdsk/y", Z, nz);
	sf->WriteArray("equil/eqdsk/xypsin", Psi, nz ,nr);
	sf->WriteScalar("equil/eqdsk/psia", psia);

	sf->Close();

	// Write wall
	sf = SFile::Create(TESTFILE_LUKE_WALL, SFILE_MODE_WRITE);
	sf->WriteArray("wall", wall, 2, nwall);
	sf->Close();
}

bool Test_MagneticFieldLUKE::Run() {
	slibreal_t
		B0       = magnetic_field_test_data_B0,
		Rm       = magnetic_field_test_data_Rm,
		zm       = magnetic_field_test_data_zm,
		rminor   = magnetic_field_test_data_rminor,
		qa_const = magnetic_field_test_data_qa_const;
	
	MagneticFieldAnalytical2D *mfa = new MagneticFieldAnalytical2D(B0, Rm, zm, rminor, MFAFS_CW, MFAFS_CCW, MFASF_CONSTANT, qa_const, 0.0);
	GenerateLUKEMF(mfa, Rm, rminor, NR, NZ, NTHETA, NWALL);
	bool success = true;

	try {
		MagneticFieldLUKE *mf = new MagneticFieldLUKE(TESTFILE_LUKE, TESTFILE_LUKE_WALL);
		try {
			ComparePoints(magnetic_field_test_data_const, mf, B0, "constant", false);
		} catch (SOFTLibException &ex) {
			this->PrintError("[MagneticFieldLUKE]: %s", ex.what());
			return false;
		}
	} catch (SFileException& ex) {
		this->PrintError("[MagneticFieldLUKE]: Unable to open '%s' and '%s'. Failing test...\n%s", TESTFILE_LUKE, TESTFILE_LUKE_WALL, ex.what());
		success = false;
	}

	if (success) {
		this->PrintOK("Magnetic field evaluation -- compared to table");

		// Delete test files
		remove(TESTFILE_LUKE_WALL);
		remove(TESTFILE_LUKE);
	}

	return success;
}

