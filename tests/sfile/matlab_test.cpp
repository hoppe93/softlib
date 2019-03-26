
#include <iostream>
#include <runtest.h>
#include <softlib/SFile_MAT.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>
#include <cmath>
#include <limits>
#include "matlab_test.h"
#include "sfile_test.h"

using namespace std;

Test_SFile_MAT::Test_SFile_MAT(const string& name) : UnitTest(name) {}

bool Test_SFile_MAT::Run() {
	bool success = true;

	try {
		SFile *sf = new SFile_MAT();
#ifdef OFFICIAL_MATLAB
		success = sfile_test(sf, "sfile_matlab_test.mat", true, false);
#else
		success = sfile_test(sf, "sfile_matlab_test.mat", true, true);
#endif
	} catch (SOFTLibException& ex) {
		this->PrintError("SFile_MAT: "+ex.whats());
		success = false;
	}

	return success;
}

