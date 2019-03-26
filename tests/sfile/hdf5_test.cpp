
#include <iostream>
#include <runtest.h>
#include <softlib/SFile_HDF5.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>
#include <cmath>
#include <limits>
#include "hdf5_test.h"
#include "sfile_test.h"
#include "hdf5.h"

using namespace std;

Test_SFile_HDF5::Test_SFile_HDF5(const string& name) : UnitTest(name) {}

bool Test_SFile_HDF5::Run() {
	bool success = true;

	try {
		SFile *sf = new SFile_HDF5();
		success = sfile_test(sf, "sfile_hdf5_test.h5", true, true);
	} catch (SOFTLibException& ex) {
		this->PrintError("[SFile_HDF5]: "+ex.whats());
		success = false;
	}

	return success;
}

