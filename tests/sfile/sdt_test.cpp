
#include <iostream>
#include <runtest.h>
#include <softlib/SFile_SDT.h>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>
#include <cmath>
#include <limits>
#include "sdt_test.h"
#include "sfile_test.h"

using namespace std;

Test_SFile_SDT::Test_SFile_SDT(const string& name) : UnitTest(name) {}

bool Test_SFile_SDT::Run() {
	bool success = true;

	try {
		SFile *sf = new SFile_SDT();
		success = sfile_test(sf, "sfile_sdt_test.sdt", false);
	} catch (SOFTLibException& ex) {
		this->PrintError("SFile_SDT: "+ex.whats());
		success = false;
	}

	return success;
}

