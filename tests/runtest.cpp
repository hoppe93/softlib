/* Unit tests */

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "runtest.h"

#include <softlib/config.h>

// Tests
#include "configuration/configuration_test.h"
#include "distribution/Test_AnalyticalAvalanche.h"
#include "distribution/Test_CODEDistributionFunction.h"
#include "distribution/Test_NumericMomentumSpaceDistribution.h"
#include "distribution/Test_SOFTDistributionFunction.h"
#include "magnetic/domain.h"
#include "integrator/rkdp45.h"
#include "magnetic/analytical2d.h"
#include "magnetic/numeric2d.h"
#ifdef SOFT_HDF5
#	include "sfile/hdf5_test.h"
#endif
#include "sfile/matlab_test.h"
#include "sfile/sdt_test.h"

using namespace std;

vector<UnitTest*> tests;

void add_test(UnitTest *t) {
	tests.push_back(t);
}
void init() {
	add_test(new Test_Configuration("configuration"));
    add_test(new Test_AnalyticalAvalanche("distribution_analytical"));
    add_test(new Test_CODEDistributionFunction("distribution_code"));
    add_test(new Test_NumericMomentumSpaceDistribution("distribution_numericms"));
    add_test(new Test_SOFTDistributionFunction("distribution_soft"));
	add_test(new Test_Domain("domain"));
	add_test(new Test_RKDP45("int_rkdp45"));
	add_test(new Test_MagneticFieldAnalytical2D("magnetic_analytical2d"));
	add_test(new Test_MagneticFieldNumeric2D("magnetic_numeric2d"));
#ifdef SOFT_HDF5
	add_test(new Test_SFile_HDF5("sfile_hdf5"));
#endif
	add_test(new Test_SFile_MAT("sfile_matlab"));
    add_test(new Test_SFile_SDT("sfile_sdt"));
}

int has_test(const char *name) {
	size_t i;
	for (i = 0; i < tests.size(); i++) {
		if (tests[i]->HasName(name))
			return i;
	}

	return -1;
}

int run_test(int index) {
    cout << "\e[1m:: " << tests[index]->GetName() << "\e[0m" << endl;
    try {
        if (tests[index]->Run()) {
            cout << "   \e[1;32m[SUCCESS]\e[0m Test '" << tests[index]->GetName() << "' completed successfully." << endl;
            return 0;
        } else {
            cout << "   \e[1;31m[FAIL]\e[0m    Test '" << tests[index]->GetName() << "' failed." << endl;
            return 1;
        }
    } catch (SOFTLibException& e) {
        cout << "   \e[1;31m[ERROR]\e[0m   --> " << e.whats() << endl;
        cout << "   \e[1;31m[FAIL]\e[0m    Test '" << tests[index]->GetName() << "' unexpectedly threw exception." << endl;
    }

    return 1;
}
int run_all_tests() {
	size_t failed = 0, i;

	for (i = 0; i < tests.size(); i++) {
		failed += run_test(i);
	}

	return failed;
}

UnitTest::UnitTest(const string& name) { this->name = name; }
string& UnitTest::GetName() { return this->name; }
bool UnitTest::HasName(const string& cmp) { return (this->name==cmp); }
void UnitTest::PrintError(const string& s, ...) {
	//cerr << "\e[1;31m[ERROR]\e[0m   " << s << endl;
    va_list args;
    va_start(args, s);

    fprintf(stderr, "   \e[1;31m[ERROR]\e[0m   ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
void UnitTest::PrintOK(const string& s, ...) {
	//cout << "\e[1;93m[OK]\e[0m      --> " << s << endl;
    va_list args;
    va_start(args, s);

    fprintf(stderr, "   \e[1;92m[OK]\e[0m      --> ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
void UnitTest::PrintStatus(const string& s, ...) {
	//cout << s << endl;
    va_list args;
    va_start(args, s);

    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
void UnitTest::PrintWarning(const string& s, ...) {
	//cout << "\e[1;93m[WARNING]\e[0m " << s << endl;
    va_list args;
    va_start(args, s);

    fprintf(stderr, "   \e[1;93m[WARNING]\e[0m ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
//bool UnitTest::Run() { return false; }

/**
 * Random number generator.
 */
void UnitTest::InitRand() {
    srand(time(NULL));
}
slibreal_t UnitTest::Rand() {
    return ((slibreal_t)rand())/((slibreal_t)RAND_MAX);
}
slibreal_t UnitTest::Rand(slibreal_t min, slibreal_t max) {
    return (min + (max-min)*Rand());
}

/**
 * Prints a brief help message.
 */
void help() {
    printf(
        "Unit tests for SOFTLib.\n\n"

        "Usage:\n"
        "    softlib_tests          Show this help message.\n"
        "    softlib_tests all      Run all tests.\n"
        "    softlib_tests [test1 [test2 [...]]]\n"
        "                           Runs the tests with names 'test1', 'test2' etc.\n\n"

        "Available tests:\n"
    );

	for (unsigned int i = 0; i < tests.size(); i++) {
        printf("    %s\n", tests[i]->GetName().c_str());
	}

    printf("\n");
}

int main(int argc, char *argv[]) {
	int i, t, failed, total;

	init();

	if (argc == 1) {
        help();
        return 0;
    } else if (argc == 2 && !strcmp(argv[1], "all")) {
		failed = run_all_tests();
        total = tests.size();
	} else {
        total = argc-1;
		for (i = 1, failed=0; i < argc; i++) {
			if ((t=has_test(argv[i]))>=0) {
				failed += run_test(t);
			} else {
				cerr << "\e[1;31m[ERROR]\e[0m   Unrecognized test: '" << argv[i] << "'." << endl;
				failed++;
			}
		}
	}

	cout << endl << (total-failed) << " of " << total << " tests passed." << endl;

	return failed;
}
