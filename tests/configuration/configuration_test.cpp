
#include <iostream>
#include <runtest.h>
#include <cmath>
#include <limits>

#include <softlib/Configuration.h>
#include <softlib/SOFTLibException.h>

#include "configuration_test.h"

using namespace std;

Test_Configuration::Test_Configuration(const string& name) : UnitTest(name) {}

const unsigned int Test_Configuration_NCONFIGS=3;
const string Test_Configuration_configs[Test_Configuration_NCONFIGS] = {
    // Config 1 (very basic)
    "option1 = val1;\n"
    "option2 = val2;",
    // Config 2 (basic with blocks)
    "option1 = block1;\n"
    "option2 = block2;\n\n"

    "@Block block1 {\n"
    "    option1 = val1;\n"
    "    option2 = val2;\n"
    "}\n\n"

    "@Block block2 {\n"
    "    option1 = val1;\n"
    "    option2 = val2;\n"
    "}",
    // Config 3 (include files)
    "option1 = block1;\n"
    "option2 = val2;\n\n"
    
    "<../tests/configuration/example1>\n\n"

    "@Block block1 {\n"
    "    option1 = val1;\n"
    "    option2 = val2;\n"
    "    <../tests/configuration/example2>\n"
    "}"
};

bool Test_Configuration::Run() {
	bool success = true;

    for (unsigned int i = 0; i < Test_Configuration_NCONFIGS && success; i++) {
        Configuration *cnf = new Configuration();
        cnf->RegisterBlockType("@Block");

        try {
            cnf->FromString(Test_Configuration_configs[i]);

            if (cnf->HasError()) {
                success = false;
                this->PrintError("Configuration test #%u failed.", i+1);
            }
        } catch (SOFTLibException& ex) {
            this->PrintError("Configuration test #%u failed with error '%s'", i+1, ex.what());
            success = false;
        }

        delete cnf;
    }

	return success;
}

