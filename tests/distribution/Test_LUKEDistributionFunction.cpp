/**
 * Implementation of unit tests for the 'LUKEDistributionFunction'
 * class.
 */

#include <cmath>
#include <ctime>

#include <softlib/config.h>
#include <softlib/DistributionFunction/LUKEDistributionFunction.h>
#include "Test_NumericDistributionFunction.h"
#include "Test_LUKEDistributionFunction.h"

/**
 * Run the tests.
 */
bool Test_LUKEDistributionFunction::Run() {
    return RunInternal<LUKEDistributionFunction>(this->inputfilename);
}

