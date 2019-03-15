/**
 * Implementation of unit tests for the 'SOFTDistributionFunction'
 * class.
 */

#include <cmath>
#include <ctime>

#include <softlib/config.h>
#include <softlib/DistributionFunction/SOFTDistributionFunction.h>
#include "Test_NumericDistributionFunction.h"
#include "Test_SOFTDistributionFunction.h"

/**
 * Run the tests.
 */
bool Test_SOFTDistributionFunction::Run() {
    return RunInternal<SOFTDistributionFunction>(this->inputfilename);
}

