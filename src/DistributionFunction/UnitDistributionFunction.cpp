/**
 * Implementation of unit distribution function.
 */

#include <softlib/DistributionFunction/UnitDistributionFunction.h>

/**
 * Clone this distribution function.
 */
UnitDistributionFunction *UnitDistributionFunction::MinClone() {
    return new UnitDistributionFunction();
}
