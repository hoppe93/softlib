/**
 * Implementation of a distribution function.
 */

#include <cmath>
#include <cstdio>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/RadialProfile.h>
#include <softlib/SOFTLibException.h>

/**
 * Create a new radial distribution function (i.e. a distribution
 * function with both radial and momentum dimensions) from this
 * momentum-space distribution function with the given radial profile.
 * If no radial profile is given, a uniform radial profile is used.
 *
 * (optional) rp: Radial profile object that modulates the distribution
 *    function in radius.
 *
 * RETURNS a newly allocated RadialDistributionFunction object which
 * contains pointers to this MomentumSpaceDistributionFunction as well
 * as the radial profile object given as argument, or a completely new
 * radial profile object if none was passed to this method.
 */
RadialDistributionFunction *MomentumSpaceDistributionFunction::ToRadialDistribution() {
    RadialDistributionFunction *rdf = new RadialDistributionFunction();
    UniformRadialProfile *urp = new UniformRadialProfile();

    rdf->Initialize(urp, this);
    return rdf;
}
RadialDistributionFunction *MomentumSpaceDistributionFunction::ToRadialDistribution(RadialProfile *rp) {
    RadialDistributionFunction *rdf = new RadialDistributionFunction();

    rdf->Initialize(rp, this);
    return rdf;
}
