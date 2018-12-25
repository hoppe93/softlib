/**
 * Implementation of non-abstract methods
 * in the MagneticField2D class.
 */

#include <cmath>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructor
 */
MagneticField2D::MagneticField2D() {}
/**
 * Destructor
 */
MagneticField2D::~MagneticField2D() {}

/**
 * Checks whether the particle is inside the domain.
 * "Inside" is defined as the side of the domain
 * on which the magnetic axis lies, so that if the
 * line from the magnetic axis to the given point intersects
 * the domain, the particle has left the device.
 *
 * xyz/x,y,z: Point to check.
 *
 * RETURNS true if the point lies inside the domain.
 *   false otherwise.
 */
bool MagneticField2D::CrossesDomain(slibreal_t *xyz) {
    return CrossesDomain(magnetic_axis[0], 0.0, magnetic_axis[1], xyz[0], xyz[1], xyz[2]);
}
bool MagneticField2D::CrossesDomain(slibreal_t *xyz1, slibreal_t *xyz2) {
	return CrossesDomain(xyz1[0], xyz1[1], xyz1[2], xyz2[0], xyz2[1], xyz2[2]);
}
bool MagneticField2D::CrossesDomain(slibreal_t x, slibreal_t y, slibreal_t z) {
	return CrossesDomain(magnetic_axis[0], 0.0, magnetic_axis[1], x, y, z);
}
bool MagneticField2D::CrossesDomain(
	slibreal_t x1, slibreal_t y1, slibreal_t z1,
	slibreal_t x2, slibreal_t y2, slibreal_t z2
) {
    slibreal_t r = hypot(x2, y2),
        x00 = hypot(x1,y1), y00=z1,
        x01 = r - x00, y01 = z2 - y00,
        x10, y10, x11, y11, s, t, det;

    unsigned int i;
    for (i = 0; i < ndomain-1; i++) {
        x10 = rdomain[i];
        x11 = rdomain[i+1] - x10;
        y10 = zdomain[i];
        y11 = zdomain[i+1] - y10;

        /* Check if matrix is zero */
        if (x00-x10 == 0 && y00-y10 == 0) return false;
        
        /* Compute determinant and check if zero */
        det = x11*y01 - x01*y11;
        if (det == 0) return true;

        s = (1/det) * ((x00-x10)*y01 - (y00-y10)*x01);
        t =-(1/det) * ((y00-y10)*x11 - (x00-x10)*y11);

        /* If s and t are between 0 and 1 => intersection */
        if (s>=0 && s<=1 && t>=0 && t<=1) return false;
    }
    
    return true;
}

/**
 * Short-hands for evaluating magnetic field.
 *
 * xyz: Cartesian coordinate vector representing
 *    the point in which to evaluate the magnetic field.
 */
slibreal_t *MagneticField2D::Eval(slibreal_t *xyz) {
	return Eval(xyz[0], xyz[1], xyz[2]);
}
slibreal_t *MagneticField2D::Eval(Vector<3>& xyz) {
	return Eval(xyz[0], xyz[1], xyz[2]);
}

/**
 * Short-hands for evaluating derivatives of
 * the magnetic field.
 *
 * xyz: Cartesian coordinate vector representing
 *    the point in which to evaluate the magnetic field.
 */
struct magnetic_field_data& MagneticField2D::EvalDerivatives(slibreal_t *xyz) {
	return EvalDerivatives(xyz[0], xyz[1], xyz[2]);
}
struct magnetic_field_data& MagneticField2D::EvalDerivatives(Vector<3>& xyz) {
	return EvalDerivatives(xyz[0], xyz[1], xyz[2]);
}

/**
 * Short-hands for evaluating magnetic flux.
 *
 * xyz: Cartesian coordinate vector representing
 *      the point in which to evaluate the flux.
 */
slibreal_t MagneticField2D::EvalFlux(slibreal_t *xyz) {
    return EvalFlux(xyz[0], xyz[1], xyz[2]);
}
slibreal_t MagneticField2D::EvalFlux(Vector<3>& xyz) {
    return EvalFlux(xyz[0], xyz[1], xyz[2]);
}

/**
 * Short-hands for evaluating magnetix flux derivatives.
 *
 * xyz: Cartesian coordinate vector representing
 *      the point in which to evaluate the flux.
 */
struct flux_diff *MagneticField2D::EvalFluxDerivatives(slibreal_t *xyz) {
    return EvalFluxDerivatives(xyz[0], xyz[1], xyz[2]);
}
struct flux_diff *MagneticField2D::EvalFluxDerivatives(Vector<3> &xyz) {
    return EvalFluxDerivatives(xyz[0], xyz[1], xyz[2]);
}

/**
 * Returns the edges of a square which exactly bounds the entire
 * domain. This can be used as a first, cheaper domain check,
 * instead of using the relatively expensive 'CrossesDomain()'
 * routine.
 *
 * The function returns the bounds in the four input parameters:
 *
 * rmin: Minimum radius of domain.
 * rmax: Maximum radius of domain.
 * zmin: Minimum height of domain.
 * zmax: Maximum height of domain.
 */
void MagneticField2D::GetDomainBounds(slibreal_t &rmin, slibreal_t &rmax, slibreal_t &zmin, slibreal_t &zmax) {
    if (this->bounds_set) {
        rmin = this->bounds.rmin;
        rmax = this->bounds.rmax;
        zmin = this->bounds.zmin;
        zmax = this->bounds.zmax;
        return;
    }
        
    rmin = rmax = GetMagneticAxisR();
    zmin = zmax = GetMagneticAxisZ();

    for (unsigned int i = 0; i < ndomain; i++) {
        if (rdomain[i] < rmin)
            rmin = rdomain[i];
        else if (rdomain[i] > rmax)
            rmax = rdomain[i];

        if (zdomain[i] < zmin)
            zmin = zdomain[i];
        else if (zdomain[i] > zmax)
            zmax = zdomain[i];
    }

    this->bounds.rmin = rmin;
    this->bounds.rmax = rmax;
    this->bounds.zmin = zmin;
    this->bounds.zmax = zmax;
    this->bounds_set = true;
}

/**
 * Returns the radial position of the separatrix
 * (or wall, if the separatrix is not available)
 * at the vertical level of the magnetic axis.
 * This is the maximum radius at which a particle
 * could be launched in the plasma.
 */
slibreal_t MagneticField2D::GetMaxRadius() {
	if (this->maxradius > 0)
		return this->maxradius;
	else
		return this->FindMaxRadius();
}
slibreal_t MagneticField2D::FindMaxRadius() {
	return this->FindMaxRadius(this->ndomain, this->rdomain, this->zdomain);
}
slibreal_t MagneticField2D::GetMinRadius() {
	if (this->minradius >= 0)
		return this->minradius;
	else
		return this->FindMinRadius();
}
slibreal_t MagneticField2D::FindMinRadius() {
	return this->FindMinRadius(this->ndomain, this->rdomain, this->zdomain);
}

/**
 * Finds the maximum and minimum allowed particle
 * radii, given a domain. Returns only the former.
 *
 * n: Number of points in the domain.
 * r: R-coordinates of the domain.
 * z: Z-coordinates of the domain.
 *
 * RETURNS the maximum radius at which a particle
 * can be launched, on the vertical level of the
 * magnetic axis.
 */
slibreal_t MagneticField2D::FindMaxRadius(unsigned int n, slibreal_t *r, slibreal_t *z) {
	#define MAGFIELD2D_MAXRPOINTS 2
	unsigned int i, ir;
	slibreal_t za = this->GetMagneticAxisZ();
	slibreal_t rpoints[MAGFIELD2D_MAXRPOINTS];

	for (i = ir = 0; i < n-1 && ir < MAGFIELD2D_MAXRPOINTS; i++) {
		if ((z[i] >= za && z[i+1] <= za) ||
			(z[i] <= za && z[i+1] >= za)) {
			if (z[i+1] == z[i])
				rpoints[ir++] = max(r[i], r[i+1]);
			else
				rpoints[ir++] = r[i] + (r[i+1]-r[i]) * (za-z[i])/(z[i+1]-z[i]);
		}
	}

	if (ir == 0)
		throw SOFTLibException("Unable to find maximum allowed drop radius. Magnetic axis not within domain.");

	// Find max & min radii
	this->maxradius = rpoints[0];
	this->minradius = rpoints[0];
	for (i = 1; i < ir; i++) {
		if (rpoints[i] > this->maxradius)
			this->maxradius = rpoints[i];
		else if (rpoints[i] < this->minradius)
			this->minradius = rpoints[i];
	}

	return this->maxradius;
}

/**
 * Find the maximum and minimum z value of the domain
 * at the given radius.
 *
 * r:    Radius at which
 * zmin: Minimum value of Z at r.
 * zmax: Maximum value of Z at r.
 */
void MagneticField2D::FindMaxZ(const slibreal_t r, slibreal_t &zmin, slibreal_t &zmax) {
    FindMaxZ(ndomain, rdomain, zdomain, r, zmin, zmax);
}

/**
 * Finds the maximum and minimum allowed vertical
 * locations of the particle, at the given radius.
 *
 * n:    Number of points in the domain.
 * r:    R-coordinates of the domain.
 * z:    Z-coordinates of the domain.
 * R:    Radius at which to find zmin/zmax.
 * zmin: Contains minimum Z on return.
 * zmax: Contains maximum Z on return.
 */
void MagneticField2D::FindMaxZ(
    unsigned int n, const slibreal_t *r, const slibreal_t *z,
    const slibreal_t R, slibreal_t &zmin, slibreal_t &zmax
) {
    FindMaxMin(n, z, r, R, zmin, zmax);
}

/**
 * Find maximum and minimum value of 'x', at 'y' = 'a'.
 * The minimum and maximum values are returned in 'mn'
 * and 'mx'.
 *
 * n:  Number of points in 'x' and 'y'.
 * x:  X-vector of points (to find min/max of)
 * y:  Y-vector of points (to act as constraint)
 * a:  Value of 'Y' to hold as constraint.
 * mn: Contains minimum value of 'x' at 'y' = 'a' on return.
 * mx: Contains maximum value of 'x' at 'y' = 'a' on return.
 */
void MagneticField2D::FindMaxMin(
    unsigned int n, const slibreal_t *x, const slibreal_t *y,
    const slibreal_t a, slibreal_t &mn, slibreal_t &mx
) {
	#define MAGFIELD2D_MAXXPOINTS 2
	unsigned int i, ix;
	slibreal_t xpoints[MAGFIELD2D_MAXXPOINTS];

	for (i = ix = 0; i < n-1 && ix < MAGFIELD2D_MAXRPOINTS; i++) {
		if ((y[i] >= a && y[i+1] <= a) ||
			(y[i] <= a && y[i+1] >= a)) {
			if (y[i+1] == y[i])
				xpoints[ix++] = max(x[i], x[i+1]);
			else
				xpoints[ix++] = x[i] + (x[i+1]-x[i]) * (a-y[i])/(y[i+1]-y[i]);
		}
	}

	if (ix == 0)
		throw SOFTLibException("Unable to find maximum point. Requested point is not within the domain.");

	// Find max & min radii
	mx = xpoints[0];
	mn = xpoints[0];
	for (i = 1; i < ix; i++) {
		if (xpoints[i] > mx)
			mx = xpoints[i];
		else if (xpoints[i] < mn)
			mn = xpoints[i];
	}
}

/**
 * Finds the maximum and minimum allowed particle
 * radii, given a domain. Only returns the former.
 *
 * n: Number of points in the domain.
 * r: R-coordinates of the domain.
 * z: Z-coordinates of the domain.
 *
 * RETURNS the minimum radius at which a particle
 * can be launched, on the vertical level of the
 * magnetic axis.
 */
slibreal_t MagneticField2D::FindMinRadius(unsigned int n, slibreal_t *r, slibreal_t *z) {
	FindMaxRadius(n, r, z);
	return this->minradius;
}

/**
 * Get the location of the magnetic axis.
 * 'GetMagneticAxis' returns both the R and Z
 * coordinates in one call. 'GetMagneticAxisR'
 * and 'GetMagneticAxisZ' returns only the R
 * and Z coordinates respectively.
 */
slibreal_t *MagneticField2D::GetMagneticAxis()  { return this->magnetic_axis; }
slibreal_t  MagneticField2D::GetMagneticAxisR() { return this->magnetic_axis[0]; }
slibreal_t  MagneticField2D::GetMagneticAxisZ() { return this->magnetic_axis[1]; }

/**
 * Checks whether the given vector intersects
 * the domain in 3D.
 *
 * xyz0/x0,y0,z0: Origin of vector
 * xyz1/x1,y1,y1: Vector coordinates
 */
bool MagneticField2D::IntersectsDomain3D(slibreal_t *xyz0, slibreal_t *xyz1, bool has_outer_wall) {
    return IntersectsDomain3D(
        xyz0[0], xyz0[1], xyz0[1],
        xyz1[0], xyz1[1], xyz1[2],
        has_outer_wall
    );
}
bool MagneticField2D::IntersectsDomain3D(
    slibreal_t x0, slibreal_t y0, slibreal_t z0,
    slibreal_t x1, slibreal_t y1, slibreal_t z1,
    bool has_outer_wall
) {
    bool retval = false;
    slibreal_t dr, drp, dz, dzp;
    unsigned int i;
    for (i = 0; i < ndomain; i++) {
        if (i+1 < ndomain) {
            dr = rdomain[i], drp = rdomain[i+1];
            dz = zdomain[i], dzp = zdomain[i+1];
        } else {
            dr = rdomain[i], drp = rdomain[0];
            dz = zdomain[i], dzp = zdomain[0];
        }

        /* Ignore all outer wall tiles if requested */
        if (!has_outer_wall && dr > magnetic_axis[0] && drp > magnetic_axis[0])
            continue;

        if ((dz < z1 && dzp < z1 && dz < z0 && dzp < z0) ||
            (dz > z1 && dzp > z1 && dz > z0 && dzp > z0))
            continue;

        if (__IntersectsDomain3D(
            drp, dzp, dr, dz,
            x1, y1, z1,
            x0, y0, z0))
            retval = true;
    }

    return retval;
}

/**
 * Internal routine for check if a vector
 * intersects a part of the domain in 3D.
 */
bool MagneticField2D::__IntersectsDomain3D(
    slibreal_t xnp, slibreal_t znp, slibreal_t xn, slibreal_t zn,
    slibreal_t xp, slibreal_t yp, slibreal_t zp,
    slibreal_t xc, slibreal_t yc, slibreal_t zc
) {
    slibreal_t t1=0, t2=0, s1=0, s2=0,
        ic, term1, term2, square, sqr,
        izn, kz, mz, a, b, c, m, abc, iabc, term, q;

    if (znp == zn) {
        /* Equation becomes singular => line if sight does not intersect domain */
        if (zc == zp) return false;

        t1 = t2 = (zn-zp)/(zc-zp);
        ic = 1/(xn-xnp);
        term1 = xp + (xc-xp)*t1, term2 = yp + (yc-yp)*t1;
        square = term1*term1 + term2*term2;
        sqr = sqrt(square);

        s1 = ic * (xn + sqr);
        s2 = ic * (xn + sqr);

        if ((s1 >= 0 && s1 <= 1) || (s2 >= 0 && s2 <= 1))
            return true;
        else return false;
    } else {
        izn = 1 / (znp-zn);
        kz = (zc-zp)*izn;
        mz = (zp-zn)*izn;
        a=xc-xp, b=yc-yp, c=(xnp-xn)*kz;
        m = xn+mz*(xnp-xn);

        abc = a*a + b*b - c*c;
        if (abc == 0) return true;
        iabc = 1/abc;

        term = (c*m - a*xp - b*yp) * iabc;
        q = (m*m - xp*xp - yp*yp)*iabc;
        square = term*term + q;

        if (square < 0) return false;

        sqr = sqrt(square);
        t1 = term + sqr, t2 = term - sqr;
        s1 = kz*t1,      s2 = kz*t2;

        if ((t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1))
            return true;
    }

    return false;
}

void MagneticField2D::SetDomain(slibreal_t *rd, slibreal_t *zd, unsigned int nd) {
    this->rdomain = rd;
    this->zdomain = zd;
    this->ndomain = nd;

	this->FindMaxRadius();
}

