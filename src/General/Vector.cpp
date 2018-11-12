/**
 * Implementation of the SOFTLib Vector type.
 */

#include <softlib/Vector.h>

/**
 * Specialization of cross-product.
 */
template<>
Vector<3> &Vector<3>::Cross(const Vector<3> &a, const Vector<3> &b, Vector<3>& c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

    return c;
}

template<>
Vector<3> Vector<3>::Cross(const Vector<3> &a, const Vector<3> &b) {
    Vector<3> c;
    return Cross(a, b, c);
}

