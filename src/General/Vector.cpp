/**
 * Implementation of the SOFTLib Vector type.
 */

#include <softlib/Vector.h>

/**
 * Specialization of cross-product.
 */
template<class T>
Vector<3,T> &Vector<3,T>::Cross(const Vector<3,T> &a, const Vector<3,T> &b, Vector<3,T>& c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

    return c;
}

template<class T>
Vector<3,T> Vector<3,T>::Cross(const Vector<3,T> &a, const Vector<3,T> &b) {
    Vector<3,T> c;
    return Cross(a, b, c);
}

