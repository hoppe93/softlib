/**
 * Implementation of the SOFTLib Vector type.
 */

#include <cmath>
#include <softlib/SOFTLibException.h>
#include <softlib/Vector.h>

/**
 * Constructor
 */
template<unsigned int N> Vector<N>::Vector() {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = 0.0;
}
template<unsigned int N> Vector<N>::Vector(const slibreal_t *init) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = init[i];
}
template<unsigned int N> Vector<N>::Vector(const Vector<N>& v) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = v[i];
}

/**
 * Destructor
 */
template<unsigned int N> Vector<N>::~Vector() { }

/**************
 * ASSIGNMENT *
 **************/
/**
 * Assign vector.
 * This method assigns values to
 * the vector component by component.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator=(const Vector<N>& rhs) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = rhs[i];
	
	return *this;
}
/**
 * Assign to vector from array.
 * The assigned value must contain
 * at least N elements (this is not
 * checked by this function). Elements
 * are copied to this object.
 *
 * rhs: Array of floating-point numbers
 *   to copy to this Vector.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator=(const slibreal_t *rhs) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = rhs[i];
	
	return *this;
}
/**
 * Assignment of scalar to 1D Vector.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator=(const slibreal_t rhs) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] = rhs;

	return *this;
}

/**
 * Get element of vector.
 *
 * i: Index of element to return.
 */
template<unsigned int N>
slibreal_t& Vector<N>::operator[](const unsigned int i) { return elems[i]; }
template<unsigned int N>
const slibreal_t& Vector<N>::operator[](const unsigned int i) const { return elems[i]; }

/**
 * Add and assign operator.
 *
 * v: Vector to add to this.
 *   or
 * s: Scalar to add to each element.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator+=(const Vector<N>& v) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] += v[i];

	return *this;
}
template<unsigned int N>
Vector<N>& Vector<N>::operator+=(const slibreal_t s) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] += s;

	return *this;
}

/**
 * Subtract and assign operator
 * 
 * v: Vector to subtract from this.
 *   or
 * s: Scalar to subtract from each element.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator-=(const Vector<N>& v) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] -= v[i];

	return *this;
}
template<unsigned int N>
Vector<N>& Vector<N>::operator-=(const slibreal_t s) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] -= s;

	return *this;
}

/**
 * Multiply by scalar and assign.
 *
 * s: Scalar value to multiply this vector with.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator*=(const slibreal_t s) {
	for (unsigned int i = 0; i < N; i++)
		elems[i] *= s;

	return *this;
}

/**
 * Divide by scalar and assign.
 *
 * s: Scalar value to multiply this vector with.
 */
template<unsigned int N>
Vector<N>& Vector<N>::operator/=(const slibreal_t s) {
	(*this) *= (1.0/s);
	return *this;
}

/**
 * Negate all elements.
 */
template<unsigned int N>
Vector<N> Vector<N>::operator-() {
    Vector<N> v;
    for (unsigned int i = 0; i < N; i++)
        v[i] = -elems[i];

    return v;
}

/**
 * Scalar product between this and v.
 *
 * v: Vector to scalar multiply with.
 */
template<unsigned int N>
slibreal_t Vector<N>::Dot(const Vector<N>& v) const {
	slibreal_t s = 0;
	for (unsigned int i = 0; i < N; i++)
		s += elems[i] * v[i];
	
	return s;
}

/**
 * Calculate the 2-norm of this vector.
 */
template<unsigned int N>
slibreal_t Vector<N>::Norm() const {
	slibreal_t s = 0;
	for (unsigned int i = 0; i < N; i++)
		s += elems[i] * elems[i];
	
	return std::sqrt(s);
}

/**
 * Normalize this vector so that its
 * norm is unity.
 */
template<unsigned int N> void Vector<N>::Normalize() { (*this) /= this->Norm(); }

/**
 * Return the number of elements in this vector.
 */
template<unsigned int N> unsigned int Vector<N>::size() const { return n; }

/**
 * Convert this class to a pure C++ array.
 * This function will assign the elements of
 * the vector to the first N elements of the
 * given array.
 *
 * arr: Array, expected to be of size N (though
 *   larger arrays are perfectly okay), to which
 *   the contents of this vector should be copied.
 */
template<unsigned int N> void Vector<N>::ToArray(slibreal_t *arr) const {
	unsigned int i;
	for (i = 0; i < N; i++)
		arr[i] = elems[i];
}


/********************************
 * RELATED NON-MEMBER FUNCTIONS *
 ********************************/
/**
 * Add two vectors.
 */
template<unsigned int N>
Vector<N> operator+(const Vector<N>& lhs, const Vector<N>& rhs) {
	Vector<N> v = Vector<N>(lhs);
	v += rhs;
	return v;
}
template<unsigned int N>
Vector<N> operator+(const Vector<N>& lhs, const slibreal_t rhs) {
	Vector<N> v = Vector<N>(lhs);
	v += rhs;
	return v;
}
template<unsigned int N>
Vector<N> operator+(const slibreal_t lhs, const Vector<N>& rhs) { return (rhs + lhs); }

/**
 * Subtract two vectors.
 */
template<unsigned int N>
Vector<N> operator-(const Vector<N>& lhs, const Vector<N>& rhs) {
	Vector<N> v = Vector<N>(lhs);
	v -= rhs;
	return v;
}
template<unsigned int N>
Vector<N> operator-(const Vector<N>& lhs, const slibreal_t rhs) {
	Vector<N> v = Vector<N>(lhs);
	v -= rhs;
	return v;
}
template<unsigned int N>
Vector<N> operator-(const slibreal_t lhs, const Vector<N>& rhs) {
	Vector<N> v;

	for (unsigned int i = 0; i < N; i++)
		v[i] = lhs - v[i];
	
	return v;
}

/**
 * Multiply vector by scalar
 */
template<unsigned int N>
Vector<N> operator*(const Vector<N>& lhs, const slibreal_t rhs) {
	Vector<N> v = Vector<N>(lhs);
	v *= rhs;
	return v;
}
template<unsigned int N>
Vector<N> operator*(const slibreal_t lhs, const Vector<N>& rhs) { return (rhs * lhs); }

/**
 * Divide vector by scalar
 */
template<unsigned int N>
Vector<N> operator/(const Vector<N>& lhs, const slibreal_t rhs) {
	Vector<N> v = Vector<N>(lhs);
	v /= rhs;
	return v;
}

/**
 * Evaluate cross-product between two vectors.
 * NOTE: Only defined for Vector<3>.
 */
template<unsigned int N>
Vector<N> Vector<N>::Cross(const Vector<N>&, const Vector<N>&) {
    throw SOFTLibException("Cross products are only defined for 3-vectors.");
}
template<unsigned int N>
Vector<N> &Vector<N>::Cross(const Vector<N>&, const Vector<N>&, Vector<N>&) {
    throw SOFTLibException("Cross products are only defined for 3-vectors.");
}

