#ifndef _VECTOR_H
#define _VECTOR_H

#include <softlib/config.h>

template<unsigned int N>
class Vector {
	protected:
		const unsigned int n = N;
		slibreal_t elems[N];
	public:
		Vector<N>();
		Vector<N>(const slibreal_t*);
		Vector<N>(const Vector<N>&);
		~Vector<N>();

		slibreal_t& operator[](const unsigned int);
		const slibreal_t& operator[](const unsigned int) const;

		/* Arithmetic */
		Vector<N>& operator=(const Vector<N>&);
		Vector<N>& operator=(const slibreal_t*);
		Vector<N>& operator=(const slibreal_t);
		Vector<N>& operator+=(const Vector<N>&);
		Vector<N>& operator+=(const slibreal_t);
		Vector<N>& operator-=(const Vector<N>&);
		Vector<N>& operator-=(const slibreal_t);
		Vector<N>& operator*=(const slibreal_t);
		Vector<N>& operator/=(const slibreal_t);
        Vector<N>  operator-();
		slibreal_t Dot(const Vector<N>&) const;

		/* Useful operations */
		slibreal_t Norm() const;
		void Normalize();
		void ToArray(slibreal_t arr[N]) const;

        static Vector<N> Cross(const Vector<N>&, const Vector<N>&);
        static Vector<N> &Cross(const Vector<N>&, const Vector<N>&, Vector<N>&);

		unsigned int size() const;
};

template<>
Vector<3> Vector<3>::Cross(const Vector<3>&, const Vector<3>&);
template<>
Vector<3> &Vector<3>::Cross(const Vector<3>&, const Vector<3>&, Vector<3>&);

/**********************
 * OPERATOR OVERLOADS *
 **********************/
template<unsigned int N> Vector<N> operator+(const Vector<N>&, const Vector<N>&);
template<unsigned int N> Vector<N> operator+(const Vector<N>&, const slibreal_t);
template<unsigned int N> Vector<N> operator+(const slibreal_t, const Vector<N>&);
template<unsigned int N> Vector<N> operator-(const Vector<N>&, const Vector<N>&);
template<unsigned int N> Vector<N> operator-(const Vector<N>&, const slibreal_t);
template<unsigned int N> Vector<N> operator-(const slibreal_t, const Vector<N>&);
template<unsigned int N> Vector<N> operator*(const Vector<N>&, const slibreal_t);
template<unsigned int N> Vector<N> operator*(const slibreal_t, const Vector<N>&);
template<unsigned int N> Vector<N> operator/(const Vector<N>&, const slibreal_t);

// Implementation
#include <softlib/General/Vector.tcc>

#endif/*_VECTOR_H*/
