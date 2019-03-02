#ifndef _VECTOR_H
#define _VECTOR_H

#include <softlib/config.h>

template<unsigned int N, class T = slibreal_t>
class Vector {
	protected:
		const unsigned int n = N;
		T elems[N];
	public:
		Vector<N,T>();
		Vector<N,T>(const T*);
		Vector<N,T>(const Vector<N,T>&);
		~Vector<N,T>();

		T& operator[](const unsigned int);
		const T& operator[](const unsigned int) const;

		/* Arithmetic */
		Vector<N,T>& operator=(const Vector<N,T>&);
		Vector<N,T>& operator=(const T*);
		Vector<N,T>& operator=(const T&);
		Vector<N,T>& operator+=(const Vector<N,T>&);
		Vector<N,T>& operator+=(const T&);
		Vector<N,T>& operator-=(const Vector<N,T>&);
		Vector<N,T>& operator-=(const T&);
		Vector<N,T>& operator*=(const T&);
		Vector<N,T>& operator/=(const T&);
        Vector<N,T>  operator-();
		T Dot(const Vector<N,T>&) const;

		/* Useful operations */
		T Norm() const;
		void Normalize();
		void ToArray(T arr[N]) const;

        static Vector<N,T> Cross(const Vector<N,T>&, const Vector<N,T>&);
        static Vector<N,T> &Cross(const Vector<N,T>&, const Vector<N,T>&, Vector<N,T>&);

		unsigned int size() const;
};

template<>
Vector<3> Vector<3>::Cross(const Vector<3>&, const Vector<3>&);
template<>
Vector<3> &Vector<3>::Cross(const Vector<3>&, const Vector<3>&, Vector<3>&);

/**********************
 * OPERATOR OVERLOADS *
 **********************/
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator+(const Vector<N,T>&, const Vector<N,T>&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator+(const Vector<N,T>&, const T&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator+(const T&, const Vector<N,T>&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator-(const Vector<N,T>&, const Vector<N,T>&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator-(const Vector<N,T>&, const T&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator-(const T&, const Vector<N,T>&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator*(const Vector<N,T>&, const T&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator*(const T&, const Vector<N,T>&);
template<unsigned int N, class T = slibreal_t> Vector<N,T> operator/(const Vector<N,T>&, const T&);

// Implementation
#include <softlib/General/Vector.tcc>

#endif/*_VECTOR_H*/
