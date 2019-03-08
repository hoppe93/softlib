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
        template<class S = T> Vector<N,T>(const Vector<N,S>&);
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

/*template<class T = slibreal_t>
Vector<3,T> Vector<3,T>::Cross(const Vector<3,T>&, const Vector<3,T>&);
template<class T = slibreal_t>
Vector<3,T> &Vector<3,T>::Cross(const Vector<3,T>&, const Vector<3,T>&, Vector<3,T>&);*/

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
