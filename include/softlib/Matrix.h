#ifndef _MATRIX_H
#define _MATRIX_H

//#define MATRIX_BOUNDS_CHECK

#include <softlib/config.h>
#include <softlib/Vector.h>

template<unsigned int M, unsigned int N, class T = slibreal_t>
class Matrix {
    protected:
        const unsigned int
            m = M, n = N;
        T elems[M*N];

    public:
        Matrix<M,N,T>();
        Matrix<M,N,T>(const T*);
        Matrix<M,N,T>(const T**);
        Matrix<M,N,T>(const Matrix<M,N,T>&);
        Matrix<M,N,T>(const Vector<M,T>&, const Vector<N,T>&);
        ~Matrix<M,N,T>();

        T& operator()(const unsigned int);
        const T& operator()(const unsigned int) const;
        T& operator()(const unsigned int, const unsigned int);
        const T& operator()(const unsigned int, const unsigned int) const;

        /* Arithmetic */
        Matrix<M,N,T>& operator=(const Matrix<M,N,T>&);
        Matrix<M,N,T>& operator+=(const Matrix<M,N,T>&);
        Matrix<M,N,T>& operator+=(const T);
        Matrix<M,N,T>& operator-=(const Matrix<M,N,T>&);
        Matrix<M,N,T>& operator-=(const T);
        Matrix<M,N,T>& operator*=(const Matrix<N,M,T>&);
        Matrix<M,N,T>& operator*=(const T);
        Matrix<M,N,T>& operator-();

        unsigned int Rows() const;
        unsigned int Cols() const;
};

/**********************
 * OPERATOR OVERLOADS *
 **********************/
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator+(const Matrix<M,N,T>&, const Matrix<M,N,T>&);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator+(const Matrix<M,N,T>&, const T);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator+(const T,              const Matrix<M,N,T>&);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator-(const Matrix<M,N,T>&, const Matrix<M,N,T>&);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator-(const Matrix<M,N,T>&, const T);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator-(const T,              const Matrix<M,N,T>&);
template<unsigned int M, unsigned int N, unsigned int K, class T> Matrix<M,K,T> operator*(const Matrix<M,N,T>&, const Matrix<N,K,T>&);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator*(const Matrix<M,N,T>&, const T);
template<unsigned int M, unsigned int N, class T> Matrix<M,N,T> operator*(const T,              const Matrix<M,N,T>&);
template<unsigned int M, unsigned int N, class T> Vector<M,T>   operator*(const Matrix<M,N,T>&, const Vector<N,T>&);
template<unsigned int M, unsigned int N, class T> Vector<N,T>   operator*(const Vector<M,T>&,   const Matrix<M,N,T>&);

#include <softlib/General/Matrix.tcc>

#endif/*_MATRIX_H*/
