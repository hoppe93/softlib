#ifndef _MATRIX_H
#define _MATRIX_H

//#define MATRIX_BOUNDS_CHECK

#include <softlib/config.h>

template<unsigned int M, unsigned int N, class T = slibreal_t>
class Matrix {
    protected:
        const unsigned int
            m = M, n = N;
        T elems[M*N];

    public:
        Matrix<M,N>();
        Matrix<M,N>(const T*);
        Matrix<M,N>(const T**);
        Matrix<M,N>(const Matrix<M,N>&);
        Matrix<M,N>(const Vector<M,T>&, const Vector<N,T>&);
        ~Matrix<M,N>();

        T& operator[](const unsigned int);
        const T& operator[](const unsigned int) const;
        T& operator[](const unsigned int, const unsigned int);
        const T& operator[](const unsigned int, const unsigned int) const;

        /* Arithmetic */
        Matrix<M,N>& operator=(const Matrix<M,N>&);
        Matrix<M,N>& operator+=(const Matrix<M,N>&);
        Matrix<M,N>& operator+=(const T);
        Matrix<M,N>& operator-=(const Matrix<M,N>&);
        Matrix<M,N>& operator-=(const T);
        Matrix<M,N>& operator*=(const Matrix<N,M>&);
        Matrix<M,N>& operator*=(const T);
        Matrix<M,N>& operator-();

        unsigned int Rows() const;
        unsigned int Cols() const;
};

/**********************
 * OPERATOR OVERLOADS *
 **********************/
template<unsigned int M, unsigned int N> Matrix<M,N> operator+(const Matrix<M,N>&, const Matrix<M,N>&);
template<unsigned int M, unsigned int N> Matrix<M,N> operator+(const Matrix<M,N>&, const T);
template<unsigned int M, unsigned int N> Matrix<M,N> operator+(const T,            const Matrix<M,N>&);
template<unsigned int M, unsigned int N> Matrix<M,N> operator-(const Matrix<M,N>&, const Matrix<M,N>&);
template<unsigned int M, unsigned int N> Matrix<M,N> operator-(const Matrix<M,N>&, const T);
template<unsigned int M, unsigned int N> Matrix<M,N> operator-(const T,            const Matrix<M,N>&);
template<unsigned int M, unsigned int N, unsigned int K> Matrix<M,K> operator*(const Matrix<M,N>&, const Matrix<N,K>&);
template<unsigned int M, unsigned int N> Matrix<M,N> operator*(const Matrix<M,N>&, const T);
template<unsigned int M, unsigned int N> Matrix<M,N> operator*(const T,            const Matrix<M,N>&);
template<unsigned int M, unsigned int N> Vector<M,T>   operator*(const Matrix<M,N>&, const Vector<N,T>&);
template<unsigned int M, unsigned int N> Vector<N,T>   operator*(const Vector<M,T>&,   const Matrix<M,N>&);

#include <softlib/General/Matrix.tcc>

#endif/*_MATRIX_H*/
