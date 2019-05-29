/**
 * Implementation of the SOFTlib Matrix type.
 */

#include <cmath>
#include <softlib/SOFTLibException.h>
#include <softlib/Matrix.h>

/**
 * Constructor
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::Matrix() {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] = 0.0;
}

template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::Matrix(const T *mat) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] = mat[i];
}

template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::Matrix(const T **mat) {
    for (unsigned int i = 0; i < M; i++)
        for (unsigned int j = 0; j < N; j++)
            elems[i*N + j] = mat[i][j];
}

template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::Matrix(const Matrix<M,N,T>& mat) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] = mat(i);
}

/**
 * Destructor.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::~Matrix() {}

/**
 * Construct the matrix from the two vectors
 * a and b according to
 *
 *   A = a b^T
 * 
 * where A is the matrix and T denotes transpose.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>::Matrix(const Vector<M,T>& a, const Vector<N,T>& b) {
    for (unsigned int i = 0; i < M; i++)
        for (unsigned int j = 0; j < N; j++)
            elems[i*N+j] = a[i] * b[j];
}

/**
 * Access a given element of the matrix.
 */
template<unsigned int M, unsigned int N, class T>
T& Matrix<M,N,T>::operator()(const unsigned int i) {
#   ifdef MATRIX_BOUNDS_CHECK
    if (i < 0 || i >= M*N)
        throw SOFTLibException("Matrix linear index is out of bounds: %d.\n", i);
#   endif
    return elems[i];
}
template<unsigned int M, unsigned int N, class T>
const T& Matrix<M,N,T>::operator()(const unsigned int i) const {
#   ifdef MATRIX_BOUNDS_CHECK
    if (i < 0 || i >= M*N)
        throw SOFTLibException("Matrix linear index is out of bounds: %d.\n", i);
#   endif
    return elems[i];
}

template<unsigned int M, unsigned int N, class T>
T& Matrix<M,N,T>::operator()(const unsigned int i, const unsigned int j) {
#   ifdef MATRIX_BOUNDS_CHECK
    if (i < 0 || i >= M)
        throw SOFTLibException("Matrix row index is out of bounds: %d.\n", i);
    else if (j < 0 || j>= N)
        throw SOFTLibException("Matrix column index is out of bounds: %d.\n", j);
#   endif

    return elems[i*N + j];
}

template<unsigned int M, unsigned int N, class T>
const T& Matrix<M,N,T>::operator()(const unsigned int i, const unsigned int j) const {
#   ifdef MATRIX_BOUNDS_CHECK
    if (i < 0 || i >= M)
        throw SOFTLibException("Matrix row index is out of bounds: %d.\n", i);
    else if (j < 0 || j>= N)
        throw SOFTLibException("Matrix column index is out of bounds: %d.\n", j);
#   endif

    return elems[i*N + j];
}

/**
 * Returns the number of rows in this matrix
 * (i.e. elements along first dimension).
 */
template<unsigned int M, unsigned int N, class T>
unsigned int Matrix<M,N,T>::Rows() const { return M; }

/**
 * Returns the number of columns in this matrix
 * (i.e. elements along second dimension).
 */
template<unsigned int M, unsigned int N, class T>
unsigned int Matrix<M,N,T>::Cols() const { return N; }

/*************************
 * ARITHMETIC OPERATIONS *
 *************************/
/**
 * Copy the given matrix to this matrix.
 *
 * mat: Matrix to copy the elements of.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator=(const Matrix<M,N,T>& mat) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] = mat[i];

    return *this;
}

/**
 * Add the given matrix to this matrix.
 *
 * mat: Matrix to add.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator+=(const Matrix<M,N,T>& mat) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] += mat[i];

    return *this;
}

/**
 * Add the given scalar to every element of
 * this matrix.
 *
 * s: Scalar to add to every element.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator+=(const T s) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] += s;

    return *this;
}

/**
 * Subtract the given matrix from this matrix.
 *
 * mat: Matrix to subtract.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator-=(const Matrix<M,N,T>& mat) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] -= mat(i);
    
    return *this;
}

/**
 * Subtract the given scalar from every
 * element of this matrix.
 *
 * s: Scalar to subtract from every element of
 *    this matrix.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator-=(const T s) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] -= s;

    return *this;
}

/**
 * Multiply the given matrix with this matrix
 * (from the rignt).
 *
 * mat: Matrix to multiply with.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator*=(const Matrix<N,M,T>& mat) {
    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < M; j++) {
            T sum = 0;
            for (unsigned k = 0; k < N; k++)
                sum += elems[i*N+k] * mat(k,j);

            elems[i*N+j] = sum;
        }
    }

    return *this;
}

/**
 * Multiply this matrix with the given scalar.
 *
 * s: Scalar to multiply every element of this matrix with.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator*=(const T s) {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] *= s;

    return *this;
}

/**
 * Negate all elements of this matrix.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T>& Matrix<M,N,T>::operator-() {
    for (unsigned int i = 0; i < M*N; i++)
        elems[i] = -elems[i];

    return *this;
}


/********************************
 * RELATED NON-MEMBER FUNCTIONS *
 ********************************/
/**
 * Add the two matrices together.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator+(const Matrix<M,N,T>& lhs, const Matrix<M,N,T>& rhs) {
    Matrix<M,N,T> C(lhs);
    C += rhs;
    return C;
}

/**
 * Add a scalar and a matrix.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator+(const Matrix<M,N,T>& lhs, const T rhs) {
    Matrix<M,N,T> C(lhs);
    C += rhs;
    return C;
}
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator+(const T lhs, const Matrix<M,N,T>& rhs) {
    Matrix<M,N,T> C(rhs);
    C += lhs;
    return C;
}

/**
 * Subtract the second matrix from the first.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator-(const Matrix<M,N,T>& lhs, const Matrix<M,N,T>& rhs) {
    Matrix<M,N,T> C(lhs);
    C -= rhs;
    return C;
}

/**
 * Subtract a scalar from the matrix.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator-(const Matrix<M,N,T>& lhs, const T rhs) {
    Matrix<M,N,T> C(lhs);
    C -= rhs;
    return C;
}
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator-(const T lhs, const Matrix<M,N,T>& rhs) {
    Matrix<M,N,T> C(rhs);
    C -= lhs;
    return C;
}

/**
 * Multiply the two matrices together.
 */
template<unsigned int M, unsigned int N, unsigned int K, class T>
Matrix<M,K,T> operator*(const Matrix<M,N,T>& l, const Matrix<N,K,T>& r) {
    Matrix<M,K,T> C;

    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < K; j++) {
            T sum = 0;
            for (unsigned int k = 0; k < N; k++)
                sum += l(i,k) * r(k,j);

            C(i,j) = sum;
        }
    }

    return C;
}

/**
 * Multiply the given matrix with the given scalar.
 */
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator*(const Matrix<M,N,T>& l, const T r) {
    Matrix<M,N,T> C(l);
    C += r;
    return C;
}
template<unsigned int M, unsigned int N, class T>
Matrix<M,N,T> operator*(const T l, const Matrix<M,N,T>& r) {
    Matrix<M,N,T> C(r);
    C += l;
    return C;
}

/**
 * Multiply the given vector with the matrix
 * (from the right).
 */
template<unsigned int M, unsigned int N, class T>
Vector<M,T> operator*(const Matrix<M,N,T>& l, const Vector<N,T>& r) {
    Vector<M,T> vec;

    for (unsigned int i = 0; i < M; i++) {
        T sum = 0;
        for (unsigned int j = 0; j < N; j++)
            sum += l(i,j) * r[j];

        vec[i] = sum;
    }

    return vec;
}

/**
 * Multiply the given vector with the matrix
 * (from the left).
 */
template<unsigned int M, unsigned int N, class T>
Vector<N,T> operator*(const Vector<M,T>& l, const Matrix<M,N,T>& r) {
    Vector<N,T> vec;

    for (unsigned int i = 0; i < N; i++) {
        T sum = 0;
        for (unsigned int j = 0; j < M; j++)
            sum += l(j,i) * r[j];

        vec[i] = sum;
    }

    return vec;
}

