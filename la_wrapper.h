#ifndef LA_WRAPPER
#define LA_WRAPPER

// #include <gsl/gsl_cblas.h>
#include <cblas.h>
#include <complex>
#include <lapack_header.h>

namespace blas_wrapper
{

//! @param n : Anzahl Elemente.
//! @param x : Quellvektor x.
//! @param incx : Speicher Abstand zwischen Elemente in Vector x.
//! @param y : Zielvektor (y= alpha*x+y).
//! @param incy: Speicher Abstand zwischen Elemente in Vector y.

static void
copy(int n, const float *x, int incx, float *y, int incy)
{
    cblas_scopy(n, x, incx, y, incy);
}

static void
copy(int n, const double *x, int incx, double *y, int incy)
{
    cblas_dcopy(n, x, incx, y, incy);
}

static void
copy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy)
{
    cblas_zcopy(n, reinterpret_cast<const double*>(x), incx, reinterpret_cast<double*>(y), incy);
}

static void
copy(int n, const std::complex<float> *x, int incx, std::complex<float> *y, int incy)
{
    cblas_ccopy(n, reinterpret_cast<const float*>(x), incx, reinterpret_cast<float*>(y), incy);
}


static void
copy(int n, const int *x, int incx, int *y, int incy)
{
    for(int i=0; i<n; i++)
    {
        y[i*incy]=x[i*incx];
    }
}


//! @param n : Anzahl Elemente.
//! @param alpha: Skalar.
//! @param x : Quellvektor x (auch output).
//! @param incx : Speicher Abstand zwischen Elemente in Vector x.


static void
scal (int n, float alpha, float *x, int incx)
{
    cblas_sscal (n, alpha, x, incx);
}

static void
scal (int n, double alpha, double *x, int incx)
{
    cblas_dscal (n, alpha, x, incx);
}

static void
scal (int n, std::complex<float> alpha, std::complex<float> *x, int incx)
{
   cblas_cscal(n, reinterpret_cast<const float*>(&alpha), reinterpret_cast<float*>(x), incx);
}

static void
scal (int n, std::complex<double>  alpha, std::complex<double> *x, int incx)
{
   cblas_zscal(n, reinterpret_cast<const double*>(&alpha), reinterpret_cast<double*>(x), incx);
}


//! @param trans : gibt an ob A transponiert ist oder nicht. Sei trans = 'N' oder 'n' so ist op(A)= A, sei trans = 'T', 't','C' oder 'c' so ist op(A)= trans(A)
//! @param m : Anzahl Zeilen in Matrix A.
//! @param n : Anzahl Spalten in Matrix A.
//! @param alpha: Skalar fuer A.
//! @param A : Matrix A
//! @param lda : leading dimension von A.
//! @param x : Vektor mit der laenge von mindestens (1+(n-1)*abs(incx)) falls trans = 'N' oder 'n', sonst mindestens der laenge (1+(m-1)*abs(incx)).
//! @param incx : Speicher Abstand zwischen Elemente in Vector x.
//! @param beta : Skalar fuer Vektor y.
//! @param y : Vektor mit der laenge von mindestens (1+(n-1)*abs(incy)) falls trans = 'N' oder 'n', sonst mindestens der laenge (1+(m-1)*abs(incy)).
//! @param incy: Speicher Abstand zwischen Elemente in Vector y.

static void gemv (char trans, int m, int n, float alpha,
                  const float * const A, int lda,
                  const float * const x,  int incx, float beta,
                  float *y, int incy)
{
    CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
    cblas_sgemv (CblasColMajor, tr, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

static void gemv (char trans, int m, int n, double alpha,
                  const double * const A, int lda,
                  const double * const x,  int incx, double beta,
                  double *y, int incy)
{
    CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
    cblas_dgemv (CblasColMajor, tr, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

static void gemv (char trans, int m, int n, std::complex<float> & alpha,
                  const std::complex<float> * const A, int lda,
                  const std::complex<float> * const x,  int incx, std::complex<float> & beta,
                  std::complex<float> *y, int incy)
{
    CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
    cblas_cgemv (CblasColMajor, tr, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
}

static void gemv (char trans, int m, int n, std::complex<double> & alpha,
                  const std::complex<double> * const A, int lda,
                  const std::complex<double> * const x,  int incx, std::complex<double>  & beta,
                  std::complex<double> *y, int incy)
{
    CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
    cblas_zgemv (CblasColMajor, tr, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
}


//! @param transa : gibt an ob A transponiert ist oder nicht. Sei transa = 'N' oder 'n' so ist op(A)= A, sei transa = 'T' oder 't' so ist op(A)= trans(A), sei transa = 'C' oder 'c' so ist op(A)=adjoint(A)
//! @param transb : gibt an ob B transponiert ist oder nicht. Sei transb = 'N' oder 'n' so ist op(B)= A, sei transb = 'T' oder 't' so ist op(B)= trans(B), sei transb = 'C' oder 'c' so ist op(B)=adjoint(B)
//! @param m : Anzahl Zeilen in Matrix A und Matrix C.
//! @param n : Anzahl Spalten in Matrix B und Matrix C.
//! @param k : Anzahl Spalten in Matrix A und Zeilen in Matrix B.
//! @param alpha: Skalar fuer op(A)*op(B).
//! @param A : Matrix A
//! @param lda : leading dimension von A.
//! @param B : Matrix B.
//! @param ldb : leading dimension von B.
//! @param beta : Skalar fuer Matrix C.
//! @param C : Matrix C.
//! @param ldc : leading dimension von C.
static void gemm(char transa, char transb, int m, int n, int k, float alpha,
                 const float * const A, int lda, const float * const B, int ldb,
                 float beta, float * C, int ldc)
{
    CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasTrans : CblasNoTrans );
    CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasTrans : CblasNoTrans );


    cblas_sgemm(CblasColMajor,
                tr_a, tr_b,
                m, n, k,
                alpha,
                A, lda,
                B, ldb,
                beta,
                C, ldc);
}


static void gemm(char transa, char transb, int m, int n, int k, double alpha,
                 const double * const A, int lda, const double * const B, int ldb,
                 double beta, double * C, int ldc)
{
    CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasTrans : CblasNoTrans );
    CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasTrans : CblasNoTrans );

    cblas_dgemm(CblasColMajor, tr_a, tr_b, m, n, k, alpha,
                A, lda, B, ldb,
                beta, C, ldc);
}

static void gemm(char transa, char transb, int m, int n, int k, const std::complex<float> alpha,
                 const std::complex<float> * const A, int lda, const std::complex<float> * const B, int ldb,
                 const std::complex<float> beta, std::complex<float> * C, int ldc)
{
    CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasConjTrans : CblasNoTrans );
    CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasConjTrans : CblasNoTrans );

    cblas_cgemm(CblasColMajor,
                tr_a, tr_b,
                m, n, k,
                reinterpret_cast<const float*>(&alpha),
                reinterpret_cast<const float*>(A), lda,
                reinterpret_cast<const float*>(B), ldb,
                reinterpret_cast<const float*>(&beta),
                reinterpret_cast<float*>(C), ldc);
}


static void gemm(char transa, char transb, int m, int n, int k, const std::complex<double> alpha,
                 const std::complex<double> * const A, int lda, const std::complex<double> * const B, int ldb,
                 const std::complex<double> beta, std::complex<double> * C, int ldc)
{
    CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasConjTrans : CblasNoTrans );
    CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasConjTrans : CblasNoTrans );

    cblas_zgemm(CblasColMajor,
                tr_a, tr_b,
                m, n, k,
                reinterpret_cast<const double*>(&alpha),
                reinterpret_cast<const double*>(A), lda,
                reinterpret_cast<const double*>(B), ldb,
                reinterpret_cast<const double*>(&beta),
                reinterpret_cast<double*>(C), ldc);
}


} // END NAMESPACE blas_wrapper

namespace lapack_wrapper
{

// @sect5{Lapack function wrappers}
// @sect6{Wrapper: SVD lapack-wrapper}
// Wraps the general matrix svd lapack function for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d8/d49/sgesvd_8f.html
inline static void gesvd(char * JOBU, char * JOBVT, const LAPACK_INT * M, const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, float * S, float * U, const LAPACK_INT * LDU, float * VT, const LAPACK_INT * LDVT, float * WORK, const LAPACK_INT * LWORK, float * /*RWORK*/, LAPACK_INT * INFO)
{
    sgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
}

inline static void gesvd(char * JOBU, char * JOBVT, const LAPACK_INT * M, const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, double * S, double * U, const LAPACK_INT * LDU, double * VT, const LAPACK_INT * LDVT, double * WORK, const LAPACK_INT * LWORK, double * /*RWORK*/, LAPACK_INT * INFO)
{
    dgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
}

inline static void gesvd(char * JOBU, char * JOBVT, const LAPACK_INT * M, const LAPACK_INT * N, std::complex<float> * A, const LAPACK_INT * LDA, float * S, std::complex<float> * U, const LAPACK_INT  * LDU, std::complex<float> * VT, const LAPACK_INT * LDVT, std::complex<float> * WORK, const LAPACK_INT * LWORK, float * RWORK, LAPACK_INT * INFO)
{
    cgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

inline static void gesvd(char * JOBU, char * JOBVT, const LAPACK_INT * M, const LAPACK_INT * N, std::complex<double> * A, const LAPACK_INT * LDA, double * S, std::complex<double> * U, const LAPACK_INT * LDU, std::complex<double> * VT, const LAPACK_INT * LDVT, std::complex<double> * WORK, const LAPACK_INT * LWORK, double * RWORK, LAPACK_INT * INFO)
{
    zgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

// @sect6{Wrapper: SVD divide and conquer lapack-wrapper}
// Wraps the general matrix sdd lapack function for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga76f797b6a9e278ad7b21aae2b4a55d76.html#ga76f797b6a9e278ad7b21aae2b4a55d76
inline static void gesdd(char * JOBZ, const LAPACK_INT * M, const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, float * S, float * U, const LAPACK_INT * LDU, float * VT, const LAPACK_INT * LDVT, float * WORK, const LAPACK_INT * LWORK, float * /*RWORK*/, LAPACK_INT * IWORK, LAPACK_INT * INFO)
{
    return sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);
}
inline static void gesdd(char * JOBZ, const LAPACK_INT * M, const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, double * S, double * U, const LAPACK_INT * LDU, double * VT, const LAPACK_INT * LDVT, double * WORK, const LAPACK_INT * LWORK, double * /*RWORK*/, LAPACK_INT * IWORK, LAPACK_INT * INFO)
{
    return dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);
}
inline static void gesdd(char * JOBZ, const LAPACK_INT * M, const LAPACK_INT * N, std::complex<float> * A, const LAPACK_INT * LDA, float * S, std::complex<float> * U, const LAPACK_INT * LDU, std::complex<float> * VT, const LAPACK_INT * LDVT, std::complex<float> * WORK, const LAPACK_INT * LWORK, float * RWORK, LAPACK_INT * IWORK, LAPACK_INT * INFO)
{
    return cgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);
}
inline static void gesdd(char * JOBZ, const LAPACK_INT * M, const LAPACK_INT * N, std::complex<double> * A, const LAPACK_INT * LDA, double * S, std::complex<double> * U, const LAPACK_INT * LDU, std::complex<double> * VT, const LAPACK_INT * LDVT, std::complex<double> * WORK, const LAPACK_INT * LWORK, double * RWORK, LAPACK_INT * IWORK, LAPACK_INT * INFO)
{
    return zgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);
}

// @sect6{Wrapper: LUD lapack-wrapper}
// Wraps the general matrix LU decomposition lapack function for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/de/de2/sgetrf_8f.html
inline static void getrf(const LAPACK_INT * M, const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, LAPACK_INT * IPIV, LAPACK_INT * INFO)
{
    sgetrf_(M, N, A, LDA, IPIV, INFO);
}
inline static void getrf(const LAPACK_INT * M, const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, LAPACK_INT * IPIV, LAPACK_INT * INFO)
{
    dgetrf_(M, N, A, LDA, IPIV, INFO);
}
inline static void getrf(const LAPACK_INT * M, const LAPACK_INT * N, std::complex<float> * A, const LAPACK_INT * LDA, LAPACK_INT * IPIV, LAPACK_INT * INFO)
{
    cgetrf_(M, N, A, LDA, IPIV, INFO);
}
inline static void getrf(const LAPACK_INT * M, const LAPACK_INT * N, std::complex<double> * A, const LAPACK_INT * LDA, LAPACK_INT * IPIV, LAPACK_INT * INFO)
{
    zgetrf_(M, N, A, LDA, IPIV, INFO);
}

// @sect6{Wrapper: LU-Inverse lapack-wrapper}
// Wraps the general matrix inversion lapack function for the two types s, d.
 // Lapack doc: http://www.netlib.org/lapack/explore-html/de/de2/sgetri_8f.html
inline static void getri (const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, const LAPACK_INT * IPIV, float * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    sgetri_(N, A, LDA, IPIV, WORK, LWORK, INFO);
}
inline static void getri (const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, const LAPACK_INT * IPIV, double * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    dgetri_(N, A, LDA, IPIV, WORK, LWORK, INFO);
}

// @sect6{Wrapper: QRF lapack-wrapper}
// Wraps the general matrix QR factorization lapack function for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/df/d97/sgeqrf_8f.html
inline static void geqrf(const LAPACK_INT * M, const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, float * TAU, float * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    sgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}
inline static void geqrf(const LAPACK_INT * M, const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, double * TAU, double * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    dgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}
inline static void geqrf(const LAPACK_INT * M, const LAPACK_INT * N, std::complex<float> * A, const LAPACK_INT * LDA, std::complex<float> * TAU, std::complex<float> * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    cgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}
inline static void geqrf(const LAPACK_INT * M, const LAPACK_INT * N, std::complex<double> * A, const LAPACK_INT * LDA, std::complex<double> * TAU, std::complex<double> * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    zgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
}

// @sect6{Wrapper: MQR lapack-wrapper}
// Wraps the lapack function of the product of elementary reflectors from the QR factorization for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d0/d98/sormqr_8f.html
inline static void xxmqr(char * SIDE, char * TRANS, const LAPACK_INT * M, const LAPACK_INT * N, const LAPACK_INT * K, float * A, const LAPACK_INT * LDA, float * TAU, float * C, const LAPACK_INT * LDC, float * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
}
inline static void xxmqr(char * SIDE, char * TRANS, const LAPACK_INT * M, const LAPACK_INT * N, const LAPACK_INT * K, double * A, const LAPACK_INT * LDA, double * TAU, double * C, const LAPACK_INT * LDC, double * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
}
inline static void xxmqr(char * SIDE, char * TRANS, const LAPACK_INT * M, const LAPACK_INT * N, const LAPACK_INT * K, std::complex<float> * A, const LAPACK_INT * LDA, std::complex<float> * TAU, std::complex<float> * C, const LAPACK_INT * LDC, std::complex<float> * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    cunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
}
inline static void xxmqr(char * SIDE, char * TRANS, const LAPACK_INT * M, const LAPACK_INT * N, const LAPACK_INT * K, std::complex<double> * A, const LAPACK_INT * LDA, std::complex<double> * TAU, std::complex<double> * C, const LAPACK_INT * LDC, std::complex<double> * WORK, const LAPACK_INT * LWORK, LAPACK_INT * INFO)
{
    zunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
}

// @sect6{Wrapper: EV lapack-wrapper}
// Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the four types s, d, c and z.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html
inline static void xxev(char * JOBZ, char * UPLO, const LAPACK_INT * N, float * A, const LAPACK_INT * LDA, float * W, float * WORK, const LAPACK_INT * LWORK, float * /*RWORK*/, LAPACK_INT * INFO)
{
    return ssyev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
}
inline static void xxev(char * JOBZ, char * UPLO, const LAPACK_INT * N, double * A, const LAPACK_INT * LDA, double * W, double * WORK, const LAPACK_INT * LWORK, double * /*RWORK*/, LAPACK_INT * INFO)
{
    return dsyev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
}
inline static void xxev(char * JOBZ, char * UPLO, const LAPACK_INT * N, std::complex<float> * A, const LAPACK_INT * LDA, float * W, std::complex<float> * WORK, const LAPACK_INT * LWORK, float * RWORK, LAPACK_INT * INFO)
{
    return cheev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO);
}
inline static void xxev(char * JOBZ, char * UPLO, const LAPACK_INT * N, std::complex<double> * A, const LAPACK_INT * LDA, double * W, std::complex<double> * WORK, const LAPACK_INT * LWORK, double * RWORK, LAPACK_INT * INFO)
{
    return zheev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO);
}

// @sect6{Wrapper: EV lapack-wrapper}
// Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the two real types s and d.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html
inline static void xtev(char * JOBZ, const LAPACK_INT * N, float * D, float * E, float * Z, const LAPACK_INT * LDZ, float * WORK, LAPACK_INT * INFO)
{
    return sstev_(JOBZ, N, D, E, Z, LDZ, WORK, INFO);
}
inline static void xtev(char * JOBZ, const LAPACK_INT * N, double * D, double * E, double * Z, const LAPACK_INT * LDZ, double * WORK, LAPACK_INT * INFO)
{
    dstev_(JOBZ, N, D, E, Z, LDZ, WORK, INFO);
}

// @sect6{Wrapper: Matrix norm lapack-wrapper}
// Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the two real types s and d.
// Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html
inline static float lange(char * NORM, const LAPACK_INT * M, const LAPACK_INT * N, const float * A, const LAPACK_INT * LDA, float * WORK)
{
    return slange_( NORM, M, N, A, LDA, WORK );
}
inline static double lange(char * NORM, const LAPACK_INT * M, const LAPACK_INT * N, const double * A, const LAPACK_INT * LDA, double * WORK)
{
    return dlange_( NORM, M, N, A, LDA, WORK );
}
inline static float lange(char * NORM, const LAPACK_INT * M, const LAPACK_INT * N, const std::complex<float> * A, const LAPACK_INT * LDA, float * WORK)
{
    return clange_( NORM, M, N, A, LDA, WORK );
}
inline static double lange(char * NORM, const LAPACK_INT * M, const LAPACK_INT * N, const std::complex<double> * A, const LAPACK_INT * LDA, double * WORK)
{
    return zlange_( NORM, M, N, A, LDA, WORK );
}

} // END NAMESPACE LAPACK_WRAPPER

#endif // LA_WRAPPER

