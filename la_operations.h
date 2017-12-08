#ifndef LA_OPERATIONS_H
#define LA_OPERATIONS_H

#include <la_objects.h>
#include <la_wrapper.h>

namespace la_operations
{

/// here we have the apply function to resolve expression and call contraction
void apply(la_objects::LAMatrix<double>& _dest, const la_operations::BinaryExpression<la_objects::LAMatrix<double>, la_objects::LAMatrix<double> >& _src_expr)
{
    la_operations::contract(_src_expr.larg, _src_expr.rarg, _dest);
}

template <typename T>
void copy_data(const la_objects::LABaseObject<T>& _src, la_objects::LABaseObject<T>& _dest)
{
    _dest.resize(_src.n_rows(), _src.n_cols());
    blas_wrapper::copy(_src.n_rows() * _src.n_cols(),
                       _src.get_data_ptr(),
                       1,
                       _dest.get_data_ptr(),
                       1);
}

template <typename T>
void scale(const T& _scale, la_objects::LABaseObject<T>& _dest)
{
    int elem_dist = 1;
    int n = _dest.n_rows() * _dest.n_cols();

    blas_wrapper::scal(n, _scale, _dest.get_data_ptr(), elem_dist);
}

template <typename T>
void contract(const la_objects::LAMatrix<T>& _larg, const la_objects::LAMatrix<T>& _rarg, la_objects::LAMatrix<T>& _dest)
{
    _dest.resize(_larg.n_rows(), _rarg.n_cols());

    blas_wrapper::gemm('n',
                       'n',
                       _larg.n_rows(),
                       _rarg.n_cols(),
                       _larg.n_cols(),
                       1.0,
                       _larg.get_data_ptr(),
                       _larg.leading_dim(),
                       _rarg.get_data_ptr(),
                       _rarg.leading_dim(),
                       0.0,
                       _dest.get_data_ptr(),
                       _dest.leading_dim());
}

// Calculate singular-value-decomposition of $m \times n$ Matrix $A$ of type s, d, c, z with $A=U\cdot S\cdot Vt$.
// \param A : input matrix, contents are overwritten during lapack function call.
// \param U : output matrix of dimension $m \times min(m,n)$, contains left singular vectors.
// \param S : output vector of dimensionn $min(m,n)\times 1$, contains singular values of real type.
// \param Vt : output matrix of dimension $min(m,n) \times n$, contains right singular vectors.
template <typename T>
void SVD(la_objects::LAMatrix<T>& _arg, la_objects::LAMatrix<T>& _U, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _S, la_objects::LAMatrix<T>& _Vt)
{
    const LAPACK_INT m = _arg.n_rows();
    const LAPACK_INT n = _arg.n_cols();
    const LAPACK_INT min_mn = std::min(m, n);

    //Resize output-matrices:
    _U.resize(m, min_mn);
    _S.resize(min_mn, 1);
    _Vt.resize(min_mn, n);

    //Prepare call to LAPACK function:
    la_objects::LAMatrix<T> WORK(1,1);
    LAPACK_INT LWORK = -1;
    LAPACK_INT INFO = 0;

    char jobu = 'S';      //The first min(m,n) columns/rows of U/V**T (the left/right singular vectors) are returned in the array U/Vt
    char jobvt = 'S';     //
    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> RWORK(5*min_mn, 1);        // Should be NULL if T is not a complex type.

    LAPACK_INT lda = _arg.leading_dim();
    LAPACK_INT ldu = _U.leading_dim();
    LAPACK_INT ldvt = _Vt.leading_dim();

    //Get optimal workspace:
    lapack_wrapper::gesvd(&jobu,
                          &jobvt,
                          &m,
                          &n,
                          _arg.get_data_ptr(),
                          &lda, \
                          _S.get_data_ptr(),
                          _U.get_data_ptr(),
                          &ldu, \
                          _Vt.get_data_ptr(),
                          &ldvt,
                          WORK.get_data_ptr(),\
                          &LWORK,
                          RWORK.get_data_ptr(),
                          &INFO);

    //Resize workspace according to workspace query:
    if(INFO == 0)
    {
        LWORK = (LAPACK_INT)(abs(WORK(0,0)));
        WORK.resize(LWORK, 1);
    }
    else //Throw ecxeption in case of error.
    {
        if (INFO > 0)
        {
            throw std::runtime_error(std::string("SVD: ") + "error calculating workspace");
        }
        else if (INFO < 0)
        {
            throw std::runtime_error(std::string("SVD: ") + "error calculating workspace, argument " + std::to_string(-INFO) + " invalid");
        }
    }

//Call calculation routine:
    lapack_wrapper::gesvd(&jobu,
                          &jobvt,
                          &m,
                          &n,
                          _arg.get_data_ptr(),
                          &lda, \
                          _S.get_data_ptr(),
                          _U.get_data_ptr(),
                          &ldu, \
                          _Vt.get_data_ptr(),
                          &ldvt,
                          WORK.get_data_ptr(),\
                          &LWORK,
                          RWORK.get_data_ptr(),
                          &INFO);

    if (INFO > 0)
    {
        throw std::runtime_error(std::string("SVD: ") + "DBDSQR not cnoverged");
    }
    else if (INFO < 0)
    {
        throw std::runtime_error(std::string("SVD: ") + "argument " + std::to_string(-INFO) + " invalid");
    }
}

// @sect5{Eigen functions}
// @sect6{Function: Eigen}
// \param A : Input/output matrix; on entry: symmetric/hermitsh matrix. On exit: contains the eigenvectors of A.
// \param EV : Contains the real eigenvalues in ascending order
template <typename T>
void Eigen(la_objects::LAMatrix<T>& _A,
           la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _EV)
{
    if (_A.n_rows() != _A.n_cols())
    {
        throw std::runtime_error(std::string("Eigen:" ) + "n_rows != n_cols");
    }

    char jobz[] = "V";
    char uplo[] = "U";

    const LAPACK_INT n = _A.n_rows();
    _EV.resize(n,1);

    //Prepare call to LAPACK function:
    la_objects::LAMatrix<T> WORK(1,1);
    LAPACK_INT LWORK = -1;
    LAPACK_INT INFO = 0;
    LAPACK_INT lda = _A.leading_dim();
    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> RWORK(3*n-2,1);        // Should be NULL if T is not a complex type.

    //Get optimal workspace:
    lapack_wrapper::xxev(jobz,
                         uplo,
                         &n,
                         _A.get_data_ptr(),
                         &lda, \
                         _EV.get_data_ptr(),
                         WORK.get_data_ptr(),
                         &LWORK,
                         RWORK.get_data_ptr(),
                         &INFO);


    //Resize workspace according to workspace query:
    if(INFO == 0)
    {
        LWORK = LAPACK_INT(abs(WORK(0,0)));
        WORK.resize(LWORK,1);
    }
    else //Throw ecxeption in case of error.
    {
        if (INFO > 0)
        {
            throw std::runtime_error(std::string("Eigen: ") + "error calculating workspace");
        }
        else if (INFO < 0)
        {
            throw std::runtime_error(std::string("Eigen: ") + "error calculating workspace, argument " + std::to_string(-INFO) + " invalid");
        }
    }

    //Call calculation routine:
    lapack_wrapper::xxev(jobz,
                         uplo,
                         &n,
                         _A.get_data_ptr(),
                         &lda, \
                         _EV.get_data_ptr(),
                         WORK.get_data_ptr(),
                         &LWORK,
                         RWORK.get_data_ptr(),
                         &INFO);
    if (INFO > 0)
    {
        throw std::runtime_error(std::string("Eigen: ") + "not cnoverged");
    }
    else if (INFO < 0)
    {
        throw std::runtime_error(std::string("Eigen: ") + "argument " + std::to_string(-INFO) + " invalid");
    }
}

// @sect6{Function: Eigen}
// \param A : Input matrix: symmetric/hermitsh matrix. As it is const, no eigenvectors are calculated.
// \param EV : Contains the real eigenvalues in ascending order
template <typename T>
void Eigen(const la_objects::LAMatrix<T>& _A,
           la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _EV)
{
    if (_A.n_rows() != _A.n_cols())
    {
        throw std::runtime_error(std::string("Eigen:" ) + "n_rows != n_cols");
    }

    la_objects::LAMatrix<T> A_tmp(_A);

    char jobz[] = "N";
    char uplo[] = "U";

    const LAPACK_INT n = _A.n_rows();

    //Prepare call to LAPACK function:
    la_objects::LAMatrix<T> WORK(1,1);
    LAPACK_INT * LWORK = -1;
    LAPACK_INT INFO = 0;
    LAPACK_INT lda = _A.leading_dim();

    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> RWORK(3*n-2, 1);        // Should be NULL if T is not a complex type.

    //Get optimal workspace:
    xxev(jobz,
         uplo,
         &n,
         A_tmp.get_data_ptr(),
         &lda,
         _EV.get_data_ptr(),
         WORK.get_data_ptr(),
         &LWORK,
         RWORK.get_data_ptr(),
         &INFO);


    //Resize workspace according to workspace query:
    if(INFO == 0)
    {
        LWORK = LAPACK_INT(abs(WORK(0,0)));
        WORK.resize(LWORK,1);
    }
    else //Throw ecxeption in case of error.
    {
        if (INFO > 0)
        {
            throw std::runtime_error(std::string("Eigen: ") + "error calculating workspace");
        }
        else if (INFO < 0)
        {
            throw std::runtime_error(std::string("Eigen: ") + "error calculating workspace, argument " + std::to_string(-INFO) + " invalid");
        }
    }

    //Call calculation routine:
    xxev(jobz,
         uplo,
         &n,
         A_tmp.get_data_ptr(),
         &lda,
         _EV.get_data_ptr(),
         WORK.get_data_ptr(),
         &LWORK,
         RWORK.get_data_ptr(),
         &INFO);

    if (INFO > 0)
    {
        throw std::runtime_error(std::string("Eigen: ") + "not cnoverged");
    }
    else if (INFO < 0)
    {
        throw std::runtime_error(std::string("Eigen: ") + "argument " + std::to_string(-INFO) + " invalid");
    }
}

// @sect6{Function: Eigen}
// \param D : Input/Output nx1 vector: n diagonal elements of the tridiagonal matrix A. On exit: Eigenvalues in ascending order.
// \param E : Input/Output (n-1)x1 vector: (n-1) sub-diagonal elements of the tridiagonal matrix A, stored in (1, N-1). On Exit: Content get destroyed.
// \param Z : Output matrix: Eigenvectors of the matrix A
template <typename T>
void TriEigen(la_objects::LAMatrix<T>& _D,
              la_objects::LAMatrix<T>& _E,
              la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _Z)
{
    if (_D.n_rows() != (_E.n_rows() - 1))
    {
        throw std::runtime_error(std::string("Tri-Eigen:" ) + "main diagonal and subdiagonal have mismatching value-count");
    }

    char jobz = 'V';
    const LAPACK_INT n = _D.n_elements_active;
    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> WORK(n==1?1:(2*n)-2,1);
    LAPACK_INT INFO = 0;

    _Z.reinit(n, n);
    LAPACK_INT lda = _Z.leading_dim();

    xtev(&jobz,
         &n,
         _D.get_data_ptr(),
         _E.get_data_ptr()(),
         _Z.get_data_ptr()(),
         &lda,
         WORK.get_data_ptr()(),
         &INFO);

    if (INFO > 0)
    {
        throw std::runtime_error(std::string("Tri-Eigen: ") + "not cnoverged");
    }
    else if (INFO < 0)
    {
        throw std::runtime_error(std::string("Tri-Eigen: ") + "argument " + std::to_string(-INFO) + " invalid");
    }
}

// @sect6{Function: Eigen}
// \param D : Input/Output nx1 vector: n diagonal elements of the tridiagonal matrix A. On exit: Eigenvalues in ascending order.
// \param E : Input/Output (n-1)x1 vector: (n-1) sub-diagonal elements of the tridiagonal matrix A, stored in (1, N-1). On Exit: Content get destroyed.
template <typename T>
void TriEigen(la_objects::LAMatrix<T>& _D,
              la_objects::LAMatrix<T>& _E)
{
    if (_D.n_rows() != (_E.n_rows() - 1))
    {
        throw std::runtime_error(std::string("Tri-Eigen:" ) + "main diagonal and subdiagonal have mismatching value-count");
    }

    char jobz[] = "N";

    const LAPACK_INT n = _D.n_elements_active;
    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> WORK(1,1);
    la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType> Z(1,1);
    LAPACK_INT INFO = 0;
    LAPACK_INT lda = Z.leading_dim();

    xtev(&jobz,
         &n,
         _D.get_data_ptr(),
         _E.get_data_ptr()(),
         Z.get_data_ptr()(),
         &lda,
         WORK.get_data_ptr()(),
         &INFO);

    if (INFO > 0)
    {
        throw std::runtime_error(std::string("Tri-Eigen: ") + "not cnoverged");
    }
    else if (INFO < 0)
    {
        throw std::runtime_error(std::string("Tri-Eigen: ") + "argument " + std::to_string(-INFO) + " invalid");
    }
}

} // END NAMESPACE la_operations

#endif // LA_OPERATIONS_H

