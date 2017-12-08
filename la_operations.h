#ifndef LA_OPERATIONS
#define LA_OPERATIONS

#include <la_wrapper.h>
#include <la_base_obj.h>
#include <la_objects.h>

namespace la_operations
{

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

template <typename T>
void svd(la_objects::LAMatrix<T>& _arg, la_objects::LAMatrix<T>& _U, la_objects::LAMatrix<T>& _S, la_objects::LAMatrix<T>& _VT)
{

    //set the JOB* parameters
    char JOBU = 'A'; //all M columns of U are returned in array U
    char JOBVT = 'A'; //all N rows of V**T are returned in the array VT

    //dimensions in LAPACK_INT type
    LAPACK_INT M = _arg.n_rows();
    LAPACK_INT N = _arg.n_cols();
    LAPACK_INT MIN_MN = std::min(M, N);

    //resize and create arrays used in computation
    _U.resize(M, MIN_MN);
    _S.resize(1, MIN_MN);
    _VT.resize(MIN_MN, N);
    la_objects::LAMatrix<T> WORK(5*MIN_MN, 1);

    LAPACK_INT LDA = _arg.leading_dim();
    LAPACK_INT INFO = 0;
    LAPACK_INT LDU = _U.leading_dim();
    LAPACK_INT LDVT = _VT.leading_dim();
    LAPACK_INT LWORK = WORK.leading_dim();


    lapack_wrapper::svd(&JOBU,
                      &JOBVT,
                      &M,
                      &N,
                      _arg.get_data_ptr(),
                      &LDA,
                      _S.get_data_ptr(),
                      _U.get_data_ptr(),
                      &LDU,
                      _VT.get_data_ptr(),
                      &LDVT,
                      WORK.get_data_ptr(),
                      &LWORK,
                      &INFO);
}

//template <typename T>
//void kronicker_product(la_objects::LAMatrix<T>& _larg, la_objects::LAMatrix<T>& _rarg, la_objects::LAMatrix<T>& _out)
//{
//    lapack_wrapper::kronicker_product(_larg, _rarg, _out);
//}

} // END NAMESPACE la_operations

#endif // LA_OPERATIONS

