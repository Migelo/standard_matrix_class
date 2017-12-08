#ifndef LA_OBJECTS_H
#define LA_OBJECTS_H

#include <la_base_obj.h>

namespace la_objects
{

template <typename T> class LAMatrix;

} // END NAMESPACE la_objects

namespace la_operations
{

template <typename T>
void SVD(la_objects::LAMatrix<T>& _arg, la_objects::LAMatrix<T>& _U, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _S, la_objects::LAMatrix<T>& _Vt);

template <typename T>
void Eigen(la_objects::LAMatrix<T>& _A, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _EV);

template <typename T>
void Eigen(const la_objects::LAMatrix<T>& _A, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _EV);

template <typename T>
void TriEigen(la_objects::LAMatrix<T>& _D, la_objects::LAMatrix<T>& _E, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<T>::RealType>& _Z);

template <typename T>
void TriEigen(la_objects::LAMatrix<T>& _D, la_objects::LAMatrix<T>& _E);

template <typename Type>
void contract(const la_objects::LAMatrix<Type>& _larg, const la_objects::LAMatrix<Type>& _rarg, la_objects::LAMatrix<Type>& _dest);

template <typename T>
using LAMatrixMatrix = BinaryExpression<la_objects::LAMatrix<T>, la_objects::LAMatrix<T>>;

/// overload *-operator to return binary expresion
template <typename T>
LAMatrixMatrix<T> operator*(const la_objects::LAMatrix<T>& _larg, const la_objects::LAMatrix<T>& _rarg)
{
    return LAMatrixMatrix<T>(_larg, _rarg);
}

/// in order to link compiler needs apply method declaration before class declaration
void apply(la_objects::LAMatrix<double>& _dest, const la_operations::BinaryExpression<la_objects::LAMatrix<double>, la_objects::LAMatrix<double> >& _src_expr);

} // END NAMESPACE la_operations

namespace la_objects
{

template <typename T>
class LAMatrix: public LABaseObject<T>
{

public:

    LAMatrix<T>(const size_t _n_rows, const size_t _n_cols)
        : LABaseObject<T>(_n_rows, _n_cols)
    {}
    LAMatrix<T>()
        : LABaseObject<T>()
    {}
    LAMatrix<T>(const LAMatrix<T>& _other)
        : LABaseObject<T>(_other)
    {}

    LAMatrix<T>& operator=(const LAMatrix<T>& _src)
    {
        la_operations::copy_data(_src, *this);
        return *this;
    }

    LAMatrix<T>& operator*=(const T& _value);

    template <typename Type>
    friend std::ostream& operator<<(std::ostream& _os, const LAMatrix<Type>& _src)
    {
        _src.print(_os);
        return _os;
    }

    template <typename Type>
    friend void la_operations::contract(const LAMatrix<Type>& _larg, const LAMatrix<Type>& _rarg, LAMatrix<Type>& _dest);

    template <typename Type>
    friend void la_operations::SVD(la_objects::LAMatrix<Type>& _arg, la_objects::LAMatrix<Type>& _U, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<Type>::RealType>& _S, la_objects::LAMatrix<Type>& _Vt);

    template <typename Type>
    friend void la_operations::Eigen(la_objects::LAMatrix<Type>& _A, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<Type>::RealType>& _EV);

    template <typename Type>
    friend void la_operations::Eigen(const la_objects::LAMatrix<Type>& _A, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<Type>::RealType>& _EV);

    template <typename Type>
    friend void la_operations::TriEigen(la_objects::LAMatrix<Type>& _D, la_objects::LAMatrix<Type>& _E, la_objects::LAMatrix<typename la_objects::NumberTypeTrait<Type>::RealType>& _Z);

    template <typename Type>
    friend void la_operations::TriEigen(la_objects::LAMatrix<Type>& _D, la_objects::LAMatrix<Type>& _E);

    template<typename X>
    LAMatrix<T>& operator=(const la_operations::Expression<X>& _expr)
    {
        la_operations::apply(*this, ~_expr);
        return *this;
    }
};

} // END NAMESPACE la_objects

#endif // LA_OBJECTS_H
