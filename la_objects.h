#ifndef LA_OBJECTS_H
#define LA_OBJECTS_H

#include <la_base_obj.h>

namespace la_objects
{

template <typename T> class LAMatrix;

} // END NAMESPACE la_objects

namespace la_operations
{

template <typename Type>
void contract(const la_objects::LAMatrix<Type>& _larg, const la_objects::LAMatrix<Type>& _rarg, la_objects::LAMatrix<Type>& _dest);

template <typename Type>
void svd(la_objects::LAMatrix<Type>& _arg, la_objects::LAMatrix<Type>& _out);

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
    friend void la_operations::svd(la_objects::LAMatrix<Type>& _arg, la_objects::LAMatrix<Type>& _U, la_objects::LAMatrix<Type>& _S, la_objects::LAMatrix<Type>& _Vt);

//    template <typename Type>
//    friend void la_operations::kronicker_product(la_objects::LAMatrix<Type>& _larg, la_objects::LAMatrix<Type>& _rarg, la_objects::LAMatrix<Type>& _out);

};

} // END NAMESPACE la_objects

template <typename T>
la_objects::LAMatrix<T> operator*(const la_objects::LAMatrix<T>& _larg, const la_objects::LAMatrix<T>& _rarg);
#endif // LA_OBJECTS_H
