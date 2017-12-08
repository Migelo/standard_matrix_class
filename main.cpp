#include <iomanip>
#include <la_objects.h>
#include <la_operations.h>
#include <la_objects.cpp>
#include <chrono>
#include <fstream>
#include <math.h>

la_objects::LAMatrix<double> kroniker(la_objects::LAMatrix<double> &sigma_x, la_objects::LAMatrix<double> &sigma_z)
{
    la_objects::LAMatrix<double> out(sigma_x.n_rows() * sigma_z.n_rows(), sigma_x.n_cols() * sigma_z.n_cols());

    for (unsigned int r = 0; r < sigma_x.n_rows(); r++)
    {
        for (unsigned int c = 0; c < sigma_x.n_cols(); c++)
        {
            for (unsigned int i = 0; i < sigma_z.n_rows(); i++)
            {
                for (unsigned int j = 0; j < sigma_z.n_cols(); j++)
                {
//                    std::cout << r * sigma_x.n_rows() + i << std::endl;
//                    std::cout << c * sigma_x.n_cols() + j << std::endl;
//                    std::cout << " WORK(" << i << ", " << j << ") = " << sigma_z(i, j) << std::endl;
                    out(r * sigma_x.n_rows() + i, c * sigma_x.n_cols() + j) = sigma_x(r, c) * sigma_z(i, j);
                }
            }
        }
    }
    return out;
}
template <typename T>
void sigma_i(la_objects::LAMatrix<T> &A,la_objects::LAMatrix<T> &out, int i, int L)
{
    la_objects::LAMatrix<T> id(2, 2);
    id(0, 0) = 1;
    id(0, 1) = 0;
    id(1, 0) = 0;
    id(1, 1) = 1;

    out = id;
    if (i + 1 < L)
    {
        for (int j = 0; j < L - i - 2; j++)
        {
            out = kroniker(id, out);
        }
    }
    out = kroniker(A, out);
    if (i > 0)
    {
        for (unsigned int j = 0; j < i; j++)
        {
            out = kroniker(id, out);
        }

    }
}

template <typename T>
void sum(la_objects::LAMatrix<T> &A, la_objects::LAMatrix<T> &B, la_objects::LAMatrix<T> &out)
{
    out.resize(A.n_rows(), A.n_cols());
    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            out(r, c) = A(r, c) + B(r, c);
        }
    }
}

int main(int /*argc*/, char **/*argv[]*/)
{
    int n = 2;
    int i = 5;
    la_objects::LAMatrix<double> sigma_x(2, 2);
    la_objects::LAMatrix<double> sigma_z(2, 2);
    la_objects::LAMatrix<double> out;
    la_objects::LAMatrix<double> id(2, 2);
    sigma_x(0, 0) = 0;
    sigma_x(0, 1) = 1;
    sigma_x(1, 0) = 1;
    sigma_x(1, 1) = 0;

    sigma_z(0, 0) = 1;
    sigma_z(0, 1) = 0;
    sigma_z(1, 0) = 0;
    sigma_z(1, 1) = -1;

    id(0, 0) = 1;
    id(0, 1) = 0;
    id(1, 0) = 0;
    id(1, 1) = 1;

    la_objects::LAMatrix<double> A(n, n);
    la_objects::LAMatrix<double> B(n, n);

    la_objects::LAMatrix<double> C;
//    la_objects::LAMatrix<double> D(2, 1);

    for (unsigned int r = 0; r < A.n_rows(); r++)
    {
        for (unsigned int c = 0; c < A.n_cols(); c++)
        {
            A(r, c) = c + r*A.leading_dim();
        }
    }

    for (unsigned int r = 0; r < B.n_rows(); r++)
    {
        for (unsigned int c = 0; c < B.n_cols(); c++)
        {
            B(r, c) = c + r*B.leading_dim();
        }
    }

    std::cout << std::setprecision(3);// << std::scientific;
//    std::cout << "A:" << std::endl << A << std::endl;
//    std::cout << "B:" << std::endl << B << std::endl;

//    A *= 2.0;

//    std::cout << "2*A:" << std::endl; A.print(std::cout); std::cout << std::endl;

//    C = A * B;

//    std::cout << "A*B:" << std::endl << C << std::endl;


//    std::cout << sigma_x << std::endl;
//    std::cout << sigma_z << std::endl;

//    std::cout << kroniker(id, id);

//    out = sigma_i(sigma_z, 2, 3) * sigma_i(sigma_z, 3, 3);


    int L = 6;
    int dim = pow(2, L);
    la_objects::LAMatrix<double> H(dim, dim);
    la_objects::LAMatrix<double> H_temp(dim, dim);
    la_objects::LAMatrix<double> H_z(dim, dim), H1, H2;


    for (int i = 0; i < L; i++)
    {
        sigma_i(sigma_z, H1, i, L);
        sigma_i(sigma_z, H2, i + 1, L);
        H_temp = H1 * H2;
        sum(H, H_temp, H);
    }
    for (int i = 0; i < L; i++)
    {
        sigma_i(sigma_x, H_temp, i, L);
        sum(H_z, H_temp, H_z);
    }
    H *= -1.0;
    H_z *= -.3;
    sum(H, H_z, H);
    std::cout << "H:" << std::endl << H << std::endl;
    la_objects::LAMatrix<double> U;
    la_objects::LAMatrix<double> S;
    la_objects::LAMatrix<double> Vt;

    la_operations::Eigen(H, U);
    std::cout << "Eigenvalues:" << std::endl;
    std::cout << U;
    return 0;
}



