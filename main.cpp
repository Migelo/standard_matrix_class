#include <iomanip>
#include <la_objects.h>
#include <la_operations.h>
#include <la_objects.cpp>

int main(int /*argc*/, char **/*argv[]*/)
{
    la_objects::LAMatrix<double> A(2, 2);
    la_objects::LAMatrix<double> B(2, 2);

    la_objects::LAMatrix<double> C;
    la_objects::LAMatrix<double> D(2, 1);

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
    std::cout << "A:" << std::endl << A << std::endl;
    std::cout << "B:" << std::endl << B << std::endl;

    A *= 2.0;

    std::cout << "2*A:" << std::endl; A.print(std::cout); std::cout << std::endl;

    C = A * B;

    std::cout << "2*A*B:" << std::endl << C << std::endl;

    la_operations::evd(A, D);

    std::cout << A;

    return 0;
}

