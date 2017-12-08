#ifndef LAPACK_WRAPPER_H
#define LAPACK_WRAPPER_H

#include<complex>

typedef signed int LAPACK_INT;

//Get externally defined LAPACK functions:
extern "C" {

//SVD routines:
void sgesvd_(
        char            * JOBU,
        char            * JOBVT,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        float           * A,
        const LAPACK_INT    * LDA,
        float           * S,
        float           * U,
        const LAPACK_INT    * LDU,
        float           * VT,
        const LAPACK_INT    * LDVT,
        float           * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void dgesvd_(
        char            * JOBU,
        char            * JOBVT,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        double          * A,
        const LAPACK_INT    * LDA,
        double          * S,
        double          * U,
        const LAPACK_INT    * LDU,
        double          * VT,
        const LAPACK_INT    * LDVT,
        double          * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void cgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<float>             * A,
        const LAPACK_INT                * LDA,
        float                           * S,
        std::complex<float>             * U,
        const LAPACK_INT                * LDU,
        std::complex<float>             * VT,
        const LAPACK_INT                * LDVT,
        std::complex<float>             * WORK,
        const LAPACK_INT            * LWORK,
        float                       * RWORK,
        LAPACK_INT                         * INFO);
void zgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<double>            * A,
        const LAPACK_INT                * LDA,
        double                          * S,
        std::complex<double>            * U,
        const LAPACK_INT                * LDU,
        std::complex<double>            * VT,
        const LAPACK_INT                * LDVT,
        std::complex<double>            * WORK,
        const LAPACK_INT            * LWORK,
        double                      * RWORK,
        LAPACK_INT                         * INFO);

//SVD divide and conquer routines:
void sgesdd_(
        char            * JOBZ,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        float           * A,
        const LAPACK_INT    * LDA,
        float           * S,
        float           * U,
        const LAPACK_INT    * LDU,
        float           * VT,
        const LAPACK_INT    * LDVT,
        float           * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * IWORK,
        LAPACK_INT             * INFO);
void dgesdd_(
        char            * JOBZ,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        double          * A,
        const LAPACK_INT    * LDA,
        double          * S,
        double          * U,
        const LAPACK_INT    * LDU,
        double          * VT,
        const LAPACK_INT    * LDVT,
        double          * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * IWORK,
        LAPACK_INT             * INFO);
void cgesdd_(
        char                        * JOBZ,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        float                       * S,
        std::complex<float>  * U,
        const LAPACK_INT                * LDU,
        std::complex<float>  * VT,
        const LAPACK_INT                * LDVT,
        std::complex<float>  * WORK,
        const LAPACK_INT            * LWORK,
        float                       * RWORK,
        LAPACK_INT                         * IWORK,
        LAPACK_INT                         * INFO);
void zgesdd_(
        char                        * JOBZ,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        double                      * S,
        std::complex<double> * U,
        const LAPACK_INT                * LDU,
        std::complex<double> * VT,
        const LAPACK_INT                * LDVT,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        double                      * RWORK,
        LAPACK_INT                         * IWORK,
        LAPACK_INT                         * INFO);


//--------------------------------------------
//LUD routines:
void sgetrf_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        float           * A,
        const LAPACK_INT    * LDA,
        LAPACK_INT             * IPIV,
        LAPACK_INT             * INFO);
void dgetrf_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        double          * A,
        const LAPACK_INT    * LDA,
        LAPACK_INT             * IPIV,
        LAPACK_INT             * INFO);
void cgetrf_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        LAPACK_INT                         * IPIV,
        LAPACK_INT                         * INFO);
void zgetrf_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        LAPACK_INT                         * IPIV,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//LU-inverse routines:
void sgetri_(
        const LAPACK_INT        * N,
        float               * A,
        const LAPACK_INT        * LDA,
        const LAPACK_INT           * IPIV,
        float               * WORK,
        const LAPACK_INT    * LWORK,
        LAPACK_INT                 * INFO);
void dgetri_(
        const LAPACK_INT        * N,
        double              * A,
        const LAPACK_INT        * LDA,
        const LAPACK_INT           * IPIV,
        double              * WORK,
        const LAPACK_INT    * LWORK,
        LAPACK_INT                 * INFO);

//--------------------------------------------
//QRF routines:
void sgeqrf_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        float           * A,
        const LAPACK_INT    * LDA,
        float           * TAU,
        float           * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void dgeqrf_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        double          * A,
        const LAPACK_INT    * LDA,
        double          * TAU,
        double          * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void cgeqrf_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        std::complex<float>  * TAU,
        std::complex<float>  * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);
void zgeqrf_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        std::complex<double> * TAU,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//QP3 / RRQR routines:
void sgeqp3_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        float           * A,
        const LAPACK_INT    * LDA,
        LAPACK_INT             * JPVT,
        float           * TAU,
        float           * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void dgeqp3_(
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        double          * A,
        const LAPACK_INT    * LDA,
        LAPACK_INT             * JPVT,
        double          * TAU,
        double          * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void cgeqp3_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        LAPACK_INT                         * JPVT,
        std::complex<float>  * TAU,
        std::complex<float>  * WORK,
        const LAPACK_INT            * LWORK,
        float                       * RWORK,
        LAPACK_INT                         * INFO);
void zgeqp3_(
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        LAPACK_INT                         * JPVT,
        std::complex<double> * TAU,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        double                      * RWORK,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//MQR routines:
void sormqr_(
        char            * SIDE,
        char            * TRANS,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        const LAPACK_INT    * K,
        float           * A,
        const LAPACK_INT    * LDA,
        float           * TAU,
        float           * C,
        const LAPACK_INT    * LDC,
        float           * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void dormqr_(
        char            * SIDE,
        char            * TRANS,
        const LAPACK_INT    * M,
        const LAPACK_INT    * N,
        const LAPACK_INT    * K,
        double          * A,
        const LAPACK_INT    * LDA,
        double          * TAU,
        double          * C,
        const LAPACK_INT    * LDC,
        double          * WORK,
        const LAPACK_INT* LWORK,
        LAPACK_INT             * INFO);
void cunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        const LAPACK_INT                * K,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        std::complex<float>  * TAU,
        std::complex<float>  * C,
        const LAPACK_INT                * LDC,
        std::complex<float>  * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);
void zunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        const LAPACK_INT                * K,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        std::complex<double> * TAU,
        std::complex<double> * C,
        const LAPACK_INT                * LDC,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//SYEV/HEEV routines:
void ssyev_(
        char                        * JOBZ,
        char                        * UPLO,
        const LAPACK_INT                * N,
        float                       * A,
        const LAPACK_INT                * LDA,
        float                       * W,
        float                       * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);
void dsyev_(
        char                        * JOBZ,
        char                        * UPLO,
        const LAPACK_INT                * N,
        double                      * A,
        const LAPACK_INT                * LDA,
        double                      * W,
        double                      * WORK,
        const LAPACK_INT            * LWORK,
        LAPACK_INT                         * INFO);
void cheev_(
        char                        * JOBZ,
        char                        * UPLO,
        const LAPACK_INT                * N,
        std::complex<float>  * A,
        const LAPACK_INT                * LDA,
        float                       * W,
        std::complex<float>  * WORK,
        const LAPACK_INT            * LWORK,
        float                       * RWORK,
        LAPACK_INT                         * INFO);
void zheev_(
        char                        * JOBZ,
        char                        * UPLO,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        double                      * W,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        double                      * RWORK,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//STEV routines:
void sstev_(
        char                        * JOBZ,
        const LAPACK_INT                * N,
        float                       * D,
        float                       * E,
        float                       * Z,
        const LAPACK_INT                * LDZ,
        float                       * WORK,
        LAPACK_INT                         * INFO);

void dstev_(
        char                        * JOBZ,
        const LAPACK_INT                * N,
        double                      * D,
        double                      * E,
        double                      * Z,
        const LAPACK_INT                * LDZ,
        double                      * WORK,
        LAPACK_INT                         * INFO);

//--------------------------------------------
//LANGE norm routines:

float slange_(
        char                        * NORM,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        const float                 * A,
        const LAPACK_INT                * LDA,
        float                       * WORK );


double dlange_(
        char                        * NORM,
        const LAPACK_INT                * M,
        const LAPACK_INT                * N,
        const double                * A,
        const LAPACK_INT                * LDA,
        double                      * WORK );


float clange_(
        char                                * NORM,
        const LAPACK_INT                        * M,
        const LAPACK_INT                        * N,
        const std::complex<float>    * A,
        const LAPACK_INT                        * LDA,
        float                               * WORK );

double zlange_(
        char                                * NORM,
        const LAPACK_INT                        * M,
        const LAPACK_INT                        * N,
        const std::complex<double>   * A,
        const LAPACK_INT                        * LDA,
        double                              * WORK );

void zgeev_(
        char                        * JOBVL,
        char                        * JOBVR,
        const LAPACK_INT                * N,
        std::complex<double> * A,
        const LAPACK_INT                * LDA,
        std::complex<double> * W,
        std::complex<double> * VL,
        const LAPACK_INT                * LDVL,
        std::complex<double> * VR,
        const LAPACK_INT                * LDVR,
        std::complex<double> * WORK,
        const LAPACK_INT            * LWORK,
        double                      * RWORK,
        LAPACK_INT                         * INFO);

}


#endif // LAPACK_WRAPPER_H

