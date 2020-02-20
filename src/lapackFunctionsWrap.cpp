/**
    Definitions of functions wrapping lapack.
    Declarations depends on compilations options PROXTV_USE_LAPACK

    @author Pablo Hernandez-Cerdan
*/

#include "lapackFunctionsWrap.h"

#ifndef PROXTV_USE_LAPACK // USE_EIGEN
#include <Eigen/Dense>
void dpttrf_plus_dpttrs_eigen( lapack_int* n, double* d, double* e, double *b)
{
    using EigenMatrix = Eigen::MatrixXd;
    using EigenVector = Eigen::VectorXd;
    using EigenVectorMap = Eigen::Map<EigenVector>;
    // Eigen has to create the full matrix from the diagonal d and subdiagonal e
    int mSize = *n;
    EigenMatrix eigenM(mSize,mSize);
    EigenVectorMap diag(d, mSize);
    EigenVectorMap subAndUpperDiag(e, mSize - 1);
    EigenVectorMap inputB_outputX(b, mSize);
    // Populate matrix
    eigenM.diagonal() = diag;
    eigenM.diagonal( 1) = subAndUpperDiag;
    eigenM.diagonal(-1) = subAndUpperDiag;

    // Factorize using ldlt (ldl is also possible, faster but less accurate)
    // A = LDL'
    Eigen::LDLT<EigenMatrix> factorization(eigenM);
    // A*X = b
    // This modifies the input/output pointer: b
    inputB_outputX = factorization.solve(inputB_outputX);
    // This modifies the input/output pointers: d and e
    EigenMatrix factorized = factorization.matrixLDLT();
    diag = factorized.diagonal();
    subAndUpperDiag = factorized.diagonal(-1);
}
#endif
