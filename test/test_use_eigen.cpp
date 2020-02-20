#include <Eigen/Dense>
#include <iostream>
#include "general.h"

using EigenMatrix = Eigen::MatrixXd;
using EigenVector = Eigen::VectorXd;

EigenMatrix oneD()
{
    EigenMatrix m(1,4);
    m(0,0) = 3;
    m(0,1) = 2.5;
    m(0,2) = -1;
    m(0,3) = 1.5;
    return m;
}

EigenMatrix lapack_oneD()
{
    auto eigenM = oneD();
    int n = eigenM.cols();
    std::cout << "Cols: " << n << std::endl;
    int nn = n ;
    int rc = nn;
    int nnp = nn;
    int one = 1;
    double *w = nullptr;
    double *aux = nullptr;
    double *aux2 = nullptr;
    w = (double*)malloc(sizeof(double)*nn);
    aux = (double*)malloc(sizeof(double)*nn);
    aux2 = (double*)malloc(sizeof(double)*nn);
    for(int i=0;i<nn;i++)
      w[i] = eigenM(0, i+1) - eigenM(0, i); /* Dy */

    for(int i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;

    dpttrf_(&nnp, aux, aux2, &rc);
    std::cout << "LDL factorization:" << std::endl;
    std::cout << "aux (Diagonal):" << std::endl;
    for(int i=0;i<nn;i++) {
      std::cout << aux[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "aux2 (SubDiagonal):" << std::endl;
    for(int i=0;i<nn;i++) {
      std::cout << aux2[i] << ", ";
    }
    std::cout << std::endl;
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);

    std::cout << "Solution w:" << std::endl;
    for(int i=0;i<nn;i++)
      std::cout << w[i] << ", ";
    std::cout << std::endl;

    EigenVector m(4);
    for(int i=0;i<nn;i++)
      m(i) = w[i];
    std::cout << m << std::endl;

    if(w) free(w);
    if(aux) free(aux);
    if(aux2) free(aux2);

    return m;
}

/* https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html */
EigenMatrix eigen_lapack_equivalent_oneD()
{
    // Solve X in A*X = B, where A is triangular definite
    EigenMatrix eigenM(4,4);
    {
        EigenVector diag(eigenM.cols());
        diag.setConstant(2);
        EigenVector subUpperDiag(eigenM.cols());
        subUpperDiag.setConstant(-1);
        eigenM.diagonal() = diag;
        eigenM.diagonal(1) = subUpperDiag;
        eigenM.diagonal(-1) = subUpperDiag;
    }
    std::cout << "Matrix eigenM:\n" << eigenM << std::endl;
    EigenVector b(eigenM.cols());
    EigenMatrix oneM = oneD();
    for (int i = 0; i < oneM.cols(); ++i) {
        b[i] = oneM(0, i+1) - oneM(0, i); /* Dy */
    }

    // EigenVector v = eigenM.llt().solve(b);
    // EigenVector v = eigenM.ldlt().solve(b);
    Eigen::LDLT<EigenMatrix> ldltOfM(eigenM);
    EigenVector v = ldltOfM.solve(b);
    std::cout << v << std::endl;
    // EigenVector v = eigenM.colPivHouseholderQr().solve(b);
    EigenMatrix ldltMatrix = ldltOfM.matrixLDLT();
    std::cout << ldltMatrix << std::endl;
    EigenMatrix lMatrix = ldltOfM.matrixL();
    std::cout << lMatrix << std::endl;
    std::cout << "Diagonal of LDLT:" << std::endl;
    std::cout << ldltMatrix.diagonal() << std::endl;
    std::cout << "SubDiagonal of LDLT:" << std::endl;
    std::cout << ldltMatrix.diagonal(-1) << std::endl;

    return v;
}

/// The idea is that input and output are raw pointers to couple with existing code
EigenMatrix eigen_interface_raw_pointers_oneD()
{
    // Create the typemaps
    using EigenMatrixMap = Eigen::Map<EigenMatrix>;
    using EigenMatrixReadOnlyMap = Eigen::Map<const EigenMatrix>;
    using EigenVectorMap = Eigen::Map<EigenVector>;
    using EigenVectorReadOnlyMap = Eigen::Map<const EigenVector>;
    // Populate values using Eigen (just because easier)
    EigenMatrix eigenM(4,4);
    {
        EigenVector diag(eigenM.cols());
        diag.setConstant(2);
        EigenVector subUpperDiag(eigenM.cols());
        subUpperDiag.setConstant(-1);
        eigenM.diagonal() = diag;
        eigenM.diagonal(1) = subUpperDiag;
        eigenM.diagonal(-1) = subUpperDiag;
    }
    std::cout << "Matrix eigenM:\n" << eigenM << std::endl;
    EigenVector b(eigenM.cols());
    EigenMatrix oneM = oneD();
    for (int i = 0; i < oneM.cols(); ++i) {
        b[i] = oneM(0, i+1) - oneM(0, i); /* Dy */
    }
    // Set input raw pointers
    int nn = eigenM.cols();
    double *w = (double*)malloc(sizeof(double)*nn);
    for(int i=0;i<nn;i++) w[i] = b[i];
    // w is b as input (Ax=b), and solution x as output in lapack code.
    // http://www.netlib.org/lapack/lapack-3.1.1/html/dpttrs.f.html
    // Create a VectorMap of w to interact with Eigen
    EigenVectorMap wMap(w, b.size());
    EigenVector v = eigenM.ldlt().solve(wMap);
    std::cout << "v:\n" << v << std::endl;
    std::cout << "b:\n" << b << std::endl;
    std::cout << "wMap:\n" << wMap << std::endl;
    wMap = eigenM.ldlt().solve(wMap);
    std::cout << "InPlace wMap:\n" << wMap << std::endl;
    v = wMap;
    if(w) free(w);
    return v;
}

int main()
{
    bool test_passed = true;
    std::cout << "-- oneD Eigen --" << std::endl;
    auto mOneD = oneD();
    std::cout << mOneD << std::endl;
    std::cout << "-- lapack (from PN_TV1_Weighted) --" << std::endl;
    auto mLapack = lapack_oneD();
    std::cout << "-- eigen equivalent cholesky --" << std::endl;
    auto mEigen = eigen_lapack_equivalent_oneD();
    if(mLapack != mEigen){
        std::cerr << "Result vectors between mLapack and mEigen are not equal"  << std::endl;
        test_passed = false;
    }
    std::cout << "-- eigen interface with raw pointers --" << std::endl;
    auto mEigenRaw = eigen_interface_raw_pointers_oneD();
    if(mLapack != mEigenRaw){
        std::cerr << "Result vectors between mLapack and mEigenRaw are not equal"  << std::endl;
        test_passed = false;
    }

    if(!test_passed) return EXIT_FAILURE;
}
