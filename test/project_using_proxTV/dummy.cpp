#include "TVopt.h"


int main()
{
// int DR2_TV(size_t M, size_t N, double*unary, double W1, double W2, double norm1, double norm2, double*s, int nThreads, int maxit, double* info);

    size_t M = 1;
    size_t N = 1;
    double* unary = new double;
    double W1 = 1.0;
    double W2 = 1.0;
    double norm1 = 1.0;
    double norm2 = 1.0;
    double* s = new double;
    int nThreads = 1;
    int maxit = 1;
    double* info = new double;
    int r = DR2_TV(M, N, unary, W1, W2, norm1, norm2, s, nThreads, maxit, info);

    delete unary;
    delete s;
    delete info;

    return r;
}
