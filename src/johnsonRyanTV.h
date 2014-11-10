#include <stdlib.h> 
using namespace std;

// Dynamic programming algorithm for the 1d fused lasso problem
// (Ryan's implementation of Nick Johnson's algorithm)
void dp(int n, double *y, double lam, double *beta);
