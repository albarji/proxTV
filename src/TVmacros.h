/**
    General macros for Total Variation algorithms
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/

// Compute dual vector x from primal vector w
// A counter i must also be provided
#define DUAL2PRIMAL(w,x,i) \
    x[0] = y[0]+w[0]; \
    for(i=1;i<nn;i++) \
        x[i] = y[i]-w[i-1]+w[i]; \
    x[nn] = y[nn]-w[nn-1];

// Compute gradient vector g from dual vector x
// A counter i must also be provided    
#define PRIMAL2GRAD(x,g,i) \
    for(i=0;i<nn;i++) \
        g[i] = x[i] - x[i+1];

// Compute objective value (val) of dual problem from primal vector w and input signal y
// A counter i must also be provided            
#define DUALVAL(w,y,val,i) \
    val = -w[0]*(0.5*(-w[0])-y[0]); \
    for(i=1;i<nn;i++) \
        val += (w[i-1]-w[i])*(0.5*(w[i-1]-w[i])-y[i]); \
    val += w[nn-1]*(0.5*w[nn-1]-y[nn]);
