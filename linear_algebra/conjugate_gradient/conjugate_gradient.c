#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"
#include "heat.h"
#include "hires_timer.c"

double vector_length(double* v, int n)
{
    double len=0.0;
    int i;

    for(i=0;i<n;i++)
    {
	len = len + (v[i]*v[i]);
    }
    len = sqrt(len);
    return len;
}

double dot(double* a, double* b, int n)
{
    double sum = 0.0;
    int i;
    for(i=0;i<n;i++)
    {
	sum = sum + (a[i] * b[i]);
    }
    return sum;
}

double* heat_source(int n)
{
    double* f;
    int i;
    f = (double*) malloc(n*sizeof(double));
    for(i=0;i<n;i++)
    {
	if(i < (n/3) || i >= (2*n/3))
	{
	    f[i] = 0.0;
	} else {
	    f[i] = 1.0;
	}
    }
    return f;
}

int main()
{
    int n=100000,i;
//    double* A;
    double* x;
    double* b;
    double* r;
    double* p;
    double* A_times_p;
    double alpha;
    double beta;
    double epsilon = pow(10,-3);

    double ts, tf, t;
    init_hires_timer();
    hires_timer(&ts);
    
    // initialize A, x, and b
    //     for conjugate gradient method x is a guess, my guess is all ones
//    A = diffusion_matrix(n);
    x = ones_vector(n);
    b = heat_source(n);
#if 0
    //print matrices
    printf("Matrices for size n=%d\n",n);
    printf("\nA:\n");
    print_matrix(A,n,n);
    printf("\nb:\n");
    print_vector(b,n);
#endif

    // memory allocation for vars that do not rely of funcs to do so
    p = (double*) malloc(n*sizeof(double));

    // initialize r to A*x, then subtract b as shown in packet
    //     also sets p to -r
    r = sparse_multiply(x,n);
    for(i=0;i<n;i++)
    {
	r[i] = r[i] - b[i];
	p[i] = -r[i];
    }

    while(vector_length(r,n) > epsilon)
    {
	A_times_p = sparse_multiply(p,n);
	// calculate alpha
	//     didn't remember how r(k)t * p(k) was supposed to be calculated
	//     made a dot product function
	alpha = - (dot(r,p,n) / dot(p,A_times_p,n));
	
	// calculate new x
	for(i=0;i<n;i++)
	{
	    x[i] = x[i] + alpha*p[i];
	}
	
	// calculate new r
	free(r);
	r = sparse_multiply(x,n);
	for(i=0;i<n;i++)
	{
	    r[i] = r[i] - b[i];
	}

	// calculate beta
	beta = dot(r,A_times_p,n) / dot(p,A_times_p,n);

	// calculate new p
	for(i=0;i<n;i++)
	{
	    p[i] = beta*p[i] - r[i];
	}
	free(A_times_p);
    }

#if 0
    printf("Final x matrix:\n");
    print_vector(x,n);

    printf("Ax to prove program works:\n");
    print_vector(matrix_vector_multiply(A,x,n,n),n);
#endif

    hires_timer(&tf);
    t = tf - ts;
    printf("%lf",t);

//    free(A);
    free(x);
    free(b);
//    free(r);
    free(p);
//    free(A_times_p);

    return 0;
}
