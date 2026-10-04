#ifndef heath
#define heath

//Returns the full diffusion matrix with rows ... 0 1 -2 1 0 ...
double* diffusion_matrix(int n)
{
    if(n <= 0)
    {
        printf("Invalid dimension sent to diffusion_matrix. NULL returned.\n");
        return NULL;
    }

    int i, j;
    double dhs = (n+1)*(n+1); //a multiply will divide by h^2
    double* A;
    A = (double*) malloc(n*n*sizeof(double));

    for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
    {
        if(i == j - 1 || i == j+1)
            A[i*n + j] = -dhs;
        else if(i == j)
            A[i*n + j] = 2 * dhs;
        else
            A[i*n + j] = 0.0;
    }

    return A;
}

//multiply x by the diffusion matrix (without the matrix)
double* sparse_multiply(double* x, int n)
{
    if(n <= 2)
    {
        printf("Invalid dimension sent to sparse_multiply. NULL returned.\n");
        return NULL;
    }

    int i;
    double dhs = (n+1)*(n+1); //a multiply will divide by h^2
    double* b;
    b = (double*) malloc(n*sizeof(double));

    b[0]   = dhs*(2*x[0]   - x[1]);

    b[n-1] = dhs*(2*x[n-1] - x[n-2]);

    for(i = 1; i < n-1; i++)
        b[i] = dhs*(2*x[i] - x[i-1] - x[i+1]); 

    return b;

}


#endif
