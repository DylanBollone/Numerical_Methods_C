/***********************************************************************
*
* This library provides stripped down, basic code supporting the kind
* of hands-on analysis of algorithms we're carrying out in this class.
* 
***********************************************************************/

#ifndef linalgh
#define linalgh

#include<time.h>

/***********************************************************************
* Calculate Ax and store the result in b
*
* Assumes A is mxn, x is dimension n, and all memory is properly
* allocated
*
* Note that memory management (eg freeing b) must be handled externally
***********************************************************************/
double* matrix_vector_multiply(double* A, double* x, int m, int n)
{
    int i, j;

    if(m <= 0 || n <= 0)
    {
        printf("Error: Invalid size passed to function matrix_vector_multiply. NULL returned.\n");
        return NULL;
    }

    double* b;
    b = malloc(m*sizeof(double));

    for(i = 0; i < m; i++)
    {
        b[i] = 0;
        for(j = 0; j < n; j++)
        {
            b[i] += A[i*n + j]*x[j];
        }
    }

    return b;
}

/***********************************************************************
* Create a matrix ones_u
*
* ones_u will be nxn with 1 in the upper triangle and zero in the lower
*
* Note that memory management (free, etc) must be handled outside of
* of the function
***********************************************************************/
double* ones_u(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function ones_u. NULL returned.\n");
        return NULL;
    }

    double* u;
    u = (double*) malloc(n*n*sizeof(double));

    int i, j;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(j>=i)
            {
                u[i*n + j] = 1.0;
            }
            else
            {
                u[i*n + j] = 0.0;
            }
        }
    }

    return u;
}

/***********************************************************************
* Create a matrix ones_l
*
* ones_l will be nxn with 1 in the lower triangle and on the diagonal
* ones_l will have zero in the upper triangle
*
* Note that memory management (free, etc) must be handled outside of
* of the function
***********************************************************************/
double* ones_l(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function ones_l. NULL returned.\n");
        return NULL;
    }

    double* l;
    l = (double*) malloc(n*n*sizeof(double));

    int i, j;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(j<=i)
            {
                l[i*n + j] = 1.0;
            }
            else
            {
                l[i*n + j] = 0.0;
            }
        }
    }

    return l;
}

/***********************************************************************
* Create a vector of ones 
*
* Note that memory management (free, etc) must be handled outside of
* of the function
***********************************************************************/
double* ones_vector(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function ones_vector. NULL returned.\n");
        return NULL;
    }

    double* x;
    x = (double*) malloc(n*sizeof(double));

    int i;

    for(i = 0; i < n; i++)
    {
        x[i] = 1.0;
    }

    return x;
}

/***********************************************************************
* Print an mxn matrix
***********************************************************************/
void print_matrix(double* A, int m, int n)
{
    int i, j;

    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            printf("%8.2f",A[i*n + j]);
        }
        printf("\n");
    }
}

/***********************************************************************
* Backwards substitution
*
* Starts from bottom right of matrix and iterates up the rows
*     no idea how to explain b_minus but it works
***********************************************************************/
double* backsub1(double* U, double* b, int n)
{
    if(n<=0)
    {
        printf("Error: Invalid size passed to function backsub. NULL returned.\n");
        return NULL;
    }
    
    int i,j;
    double* x;
    double b_minus;
    x = (double*) malloc(n*sizeof(double));

    for(i=n-1;i>-1;i--)
    {
        b_minus = 0.0;
        for(j=i;j<n-1;j++)
        {
            b_minus += U[i*n+j+1]*x[j+1];
        }
        x[i] = (b[i] - b_minus) / U[i*n+i];
    }

    return x;
}

/*********************************************************************
* Forward substitution
*
* Should have done this one first
*********************************************************************/
double* forsub2(double* L, double* b, int n)
{
    if(n<=0)
    {
        printf("ERROR: Invalid size passed to function forsub. NULL returned\n");
        return NULL;
    }

    int i,j;
    double* x;
    double b_minus;
    x = (double*) malloc(n*sizeof(double));

    for(i=0;i<n;i++)
    {
        b_minus = 0.0;
        for(j=0;j<i;j++ )
        {
            b_minus += L[i*n+j]*x[j];
        }
        
        x[i] = (b[i] - b_minus) / L[i*n+i];
    }
    return x;
}

/***********************************************************************
* Create a random nxn matrix with entries between -1.0 and 1.0 off
* the diagonal and n on the diagonal
***********************************************************************/
double* random_matrix(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function random_matrix. NULL returned.");
        return NULL;
    }
    
    int i, j;
    double* A;
    double r;
    A = (double*) malloc(n*n*sizeof(double));

    srand(time(NULL));
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            if(i == j)
            {
                A[i*n + j] = n;
            }
            else
            {
                r = rand(); //r will be an integer between 0 and RAND_MAX
                r = 2*r;    //r will be an even integer between 0 and 2*RAND_MAX
                r = r / (double) RAND_MAX; //r will be a double between 0 and 2.0
                r = (r - 1.0) * 10.0;    //r will be a double between -10 and 10
                A[i*n + j] = r;
            }
    }

    return A;
}

double* random_vector(int n)
{
    int i;
    double* v;
    double r;
    v = (double*) malloc(n*sizeof(double));

    srand(time(NULL));
    for(i=0;i<n;i++)
    {
        r = rand(); //r will be an integer between 0 and RAND_MAX
        r = 2*r;    //r will be an even integer between 0 and 2*RAND_MA  
        r = r / (double) RAND_MAX; //r will be a double between 0 and 2.0
        r = r - 1.0;    //r will be a double between -1 and 1
        v[i] = r;
    }
    return v;
}

/***********************************************************************
* Create a random nxn matrix with entries equal to zero or one with
* equal probability
*
* For small values of n, these matrices are sometimes singular
* For large values of n, they may still be singular but the probability
* of this declines exponentially
***********************************************************************/
double* random_binary_matrix(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function random_matrix. NULL returned.");
        return NULL;
    }

    int i, j;
    double* A;
    double r;
    A = (double*) malloc(n*n*sizeof(double));

    srand(time(NULL));
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            r = (double) rand() / (double) RAND_MAX;
            if(r < 0.5)
            {
                A[i*n + j] = 0.0;
            }
            else
            {
                A[i*n + j] = 1.0;
            }
        }
    }

    return A;
}

/***********************************************************************
* Create a random nxn upper triangular matrix with entries between 
* -1.0 and 1.0 off the diagonal and n on the diagonal
***********************************************************************/
double* random_u_matrix(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function random_u_matrix. NULL returned.");
        return NULL;
    }
    
    int i, j;
    double* A;
    double r;
    A = (double*) malloc(n*n*sizeof(double));

    srand(time(NULL));
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(i == j)
            {
                A[i*n + j] = n;
            }
            else if(i < j)
            {
                r = rand(); //r will be an integer between 0 and RAND_MAX
                r = 2*r;    //r will be an even integer between 0 and 2*RAND_MAX
                r = r / (double) RAND_MAX; //r will be a double between 0 and 2.0
                r = r - 1.0;    //r will be a double between -1 and 1
                A[i*n + j] = r;
            }
            else
            {
                A[i*n + j] = 0.0;
            }
        }
    }

    return A;
}

/***********************************************************************
* Create a random nxn upper triangular matrix with entries between 
* -1.0 and 1.0 off the diagonal and 1 on the diagonal
***********************************************************************/
double* random_l_matrix(int n)
{
    if(n <= 0)
    {
        printf("Error: Invalid size passed to function random_l_matrix. NULL returned.");
        return NULL;
    }
    
    int i, j;
    double* A;
    double r;
    A = (double*) malloc(n*n*sizeof(double));

    srand(time(NULL));
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(i == j)
            {
                A[i*n + j] = 1.0;
            }
            else if(i > j)
            {
                r = rand(); //r will be an integer between 0 and RAND_MAX
                r = 2*r;    //r will be an even integer between 0 and 2*RAND_MAX
                r = r / (double) RAND_MAX; //r will be a double between 0 and 2.0
                r = r - 1.0;    //r will be a double between -1 and 1
                A[i*n + j] = r;
            }
            else
            {
                A[i*n + j] = 0.0;
            }
        }
    }

    return A;
}
/****************************************************************************************
*
* Print vector
*     prints a vector of size n
*
****************************************************************************************/
void print_vector(double* v, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
	printf("%8.2f\n",v[i]);
    }
    printf("\n");
}


#endif
