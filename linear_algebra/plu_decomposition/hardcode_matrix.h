/*******************************************************************
*
* This is a library soley for the purpose of hardcoding matrices
*
* Will include 4 functions for A,L,U,x,b of size n=3
*
*******************************************************************/


#ifndef hardcode_matrixh
#define hardcode_matrixh

// n will be used many times so for ease of use changing later
// n will be global variable (I know ew right)
int n = 3;

/*******************************************************************
*
* U matrix
*
*******************************************************************/
double* U_matrix()
{
    double* U;
    U = (double*) malloc(n*n*sizeof(double));
    
    // diagonal terms
    // U[0,0]
    U[0] = 2.0;
    // U[1,1]
    U[4] = 2.0;
    // U[2,2]
    U[8] = 3.0;

    // terms below diagonal = 0
    // U[1-2,0]
    U[3] = 0.0;
    U[6] = 0.0;
    // U[2,1]
    U[7] = 0.0;

    // back sub constants
    U[5] = 1.0;
    U[2] = 2.0;
    U[1] = 4.0;

    return U;

}

/*******************************************************************
*
* b matrix
*
*******************************************************************/
double* b_matrix()
{
    double* b;
    b = (double*) malloc(n*sizeof(double));

    b[0] = 2.0;
    b[1] = 4.0;
    b[2] = 8.0;

    return b;
}


/*******************************************************************
*
* L matrix
*
*******************************************************************/
double* L_matrix()
{
    double* L;
    L = (double*) malloc(n*n*sizeof(double));

    // diagonal terms = 1.0
    // L[0,0]
    L[0] = 1.0;
    // U[1,1]
    L[4] = 1.0;
    // U[2,2]
    L[8] = 1.0;

    // terms above diagonal = 0
    L[1] = 0.0;
    L[2] = 0.0;
    L[5] = 0.0;

    // forsub constants
    L[3] = 3.0;
    L[6] = 4.0;
    L[7] = 2.0;

    return L;

}
/*******************************************************************
*
* A matrix
*
*******************************************************************/
double* A_matrix()
{
    double* A;
    A = (double*) malloc(4*4*sizeof(double));
    A[0] = 2.0;
    A[1] = 3.0;
    A[2] = 1.0;
    A[3] = 5.0;
    A[4] = 5.0;
    A[5] = -8.0;
    A[6] = -15.0;
    A[7] = 7.0;
    A[8] = -6.0;
    A[9] = 12.0;
    A[10] = 18.0;
    A[11] = 6.0;
    A[12] = 1.0;
    A[13] = 1.0;
    A[14] = 1.0;
    A[15] = 1.0;

    return A;
}

#endif
