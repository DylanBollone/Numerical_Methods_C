/******************************************************************************
* Exam One: Take Home Portion
*
* To complete the take home portion of exam one, modify this file as follows:
*   (1) Place a working back substitution routine where I've indicated
*   (2) Likewise with forward substitution
*   (3) Write and test the function lu_expansive(...) outlined below.
*       This function must compute p, L, and U from A and n;
*   (4) Write a main.c file that demos the use of LU decomposition to solve
*       random systems of equations.
*   (5) Stop in at my office and explain/demo your code.
******************************************************************************/

#include<math.h> //we need this for the fabs function

//Change the name of the header and the two lines below to match your name
#ifndef dylanlibh
#define dylanlibh

//Put your backsub routine here
double* backsub(double* U, double* y, int n)
{
    if(n<=0)
    {
        printf("Error: Invalid size passed to function backsub. NULL returned.\n");
        return NULL;
    }

    int i,j;
    double* x;
    double y_minus;
    x = (double*) malloc(n*sizeof(double));

    for(i=n-1;i>-1;i--)
    {
        y_minus = 0.0;
        for(j=i;j<n-1;j++)
        {
            y_minus += U[i*n+j+1]*x[j+1];
        }
        x[i] = (y[i] - y_minus) / U[i*n+i];
    }

    return x;
}

//Put your forsub routine here
double* forsub(double* L, double* b, int n)
{
    if(n<=0)
    {
        printf("ERROR: Invalid size passed to function forsub. NULL returned\n");
        return NULL;
    }

    int i,j;
    double* y;
    double b_minus;
    y = (double*) malloc(n*sizeof(double));

    for(i=0;i<n;i++)
    {
        b_minus = 0.0;
        for(j=0;j<i;j++ )
        {
            b_minus += L[i*n+j]*y[j];
        }

        y[i] = (b[i] - b_minus) / L[i*n+i];
    }
    return y;
}

/******************************************************************************
* Expansive LU factorization with partial pivoting
*
* "Expansive" means we will explicitly store all entries of L and U
*
* To use this:
*   Define pointers double* A; int* p; double* L; double* U;
*   Make A into the matrix you want to factor;
*   Use the function call lu_expansive(A,n,&p,&L,&U);
*
* This passes pointers to the relevant pointers, hence the double stars
* in the definition.
******************************************************************************/
void lu_expansive(double* A, int n, int** rp, double** rL, double** rU)
{
    if(n <= 0)
    {
        printf("Invalid dimension sent to lu_expansive.\n");
        return;
    }

    //Inside this function, I'm going to do stuff to p, L, and U
    int i,j,k;
    int* p;
    double* L;
    double* U;

    p = (int*) malloc(n*sizeof(int));
    L = (double*) malloc(n*n*sizeof(double));
    U = (double*) malloc(n*n*sizeof(double));

    /**************************************************************************
    * Do stuff to p, L, U
    **************************************************************************/
    // p will start as a normal counting array
    //     [0 1 2 3 ... n-1]
    for(i=0;i<n;i++)
    {
        p[i] = i;
    }

    // U will start as A
    //     not as simple as U = A
    //     will cause "aborted" error after whole program run
    //     took me two days to find out why my program was aborting    
    for(i=0;i<n;i++)
    {
	for(j=0;j<n;j++)
	{
	    U[i*n+j] = A[i*n+j];
	}
    }

    // L will start as empty identity matrix
    //     Ones on the diagonal and rest is zeros
    for(i=0;i<n;i++)
    {
	for(j=0;j<n;j++)
	{
	    if(j==i)
	    {
		L[i*n+j] = 1.0;
	    } else if(j!=i) {
		L[i*n+j] = 0.0;
	    }
	}
    }

    // these vars are used for matrix operations and recording values
    // associated with those matrix operations
    double high_abs,temp,m;
    int mrow;

    for(i=0;i<n;i++)
    {
	// in each column find the highest absolute valule
	//     start with absolute value of pivot element
	high_abs = fabs(U[i*n+i]);
	mrow = i;
	for(j=i+1;j<n;j++)
	{
	    if(fabs(U[j*n+i]) > high_abs)
	    {
		high_abs = fabs(U[j*n+i]);
		mrow = j;
	    }
	}
	// now that I which row has highest absolute value
	// swap rows, using a temp var?
	// only if a row swap needs to occur, otherwise it will just swap a row with itself
	if( mrow != i)
	{
	    // swaps values in U
	    for(j=0;j<n;j++)
	    {
	        temp = U[i*n+j];
	        U[i*n+j] = U[mrow*n+j];
	        U[mrow*n+j] = temp;
	    }
	    // in L: 
	    //     swap rows i and mrow only columns 0 through i-1
	    for(j=0;j<i;j++)
	    {
		temp = L[i*n+j];
		L[i*n+j] = L[mrow*n+j];
		L[mrow*n+j] = temp;
	    }
	    // record the swap in p
	    p[i] = p[mrow];
	    p[mrow] = i;
	}

	// elminate terms in row i below U[i,i]
	//     record the multipliers in L matrix
	for(j=i+1;j<n;j++)
	{
	    m = U[j*n+i] / U[i*n+i];
	    for(k=i;k<n;k++)
	    {
		U[j*n+k] = U[j*n+k] - (U[i*n+k] * m);
	    }
	    L[j*n+i] = m;
	}
    }


    // now that p, L, and U are primed and ready
    // start my solving process

    // initialize all neccessary vars
    double* pb;
    double* b;
    double* y;
    double* x;
    b = f_homework7(n);
    pb = (double*) malloc(n*sizeof(double));   

    // print out b
//    printf("b vector:\n");
//    print_vector(b,n);
//    printf("\n");

    // make p^-1b equal to b put with the rows swapped as acccording to p
    //     accounts for row swapping in L and U by swapping rows in b
    for(i=0;i<n;i++)
    {
        pb[i] = b[p[i]];
    }


    // forward and backward substitution functions allocate memory for me
    // all I have to do is free memory for recent vars at end of this whole function

    // solve Ly = pb using forsub
    y = forsub(L,pb,n);
//    printf("y vector:\n");
//    print_vector(y,n);
//    printf("\n");
//
    // solve Ux = y using backsub
    x = backsub(U,y,n);
    // prints final result (x)
//    printf("x vector:\n");
//    print_vector(x,n);
//    printf("\n");


    // now check if it works
    //     by multiplying A (which has not changed) by x (final calculated vector)
    //     This vector should be exactly equal to b
    //         barring any inacuacies due to doubles not technically being infinetly precise
    double* Ax;
    Ax = matrix_vector_multiply(A,x,n,n);
//    printf("Ax:\n");
//    print_vector(Ax,n);
//    printf("\n");

    /**************************************************************************
    * Ok, now p, L, U are ready. Store their pointer values in the memory
    * pointed to by rp, rL, rU (dereference the pointers to pointers and store
    * pointers to our prepared objects)
    **************************************************************************/
    *rp = p;
    *rL = L;
    *rU = U;

    // in order to solbe Ax = b I created arrays of doubles to store values
    //     must free memory in function because I am done with the information
    //     as it is already printed out for user
    free(b);
    free(pb);
    free(y);
    free(x);
    free(Ax);

}


#endif
