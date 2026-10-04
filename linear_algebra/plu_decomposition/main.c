#include <stdlib.h>
#include <stdio.h>
#include "linalg.h" // matrix printing and random matrices
#include "dylanlib.h" // my version of exam1 library
#include "hardcode_matrix.h" // testing purposes

// quick note:
//     forsub and backsub are both identical procedures in linalg.h and dylanlib.h
//     changed a few letters and used the routines in dylanlib.h

int main()
{
    // these variables are needed for function calls
    //     functions also have variables named inside and are printed out in function
    double* A;
    int* p;
    double* L;
    double* U;
    // note: I could prompt user or just set n before program runs
    int n = 5,i;

    // A is only matrix called for by lu_expansive that is not going to be assigned values in func
    A = random_matrix(n);

    // does stuff to p,L,U (factorization)
    // changes p, L, and U
    // calculates and prints b, x, and Ax (which should be equal to b)
    lu_expansive(A,n,&p,&L,&U);

    // print A matrix for demo
    printf("A matrix:\n");
    print_matrix(A,n,n);
    printf("\n");

    // print p for demo
    printf("p vector:\n");
    // quick side note:
    //     I probably should make a function that prints vectors
    for(i=0;i<n;i++)
    {
        printf("%8d\n",p[i]);
    }
    printf("\n");

    // print L matrix for demo
    printf("L matrix:\n");
    print_matrix(L,n,n);
    printf("\n");
    
    // print U matrix for demo
    printf("U matrix:\n");
    print_matrix(U,n,n);   

    free(A);
    free(p);
    free(L);
    free(U);

    return 0;
}
