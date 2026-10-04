#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x)
{
    return 1.0/x;
}

/******************************************************************
* Gauss-Kronrod Nodes and Weights
******************************************************************/
const double nodes[15] = {
-0.991455371120813,
-0.949107912342759,
-0.864864423359769,
-0.741531185599394,
-0.586087235467691,
-0.405845151377397,
-0.207784955007898,
 0.000000000000000,
 0.207784955007898,
 0.405845151377397,
 0.586087235467691,
 0.741531185599394,
 0.864864423359769,
 0.949107912342759,
 0.991455371120813
};

const double gw[7] = {
0.129484966168870,
0.279705391489277,
0.381830050505119,
0.417959183673469,
0.381830050505119,
0.279705391489277,
0.129484966168870
};

const double kw[15] = {
0.022935322010529,
0.063092092629979,
0.104790010322250,
0.140653259715525,
0.169004726639267,
0.190350578064785,
0.204432940075298,
0.209482141084728,
0.204432940075298,
0.190350578064785,
0.169004726639267,
0.140653259715525,
0.104790010322250,
0.063092092629979,
0.022935322010529
};

/**************************************************
 * Dylan Bollone                                  *
 * Gauss                                          *
 * Sums the gauss points * guass weights          *
 *************************************************/
double gauss(double *fvals)
{
    double sum = 0.0;
    int i;
    for(i=0;i<7;i++)
    {
	sum = sum + fvals[2*i+1]*gw[i];
    }
    return sum;
}

/*************************************************
 * Dylan Bollone                                 *
 * Kronrod                                       *
 * Sums all points * kronrod weights             *
 ************************************************/
double kronrod(double *fvals)
{
    double sum = 0.0;
    int i;
    for(i=0;i<15;i++)
    {
	sum = sum + fvals[i]*kw[i];
    }
    return sum;

}

/*************************************************
 * Dylan Bollone
 * gk - Gauss-Kronrod
 * Gauss-Kronrod 7-15 scheme
 ************************************************/
double gk(double a, double b, double tol, int DEPTH, int MAX)
{
    // calculate h and c
    double h = (b-a)/2.0;
    double c = (b+a)/2.0;
    double fx[15],G, K;
    int i;

    // evaluate the function on curr [a,b]
    for(i=0;i<15;i++)
    {
	fx[i] = f(c + nodes[i]*h);
    }

    // find gauss estimate
    G = h*gauss(fx);

    // find gauss-kronrod estimate
    K = h*kronrod(fx);

    // compare
    //     if meets tol then return kronrod estimate
    //     or if we are at the max depth
    if(fabs(K-G)<tol || DEPTH >= MAX)
    {
	return K;
    }

    // else do recusion calling this function again
    return gk(a,c,tol/2.0,DEPTH+1,MAX) + gk(c,b,tol/2.0,DEPTH+1,MAX);
}

int main()
{
    double estimate,x,tol;
    double a = 1.0;
    int MAX_DEPTH;

    // get x value
    printf("Enter a x value between 1 and 4: ");
    scanf("%lf",&x);
    while(x<1 || x>4)
    {
	printf("ERROR: Invalid input for x\n");
	printf("Enter a new value for x: ");
	scanf("%lf",&x);
    }


    // initialize tol and max depth
    tol = pow(10,-10);
    MAX_DEPTH = 15;

    // get estimate for ln(x)
    estimate = gk(a,x,tol,0,MAX_DEPTH);
    
    // output estimate
    printf("ln(%6.4lf) = %12.10lf\n",x,estimate);
    
    return 0;
}
