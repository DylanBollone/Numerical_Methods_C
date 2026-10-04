/****************************************
 * Dylan Bollone                        *
 * MATH 261                             *
 * Project 2                            *
 * 4/23/2025                            *
 ***************************************/

/* math.h is used for as few functions as possible (fabs and pow) */
#include <stdio.h>
#include <math.h>

/***************************************************
 * Dylan Bollone                                   *
 * Square root function from earlier this semester *
 * I believe this is from Homework 7               *
 **************************************************/
double sq_rt(double a)
{
    double x;
    double error = 1.0;
    if (a<0)
    {
        printf("ERROR: invalid input to sq_rt function.");
    } else if (a==0) {
        return a;
    }

    x = a / 2.0;

    while (fabs(error) > pow(10.0,-8.0)) {
        error = (a - x*x) / (2*x);
        x = x + error;
    }

    return x;
}

/*********************************************************
 * Dylan Bollone                                         *
 * Homework 13                                           *
 * My functions that get pi to same precision as math.h  *
 ********************************************************/
double fpi(double x)
{
    return 4.0/(1+x*x);
}
/* coarse simpsons rule */
double S(double a, double b)
{
    double h = (b-a)/2.0;
    double c = (a+b)/2.0;
    return (h/3.0) * (fpi(a) + 4*fpi(c) + fpi(b));
}

/* recursion of coarse simpsons rule */
double fine(double a, double b, double tol, int curr, int max)
{
    double Sf,c;
    double Sc = S(a,b);           /* calculate simpsons rule */

    /* if we've reached the max recursion depth return coarse simpson*/
    if (curr > max)
    {
        return Sc;
    } else {
        c = (a+b)/2.0;         /* find midpoint */
        Sf = S(a,c) + S(c,b);  /* calculate fine simpsons rule*/
        curr+=1;         /* update recursion depth */

        /* if we meet tol then return fine simpsons plus error */
        if ( fabs(Sc - Sf) < (15*tol) )
        {
            return Sf + (Sf-Sc)/15.0;
        }
        /* if not good then call function again */
        else
        {
            return fine(a,c,tol/2.0,curr,max) + fine(c,b,tol/2.0,curr,max);
        }
    }
}

/*************************************************
 * Dylan Bollone                                 *
 * e^x function                                  *
 * evaluates e^x for |x|<=700                    *
 * From Homework 12 or 13                        *
 ************************************************/
double my_exp(double x)
{
    double sum=1,t=1; // t is the term of the series that is being added to sum
    int squares=0,i; // number of times repeated squaring needs to occur
    // checks for my boundary conditions on input x
    //    absolute value of x is greater than 700 e^x is infinity or zero
    if (x>700)
    {
        return INFINITY;
    }
    else if (x<-700)
    {
        return 0;
    }

    // divides x by two and counts how many times until -1<x<1
    while (fabs(x)>1)
    {
        x = x / 2;
        squares = squares + 1;
    }
    // 20 ish terms of the taylor series approximation for e^x
    for (i=1;i<20;i++)
    {
        t = t * x / i;
        sum = sum + t;
    }
    // repeated squaring if applicable
    for (i=0;i<squares;i++)
    {
        sum = sum * sum;
    }
    return sum;
}


/*************************************************
 * Dylan Bollone                                 *
 * Normal Dist. Function                         *
 * Will be the function used for quadrature      *
 ************************************************/
double norm_dist(double x, double pi)
{
    /* calc -(x^2)/2 before calling e^x so easier to read */
    double y = -(x*x)/2.0;

    /* calculate pi for the constant, then get a value for the constant */
    /*     note - This IS NOT EFFICIENT; calcs pi every function call */
    /*double pi = fine(0,1,pow(10,-12),0,15);*/
    double c = 1.0 / sq_rt(2.0*pi);

    /* return the value for norm dist at this x value */
    return c*my_exp(y);
}

/*****************************************************
 * Dylan Bollone                                     *
 * Adaptive Quadrature Scheme                        *
 * Using a 2 point gaussian quadrature method        *
 * Adaptive method (splits interval until tol met)   *
 ****************************************************/

double coarse_gauss(double a, double b, double pie, double gauss)
{
    /* this will not be the most efficient because it calculates every time */
    /*double gauss = 1/sq_rt(3);*/

    double h = (b-a)/2.0;
    double c = (b+a)/2.0;
    return h*(norm_dist(-gauss*h + c,pie) + norm_dist(gauss*h + c,pie));
}

double fine_gauss(double a, double b, double tol, double p, double g)
{
    double c = (b+a)/2.0;
    double coarse = coarse_gauss(a,b,p,g);
    double fine = coarse_gauss(a,c,p,g) + coarse_gauss(c,b,p,g);

    if(fabs(fine-coarse)<tol)
    {
	return fine;
    }
    return fine_gauss(a,c,tol/2.0,p,g) + fine_gauss(c,b,tol/2.0,p,g);
}

double cdf(double z)
{
    double a = -15.0;
    double tol = pow(10,-5);

    /* calc gauss node and pi */
    double pi = fine(0,1,pow(10,-12),0,15);
    double node = 1/sq_rt(3);

    return fine_gauss(a,z,tol,pi,node);
}


/************************************************
* Dylan Bollone                                 *
* Main function                                 *
* Where appropriate function calls will be made *
************************************************/
int main()
{
    double z;

    while(1)
    {
        printf("What is the z value? ");
        scanf("%lf",&z);
        printf("Out: %lf\n",cdf(z));
    }
    return 0;
}
