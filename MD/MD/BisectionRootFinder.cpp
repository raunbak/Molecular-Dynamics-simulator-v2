#include "stdafx.h" // Pre-compiled header.
#include "BisectionRootFinder.h"
#include<iostream>
#include <stdio.h>
#include <math.h>
using namespace std;


void BisectionRootFinder::SetWsquare(double ratio) 
{
	wsquare = ratio;
}

double BisectionRootFinder::asinh(double x) // making invers sine hyperbolic function (not in cmath lib)
{
	return log(x + sqrt(x*x+1));
}

// The AspectRatioEquation
double BisectionRootFinder::f(double alpha) 
{
	if (alpha < 1)
		return  -2*(asinh(sqrt(fabs(pow(alpha,-2) - 1))) - alpha*sqrt(fabs(pow(alpha,-2) - 1)))/(asinh(sqrt(fabs(pow(alpha,-2) - 1))) - sqrt(fabs(pow(alpha,-2) - 1))/alpha) - wsquare;
	else
		return  -2*(asin(sqrt(fabs(pow(alpha,-2) - 1))) - alpha*sqrt(fabs(pow(alpha,-2) - 1)))/(asin(sqrt(fabs(pow(alpha,-2) - 1))) - sqrt(fabs(pow(alpha,-2) - 1))/alpha) - wsquare;
}

double BisectionRootFinder::bisection(double xa, double xb, double eps)
{
double xm, fa, fm;

fa = f(xa);
        while (xb - xa > eps)
{
        xm = (xa + xb) / 2;
        fm = f(xm);
        if (fm == 0.0)
break;
        else if (fa*fm > 0)
{
        xa = xm;
        fa = fm;
}
        else {
        xb = xm;
}
}
        return xm;
}