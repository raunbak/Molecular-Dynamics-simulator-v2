#ifndef BISECTIONROOTFINDER_H_
#define BISECTIONROOTFINDER_H_


class BisectionRootFinder;

class BisectionRootFinder
{
	double wsquare;
public:
	double bisection(double xa, double xb, double eps);
	double asinh(double x);
	double f (double alpha);
	void SetWsquare(double ratio);
};

#endif