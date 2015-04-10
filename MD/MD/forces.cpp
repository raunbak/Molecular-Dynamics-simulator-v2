//forces.cpp
#include "stdafx.h" // Pre-compiled header.
#include <cmath>
#include <iostream>
#include "forces.h"
#include "constants.h"


using namespace std;


double Ftot(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend,double dt, Trap trap)
{
	return  Fcoulumb(ensemble, N, dim) + Ftrap(ensemble, N, TimeStep, dim, Vrf, Vend,dt,trap);
}




double Fcoulumb(FastEnsemble & ensemble, int N, int dim)
{
	double sum=0;
	for(int n = 0; n < ensemble.NumberOfIons; n++)
		if (n != N)
			sum += (ensemble.ions[N].Pos[dim]-ensemble.ions[n].Pos[dim])/pow(Distance(ensemble.ions[N], ensemble.ions[n]), 3);
	return e*e/(4*PI*eps0)*sum;
}



double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend, double dt, Trap trap) // fix units, do some thing smart with dt ie make a table
{
	if (dim == 0) // x
		return (e*ensemble.ions[N].charge)*((trap.eta*Vend/pow(trap.z0,2) - Vrf*cos(trap.OmegaRF*TimeStep*dt)/pow(trap.r0,2))*ensemble.ions[N].Pos[dim]);

	if (dim == 1) // y
		return (e*ensemble.ions[N].charge)*(trap.eta*Vend/pow(trap.z0,2) + Vrf/pow(trap.r0,2)*cos(trap.OmegaRF*TimeStep*dt))*ensemble.ions[N].Pos[dim];

	if (dim == 2) // z
		return -2*(e*ensemble.ions[N].charge)*trap.eta*Vend/pow(trap.z0,2)*ensemble.ions[N].Pos[dim]; // Er det 2 gange ladningen ? 

	return NULL;
}


double Ffriction(FastEnsemble & ensemble, int N, int dim)
{
	if (dim == 0 || dim == 1 || dim == 2)
		return -beta*ensemble.ions[N].Vel[dim];
	else
		return 0;
}



double Fpseudo(FastEnsemble & ensemble, int N, int dim, double Vrf, double Vend, Trap trap)
{
	double wz2=2*trap.eta*e*Vend/pow(trap.z0,2)/ensemble.Mass(N); // ERROR in formula from Magnus' thesis! he means potential and not electric potential.
	double wr2=pow(e*Vrf/ensemble.Mass(N)/trap.OmegaRF,2)/2/pow(trap.r0,4)-trap.eta*e*Vend/ensemble.Mass(N)/pow(trap.z0,2);

	if (dim == 0) // x
		return -ensemble.Mass(N)*wr2*ensemble.Position(0,N);

	if (dim == 1) // y
		return -ensemble.Mass(N)*wr2*ensemble.Position(1,N);

	if (dim == 2) // z
		return -ensemble.Mass(N)*wz2*ensemble.Position(2,N);

	return NULL;
}
