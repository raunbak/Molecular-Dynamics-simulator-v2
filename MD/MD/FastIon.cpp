//FastIon.cpp

#include "stdafx.h" // Pre-compiled header.

#include <cmath>
#include <cstdlib>
#include "FastIon.h"
#include "constants.h"

using namespace std;

// member functions

void FastIon::CleanUpIon()
{
	delete Pos;
	delete Vel;
}

void FastIon::initialize(int mass, double chargeIne) // setting mass and allocating memory
{
	m = mass*u2kg;  // mass should always be in u
	charge = chargeIne;
	Pos = new double [3];
	Vel = new double [3];
}

void FastIon::SetPosition(int dim, double Val) // set position
{
	Pos[dim] = Val;
}

void FastIon::SetVelocity(int dim, double Val) // set position
{
	Vel[dim] = Val;
}

double FastIon::GetMass()
{
	return m; // In kg
}

double FastIon::Position(int dim) // change name to get at some point
{
	return Pos[dim];
}

double FastIon::Velocity(int dim)
{
	return Vel[dim];
}

double FastIon::Velocity()
{
	return sqrt(pow(Vel[0],2) + pow(Vel[1],2) + pow(Vel[2],2));
}

double FastIon::Ekin()
{
	return 0.5*m*pow(pow(Vel[0],2) + pow(Vel[1],2) + pow(Vel[2],2),2);
}

// Friend functions
double DistanceSquar(FastIon & ion1, FastIon & ion2)
{
	return pow(ion1.Pos[0]-ion2.Pos[0],2) + pow(ion1.Pos[1]-ion2.Pos[1],2) + pow(ion1.Pos[2]-ion2.Pos[2],2);
}

// Friend functions
double Distance(FastIon & ion1, FastIon & ion2)
{
	return sqrt(pow(ion1.Pos[0]-ion2.Pos[0],2) + pow(ion1.Pos[1]-ion2.Pos[1],2) + pow(ion1.Pos[2]-ion2.Pos[2],2));
}
