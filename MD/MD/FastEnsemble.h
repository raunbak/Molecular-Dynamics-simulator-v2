//FastEnsemble.h

#ifndef FASTENSEMBLE_H_
#define FASTENSEMBLE_H_

#include <cmath>
#include <cstring>
#include "integrator.h"
#include <iostream>
#include <fstream>
#include "Trap.h"

class FastIon; // dealing with circular dependens same in Ion header.

class FastEnsemble
{
public:
	int NumberOfIons;               // Number of ions in ensemble.
	FastIon *ions;					// pointer to array of ion objects.
	
	long int ***histogram;			// pointer to 3d histogram array of size ?
	long double ***VelHistogram;    // pointer to 3d histogram array of size ?
	long int ***CountHistogram;     // pointer to 3d histogram array of size ?
	double VzSecRMS; // RMS secular z-velocity of ions - used to rescale velocity distribution
	double VrSecRMS; // RMS secular radial-velocity of ions - used to rescale velocity distribution
	double SteadyStateTemperature;  // The temperature simulated.
	double Radius;					// Radius of crystal.
	double Length;					// Length of the crystal
	double PixelToDistance;
	
	int IonOneN;
	int IonTwoN;
	double IonStartVel;
	double IonsTypeOneCharge;
	double IonsTypeTwoCharge;
	double ReducedMass;
	Trap trap;


	int HistNx;
	int HistNy;  // x and y should be the same...
	int HistNz;


	// For storing data in the two ion type simulation (stores all of the information)
	long double ****histograms;

//public:
	// constructor
	FastEnsemble(int m1, int n1,double ioncharge1, int m2, int n2,double ioncharge2, double PixelToD, int Histx, int Histy, int Histz,Trap Thetrap);
	FastEnsemble(int m1, int n1,double ioncharge, double PixelToD, int Histx, int Histy, int Histz,Trap Thetrap);
	// Member functions
	void CrystalGenerator(double Vrf, double Vend); // set ions initial positions in grid and set initial velocities using the plasma model in bcc structure - under construction
	void RescaleVelocityXYZ(double Total_V_x_rms, double Total_V_y_rms,double Total_V_z_rms);
	void CleanUpEnsemble();
	
	// Ememble set and gets
	void SetSteadyStateTemperature(double Val);
	int GetNumberOfIons();
	double Mass(int N); // returning mass of ion N
	double Position(int dim, int N); // returning position
	double Velocity(int dim, int N); // returning velocity
	double getRho0(double Vrf);
	double GetCurrentTemperature();
	double Ekin(); // return total kinetic energy of crystal for given time step
	double Ttot(); // return total temperature of crystal for given time step

	// Histogram and data functions
	void SavePositionToFile();
	void SaveIonDataToFile();
	void SaveTwoIonPositionToFile();

	void MyUpdateVelocityHistogram();
	void InitialiseHistogram();
	void InitialiseVelocityHistogram();
	void InitialiseCountHistogram();
	void UpdateHistogram();
	void UpdateCountHistogram();
	double ReturnHist(int i, int j, int k); //returning value of bin (i,j,k) in histogram
	double ReturnVelHist(int i, int j, int k); //returning value of bin (i,j,k) in histogram
	double ReturnCountHist(int i, int j, int k);


	// For two ion types, they are only for use in the two ion type case!
	void InitialiseHistogramsForTwoIonTypes();
	void UpdateVelandCountHistograms();
	void UpdateHistograms(); 
	double ReturnFromTwoTypeIonHist(int i, int j, int k, int dim);
	void RescaleVelocityXYZTwoSpecies(double Total_V_x_rms1,double Total_V_y_rms1,double Total_V_z_rms1,double Total_V_x_rms2,double Total_V_y_rms2,double Total_V_z_rms2);


	// friend function
	friend double DistanceSquar(FastIon & ion1, FastIon & ion2);
	friend double Distance(FastIon & ion1, FastIon & ion2);
	friend double Ftrap(FastEnsemble & ensemble, int N, int TimeStep, int dim, double Vrf, double Vend, double dt);
	friend double Fcoulumb(FastEnsemble & ensemble, int N, int dim);
	friend double Ffriction(FastEnsemble & ensemble, int N, int dim);
};
#endif

