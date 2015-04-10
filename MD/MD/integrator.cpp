// integrator.cpp
#include "stdafx.h" // Pre-compiled header.
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "integrator.h"
#include "forces.h"
#include "constants.h"
#include "cudaforces.cuh"

#include <time.h>

using namespace std;

void TemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap)
{

	double dt = 1/(trap.OmegaRF/2/PI)/StepsPrPeriode; 
	int StartRecordingHistogram = StartHistograms;
	bool PostionSaved = false;
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();
	ensemble.InitialiseCountHistogram();
	
	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ,*Charge;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	Charge = new float [ensemble.GetNumberOfIons()];
	
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	 {
		 Charge[i] = ensemble.ions[i].charge;
	}

	// Declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d, *Charge_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, &Charge_d, ensemble.GetNumberOfIons());


	// Making tempeary position and velocities ( for use in calculation of new pos and vel after calculation on CUDA)
	double **TempPos;
	double **TempVel;

	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	// Array for storing the last position.
	double **PosOneTauAgo;
	PosOneTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosOneTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosOneTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	double **PosTwoTauAgo;
	PosTwoTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosTwoTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosTwoTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	// Total Vrms of the terminal speed
	double Total_V_r_rms = 0.0;
	double Total_V_z_rms = 0.0;
	double Total_V_x_rms = 0.0;
	double Total_V_y_rms = 0.0;
	// 
	ofstream Temperaturefile (trap.dirName+"\\TemperatureData.txt");
	Temperaturefile << "First line of output" << endl ;
	ofstream Periodefile (trap.dirName+"\\PeriodeData.txt");
	Periodefile << "First line of output" << endl ;

	for(int t=1; t <= TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.GetNumberOfIons();N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}
		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d,Charge,Charge_d, ensemble.GetNumberOfIons());
		// calculating coulomb force on GPU with cuda

		//  For output of temperatures.
		double TotalV = 0.0;
		double TotalVz = 0.0;
		double TotalVr = 0.0;
		double TotalVx = 0.0;
		double TotalVy = 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend,dt,trap) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend,dt,trap) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend,dt,trap) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz += pow(TempVel[N][2],2);
			TotalVr += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
			TotalVx += pow(TempVel[N][0],2);
			TotalVy += pow(TempVel[N][1],2);
		}

		TotalV  = sqrt(TotalV / ensemble.GetNumberOfIons()); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVz = sqrt(TotalVz / ensemble.GetNumberOfIons());
		TotalVr = sqrt(TotalVr / ensemble.GetNumberOfIons());
		TotalVx = sqrt(TotalVx / ensemble.GetNumberOfIons());
		TotalVy = sqrt(TotalVy / ensemble.GetNumberOfIons());
		
		// Saving to temperaturefile.
		Temperaturefile << t <<  ",  "
						<< (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVz,2)* ensemble.ions[0].m)/(Kb)	<<  ",  "
						<< (pow(TotalVx,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVy,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< endl;

		// saving time step
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}
		

		// This is when 2*PI/omegaRF has passed. Which is when we calculate termalspeed.
		/*
		if ((t) % ((int) (StepsPrPeriode + 28))) == 0)
		{
			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau;
				double TermalVy	= (PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosOneTauAgo[N][2] - ensemble.ions[N].Pos[2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}
		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);
		}
		*/

		double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);
		//if (((int)(t+ LowestPartOfPeriode) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0) 
		//if (((t) % ((int) (LowestPartOfPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)
	    if (((int)(t+ LowestPartOfPeriode-1) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)  // for 105	
		{
			// Storing the postion, for next calculation.
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosTwoTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}


			for(int N = 0; N < ensemble.GetNumberOfIons(); N++)
			{
				if (Periodefile.is_open())
				{
				Periodefile << t								<< ",  " <<
								N + 1							<< ",  " <<
								ensemble.ions[N].GetMass()		<< ",  " <<
								ensemble.ions[N].Pos[0]			<< ",  " <<
								ensemble.ions[N].Pos[1]			<< ",  " <<
								ensemble.ions[N].Pos[2]			<< ",  " <<
								ensemble.ions[N].Vel[0]			<< ",  " <<
								ensemble.ions[N].Vel[1]			<< ",  " <<
								ensemble.ions[N].Vel[2]			<< ",  " <<
								ensemble.ions[N].Velocity()		<< ",  " <<
								(pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<< ",  " <<
								endl;
				}   
			}
			}

		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{

			/*
			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau;
				double TermalVy	= (PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosOneTauAgo[N][2] - ensemble.ions[N].Pos[2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosOneTauAgo[N][0] - ensemble.ions[N].Pos[0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosOneTauAgo[N][1] - ensemble.ions[N].Pos[1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}
			*/

			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / trap.OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau;
				double TermalVy	= (PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosTwoTauAgo[N][2] - PosOneTauAgo[N][2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_x_rms = Total_V_x_rms / ensemble.GetNumberOfIons();
			Total_V_x_rms = sqrt(Total_V_x_rms);

			Total_V_y_rms = Total_V_y_rms / ensemble.GetNumberOfIons();
			Total_V_y_rms = sqrt(Total_V_y_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = PosTwoTauAgo[N][dim];
				}
			}

	

			// Rescale velocities.
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);

		}	
		// update 3d histogram
		
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistogram();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}

	}

	Periodefile.close();
	Temperaturefile.close();

	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d,Charge_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	delete Charge;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;
}

void TwoIonTypesTemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap)
{

	double dt = 1/(trap.OmegaRF/2/PI)/StepsPrPeriode; 
	int StartRecordingHistogram = StartHistograms;
	bool PostionSaved = false;
	int numberOfRescales = 0;
	// allocate memory to the histogram
	ensemble.InitialiseHistogramsForTwoIonTypes();
	
	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ, *Charge;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	Charge = new float [ensemble.GetNumberOfIons()];
	
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	 {
		 Charge[i] = ensemble.ions[i].charge;
	}
 


	// Declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d, *Charge_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d,&Charge_d, ensemble.GetNumberOfIons());


	// Making tempeary position and velocities ( for use in calculation of new pos and vel after calculation on CUDA)
	double **TempPos;
	double **TempVel;

	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	// Array for storing the last position.
	double **PosOneTauAgo;
	PosOneTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosOneTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosOneTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	double **PosTwoTauAgo;
	PosTwoTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosTwoTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosTwoTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	// Total Vrms of the terminal speed
	double Total_V_r_rms = 0.0;
	double Total_V_z_rms = 0.0;
	double Total_V_x_rms = 0.0;
	double Total_V_y_rms = 0.0;
	// 
	ofstream TemperaturefileOne (trap.dirName+"\\TemperatureDataIonTypeOne.txt");
	ofstream TemperaturefileTwo (trap.dirName+"\\TemperatureDataIonTypeTwo.txt");
	ofstream TemperaturefileBoth (trap.dirName+"\\TemperatureDataIons.txt");
	TemperaturefileOne << "First line of output" << endl;
	TemperaturefileTwo << "First line of output" << endl;
	TemperaturefileBoth << "First line of output" << endl;
	
	
	
	//ofstream Temperaturefile (trap.dirName+"\\TemperatureData.txt");

	//Temperaturefile << "First line of output" << endl ;
	ofstream Periodefile (trap.dirName+"\\PeriodeData.txt");
	Periodefile << "First line of output" << endl ;

	for(int t=1; t <= TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.GetNumberOfIons();N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}
		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d,Charge,Charge_d ,ensemble.GetNumberOfIons());
		// calculating coulomb force on GPU with cuda
		//  For output of temperatures.

		double TotalVOne	= 0.0;
		double TotalVzOne	= 0.0;
		double TotalVrOne	= 0.0;
		double TotalVxOne	= 0.0;
		double TotalVyOne	= 0.0;

		double TotalVTwo	= 0.0;
		double TotalVzTwo	= 0.0;
		double TotalVrTwo	= 0.0;
		double TotalVxTwo	= 0.0;
		double TotalVyTwo	= 0.0;

		double TotalV		= 0.0;
		double TotalVz		= 0.0;
		double TotalVr		= 0.0;
		double TotalVx		= 0.0;
		double TotalVy		= 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend,dt,trap) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend,dt,trap) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend,dt,trap) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz += pow(TempVel[N][2],2);
			TotalVr += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
			TotalVx += pow(TempVel[N][0],2);
			TotalVy += pow(TempVel[N][1],2);
			if (N == (ensemble.IonOneN)-1) 
			{
				TotalVOne	= sqrt( TotalV / ensemble.IonOneN);
				TotalVzOne	= sqrt( TotalVz / ensemble.IonOneN);
				TotalVrOne	= sqrt( TotalVr / ensemble.IonOneN);
				TotalVxOne	= sqrt (TotalVx / ensemble.IonOneN);
				TotalVyOne	= sqrt( TotalVy / ensemble.IonOneN);

			}
			if (N == (ensemble.IonOneN))
			{
				TotalVTwo	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
				TotalVzTwo  += pow(TempVel[N][2],2);
				TotalVrTwo  += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
				TotalVxTwo  += pow(TempVel[N][0],2);
				TotalVyTwo  += pow(TempVel[N][1],2);	
			}
		}

		TotalV  = sqrt(TotalV / ensemble.GetNumberOfIons()); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVz = sqrt(TotalVz / ensemble.GetNumberOfIons());
		TotalVr = sqrt(TotalVr / ensemble.GetNumberOfIons());
		TotalVx = sqrt(TotalVx / ensemble.GetNumberOfIons());
		TotalVy = sqrt(TotalVy / ensemble.GetNumberOfIons());

		TotalVTwo  = sqrt(TotalVTwo  / ensemble.IonTwoN); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVzTwo = sqrt(TotalVzTwo / ensemble.IonTwoN);
		TotalVrTwo = sqrt(TotalVrTwo / ensemble.IonTwoN);
		TotalVxTwo = sqrt(TotalVxTwo / ensemble.IonTwoN);
		TotalVyTwo = sqrt(TotalVyTwo / ensemble.IonTwoN);




/*
		//  For output of temperatures.
		double TotalV = 0.0;
		double TotalVz = 0.0;
		double TotalVr = 0.0;
		double TotalVx = 0.0;
		double TotalVy = 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t-1, 0, Vrf, Vend,dt,trap) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t-1, 1, Vrf, Vend,dt,trap) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t-1, 2, Vrf, Vend,dt,trap) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz += pow(TempVel[N][2],2);
			TotalVr += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
			TotalVx += pow(TempVel[N][0],2);
			TotalVy += pow(TempVel[N][1],2);
		}

		TotalV  = sqrt(TotalV / ensemble.GetNumberOfIons()); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVz = sqrt(TotalVz / ensemble.GetNumberOfIons());
		TotalVr = sqrt(TotalVr / ensemble.GetNumberOfIons());
		TotalVx = sqrt(TotalVx / ensemble.GetNumberOfIons());
		TotalVy = sqrt(TotalVy / ensemble.GetNumberOfIons());
*/		
		// Saving to temperaturefile.

		TemperaturefileBoth << t <<  ",  "
						<< (pow(TotalV,2) * ensemble.ReducedMass)/(3*Kb)	<<  ",  "
						<< (pow(TotalVz,2)* ensemble.ReducedMass)/(Kb)	<<  ",  "
						<< (pow(TotalVx,2)* ensemble.ReducedMass)/(Kb)	<<  ",  " 
						<< (pow(TotalVy,2)* ensemble.ReducedMass)/(Kb)	<<  ",  " 
						<< endl;

		TemperaturefileOne << t <<  ",  "
						<< (pow(TotalVOne,2) * ensemble.ions[ensemble.IonOneN-1].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVzOne,2)* ensemble.ions[ensemble.IonOneN-1].m)/(Kb)	<<  ",  "
						<< (pow(TotalVxOne,2)* ensemble.ions[ensemble.IonOneN-1].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVyOne,2)* ensemble.ions[ensemble.IonOneN-1].m)/(Kb)	<<  ",  " 
						<< endl;

		TemperaturefileTwo << t <<  ",  "
						<< (pow(TotalVTwo,2) * ensemble.ions[ensemble.IonTwoN-1].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVzTwo,2)* ensemble.ions[ensemble.IonTwoN-1].m)/(Kb)	<<  ",  "
						<< (pow(TotalVxTwo,2)* ensemble.ions[ensemble.IonTwoN-1].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVyTwo,2)* ensemble.ions[ensemble.IonTwoN-1].m)/(Kb)	<<  ",  " 
						<< endl;

		/*Temperaturefile << t <<  ",  "
						<< (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVz,2)* ensemble.ions[0].m)/(Kb)	<<  ",  "
						<< (pow(TotalVx,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVy,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< endl;
		*/
		// saving time step
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}
		//-----------------------------------------------------------------------------------------------------------
		if (t == (53+105*numberOfRescales)-1)
		{
			numberOfRescales++;

			if (t >= StartRecordingHistogram)
			{
				ensemble.UpdateVelandCountHistograms();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" <<  ((pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) ) << endl;
			cout << "Tone : " << (pow(TotalVOne,2) * ensemble.ions[ensemble.IonOneN-1].m)/(3*Kb) << endl;
			cout << "Ttwo : " << (pow(TotalVTwo,2) * ensemble.ions[ensemble.IonTwoN-1].m)/(3*Kb)<< endl;
			//ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);
			ensemble.RescaleVelocityXYZTwoSpecies(TotalVxOne,TotalVyOne,TotalVzOne,TotalVxTwo,TotalVyTwo,TotalVzTwo);
			

		}
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistograms();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}
	}


	/*
		//------------------------------------------------------------------------------------------------------------
		//double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);

	    //if (((int)(t+ LowestPartOfPeriode-1) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)  // for 105	
		if (t == (53+105*numberOfRescales)-1)
		{
			numberOfRescales++;
			// Storing the postion, for next calculation.
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosTwoTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}

			if (t >= StartRecordingHistogram)
			{
				ensemble.UpdateVelandCountHistograms();
				//ensemble.MyUpdateVelocityHistogram();
				//ensemble.UpdateCountHistogram();
				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			


			for(int N = 0; N < ensemble.GetNumberOfIons(); N++)
			{
				if (Periodefile.is_open())
				{
				Periodefile << t								<< ",  " <<
								N + 1							<< ",  " <<
								ensemble.ions[N].GetMass()		<< ",  " <<
								ensemble.ions[N].Pos[0]			<< ",  " <<
								ensemble.ions[N].Pos[1]			<< ",  " <<
								ensemble.ions[N].Pos[2]			<< ",  " <<
								ensemble.ions[N].Vel[0]			<< ",  " <<
								ensemble.ions[N].Vel[1]			<< ",  " <<
								ensemble.ions[N].Vel[2]			<< ",  " <<
								ensemble.ions[N].Velocity()		<< ",  " <<
								(pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<< ",  " <<
								endl;
				}   
			}
			}

		}


		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{

			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / trap.OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau;
				double TermalVy	= (PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosTwoTauAgo[N][2] - PosOneTauAgo[N][2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_x_rms = Total_V_x_rms / ensemble.GetNumberOfIons();
			Total_V_x_rms = sqrt(Total_V_x_rms);

			Total_V_y_rms = Total_V_y_rms / ensemble.GetNumberOfIons();
			Total_V_y_rms = sqrt(Total_V_y_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = PosTwoTauAgo[N][dim];
				}
			}

	

			// Rescale velocities.
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);

		}	
		// update 3d histogram
		
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistograms();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}	
	}
	*/

	Periodefile.close();
	TemperaturefileOne.close();
	TemperaturefileTwo.close();
	TemperaturefileBoth.close();

	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d,Charge_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	delete Charge;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;
}



void NewTemperatureRescaleLeFrogintegrator(FastEnsemble & ensemble, int TimeSteps, int StartHistograms, double Vrf, double Vend, double StepsPrPeriode, Trap & trap)
	{

	int RescaleTime = 102;
	int RescaleOffset[10] = {1,0,0,-1,-1,-1,0,0,1,1};
	int OffsetIndex = 0;

	int const number = 20;
	int tlist [number] = {102, 204,305,406,506,606,706,807,908,1010,1112,1214,1315,1416,1516,1616,1716,1817,1918,2020};
	double dt = 1/(trap.OmegaRF/2/PI)/StepsPrPeriode; 
	int StartRecordingHistogram = StartHistograms;
	bool PostionSaved = false;
	// allocate memory to the histogram
	ensemble.InitialiseHistogram();
	ensemble.InitialiseVelocityHistogram();
	ensemble.InitialiseCountHistogram();
	
	// make Cuda Pos and Force arrays - its float because we dont want to overflow the GPU
	float *PosX, *PosY, *PosZ, *ForceX, *ForceY, *ForceZ,*Charge;

	PosX = new float [ensemble.GetNumberOfIons()];
	PosY = new float [ensemble.GetNumberOfIons()];
	PosZ = new float [ensemble.GetNumberOfIons()];
	ForceX = new float [ensemble.GetNumberOfIons()];
	ForceY = new float [ensemble.GetNumberOfIons()];
	ForceZ = new float [ensemble.GetNumberOfIons()];

	Charge = new float [ensemble.GetNumberOfIons()];
	
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	 {
		 Charge[i] = ensemble.ions[i].charge;
	}

	// Declare and allocate memory on GPU
	float *PosX_d, *PosY_d, *PosZ_d, *ForceX_d, *ForceY_d, *ForceZ_d, *Charge_d;
	CudaCoulombAlloc(&PosX_d, &PosY_d, &PosZ_d, &ForceX_d, &ForceY_d, &ForceZ_d, &Charge_d, ensemble.GetNumberOfIons());


	// Making tempeary position and velocities ( for use in calculation of new pos and vel after calculation on CUDA)
	double **TempPos;
	double **TempVel;

	TempPos = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempPos[i] = new double [3];

	TempVel = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		TempVel[i] = new double [3];


	// Array for storing the last position.
	double **PosOneTauAgo;
	PosOneTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosOneTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosOneTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	double **PosTwoTauAgo;
	PosTwoTauAgo = new double *[ensemble.GetNumberOfIons()];
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
	{
		PosTwoTauAgo[i] = new double [3];
		for(int dim=0; dim < 3; dim++)
			{
				PosTwoTauAgo[i][dim] = ensemble.ions[i].Pos[dim];
			}
	}

	// Total Vrms of the terminal speed
	double Total_V_r_rms = 0.0;
	double Total_V_z_rms = 0.0;
	double Total_V_x_rms = 0.0;
	double Total_V_y_rms = 0.0;
	// 
	ofstream Temperaturefile (trap.dirName+"\\TemperatureData.txt");
	Temperaturefile << "First line of output" << endl ;
	//ofstream Periodefile (trap.dirName+"\\PeriodeData.txt");
	//Periodefile << "First line of output" << endl ;

	int offset = 0;
	int numberOfRescales = 0;

	for(int t=0; t < TimeSteps ; t++)
	{

		// Prepare Data for CUDA transfer
		for(int N = 0; N < ensemble.GetNumberOfIons();N++)
		{
			PosX[N] = ((float) ensemble.ions[N].Pos[0]);
			PosY[N] = ((float) ensemble.ions[N].Pos[1]);
			PosZ[N] = ((float) ensemble.ions[N].Pos[2]);

		}
		// updating position, looping through ions and x-y-z
		//CoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ensemble.GetNumberOfIons());
		FastCoulombWrapper( ForceX, ForceY, ForceZ, PosX, PosY, PosZ, ForceX_d, ForceY_d, ForceZ_d, PosX_d, PosY_d, PosZ_d,Charge,Charge_d, ensemble.GetNumberOfIons());
		// calculating coulomb force on GPU with cuda

		//  For output of temperatures.
		double TotalV = 0.0;
		double TotalVz = 0.0;
		double TotalVr = 0.0;
		double TotalVx = 0.0;
		double TotalVy = 0.0;

		// updating position, looping through ions and x-y-z
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			// This is dim X
			double Ftot = Ftrap(ensemble, N, t, 0, Vrf, Vend,dt,trap) + ((double) ForceX[N]);
			TempPos[N][0]  = ensemble.ions[N].Pos[0] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[0];
			TempVel[N][0]  = ensemble.ions[N].Vel[0] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Y
			Ftot = Ftrap(ensemble, N, t, 1, Vrf, Vend,dt,trap) + ((double) ForceY[N]);
				
			TempPos[N][1]  = ensemble.ions[N].Pos[1] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[1];
			TempVel[N][1]  = ensemble.ions[N].Vel[1] + dt/ensemble.ions[N].m*Ftot;

			// This is dim Z
			Ftot = Ftrap(ensemble, N, t, 2, Vrf, Vend,dt,trap) + ((double) ForceZ[N]);
			
			TempPos[N][2]  = ensemble.ions[N].Pos[2] + dt*dt*Ftot/ensemble.ions[N].m/2 + dt*ensemble.ions[N].Vel[2];
			TempVel[N][2]  = ensemble.ions[N].Vel[2] + dt/ensemble.ions[N].m*Ftot;

			TotalV	+= pow(TempVel[N][0],2) + pow(TempVel[N][1],2) + pow(TempVel[N][2],2);
			TotalVz += pow(TempVel[N][2],2);
			TotalVr += pow(TempVel[N][0],2) + pow(TempVel[N][1],2);
			TotalVx += pow(TempVel[N][0],2);
			TotalVy += pow(TempVel[N][1],2);
		}

		TotalV  = sqrt(TotalV / ensemble.GetNumberOfIons()); // Now TotalV is the root mean squared of the total speed. (So avg speed of 1 ion)
		TotalVz = sqrt(TotalVz / ensemble.GetNumberOfIons());
		TotalVr = sqrt(TotalVr / ensemble.GetNumberOfIons());
		TotalVx = sqrt(TotalVx / ensemble.GetNumberOfIons());
		TotalVy = sqrt(TotalVy / ensemble.GetNumberOfIons());
		
		// Saving to temperaturefile.
		Temperaturefile << t <<  ",  "
						<< (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<<  ",  "
						<< (pow(TotalVz,2)* ensemble.ions[0].m)/(Kb)	<<  ",  "
						<< (pow(TotalVx,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< (pow(TotalVy,2)* ensemble.ions[0].m)/(Kb)	<<  ",  " 
						<< endl;

		// saving time step
		for(int N=0; N < ensemble.GetNumberOfIons(); N++)
		{
			for(int dim=0; dim < 3; dim++)
			{
				ensemble.ions[N].Pos[dim] = TempPos[N][dim];
				ensemble.ions[N].Vel[dim] = TempVel[N][dim];
			}
		}


		//double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);
		//if (((int)(t+ LowestPartOfPeriode) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0) 
		//if (((t) % ((int) (LowestPartOfPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)
	    //if (((int)(t+ LowestPartOfPeriode-1) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)  // for 105	


		if (t == (53+105*numberOfRescales)-1)
		{
			numberOfRescales++;

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" <<  ((pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) ) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);

		

		}


		/*
		if (t==RescaleTime)//(cos(t*(2*PI)/StepsPrPeriode) == 1)
		{
			cout << "cos^2 == 1 : "<< t << endl;

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" << (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);

			RescaleTime = RescaleTime + StepsPrPeriode + RescaleOffset[OffsetIndex];
			OffsetIndex++;

			if (OffsetIndex == 10)
			{
				OffsetIndex = 0;
			}



		}
		*/
		/*
		int *foo = std::find(tlist, tlist+number, t);

		if (foo != tlist+number)//(cos(t*(2*PI)/StepsPrPeriode) == 1)
		{
			cout << "cos^2 == 1 : "<< t << endl;

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" << (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);

		}
		
		if ((t+offset)%102 == 0 && t != 0 && t < 450)//(cos(t*(2*PI)/StepsPrPeriode) == 1)
		{
			numberOfRescales = numberOfRescales +1;
			cout << "cos^2 == 1 : "<< t << endl;

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" << (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);

			if (numberOfRescales == 2) 
			{
				offset = 1;
			}
			if (numberOfRescales == 3) 
			{
				offset = 2;
			}


		}

		if ((t-6)%100 == 0 && t>450 && t < 750)
		{
			cout << "cos^2 == 1 : "<< t << endl;
			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" << (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);
		}

		if ((t-7)%100 == 0 && t>750 )
		{
			cout << "cos^2 == 1 : "<< t << endl;
			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}
			}

			cout << "T :" << (pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb) << endl;
			ensemble.RescaleVelocityXYZ(TotalVx,TotalVy,TotalVz);
		}

		*/
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistogram();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}


		//if (cos(t*(2*PI)/StepsPrPeriode) == -1)
		//{
			//cout << "cos^2 == -1 : "<< t << endl;


		//}
	}


		/*

		double LowestPartOfPeriode = floor((StepsPrPeriode/2) + 0.5);
		//if (((int)(t+ LowestPartOfPeriode) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0) 
		//if (((t) % ((int) (LowestPartOfPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)
	    if (((int)(t+ LowestPartOfPeriode-1) % ((int) (StepsPrPeriode))) == 0 && ((t) % ((int) (StepsPrPeriode))) != 0)  // for 105	
		{
			// Storing the postion, for next calculation.
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosTwoTauAgo[N][dim] = ensemble.ions[N].Pos[dim];
				}
			}

			if (t >= StartRecordingHistogram)
			{
				ensemble.MyUpdateVelocityHistogram();
				ensemble.UpdateCountHistogram();

				if (PostionSaved == false)
				{
				ensemble.SavePositionToFile();
				PostionSaved = true;
				}


			for(int N = 0; N < ensemble.GetNumberOfIons(); N++)
			{
				if (Periodefile.is_open())
				{
				Periodefile << t								<< ",  " <<
								N + 1							<< ",  " <<
								ensemble.ions[N].GetMass()		<< ",  " <<
								ensemble.ions[N].Pos[0]			<< ",  " <<
								ensemble.ions[N].Pos[1]			<< ",  " <<
								ensemble.ions[N].Pos[2]			<< ",  " <<
								ensemble.ions[N].Vel[0]			<< ",  " <<
								ensemble.ions[N].Vel[1]			<< ",  " <<
								ensemble.ions[N].Vel[2]			<< ",  " <<
								ensemble.ions[N].Velocity()		<< ",  " <<
								(pow(TotalV,2) * ensemble.ions[0].m)/(3*Kb)	<< ",  " <<
								endl;
				}   
			}
			}

		}

		if (((t) % ((int) (StepsPrPeriode))) == 0)
		{

			// Calculating the root mean square of termal vel's.
			double tau = (2*PI) / trap.OmegaRF; 

			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				double TermalVx = (PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau;
				double TermalVy	= (PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau;
				
				double TermalVr	= sqrt( pow(TermalVx,2) + pow(TermalVy,2));
				Total_V_r_rms	+= pow(TermalVr,2);

				Total_V_z_rms	+= pow((PosTwoTauAgo[N][2] - PosOneTauAgo[N][2] ) / tau,2); 
				Total_V_x_rms	+= pow((PosTwoTauAgo[N][0] - PosOneTauAgo[N][0] ) / tau,2); 
				Total_V_y_rms	+= pow((PosTwoTauAgo[N][1] - PosOneTauAgo[N][1] ) / tau,2); 
			}

			Total_V_z_rms = Total_V_z_rms / ensemble.GetNumberOfIons();
			Total_V_z_rms = sqrt(Total_V_z_rms);

			Total_V_x_rms = Total_V_x_rms / ensemble.GetNumberOfIons();
			Total_V_x_rms = sqrt(Total_V_x_rms);

			Total_V_y_rms = Total_V_y_rms / ensemble.GetNumberOfIons();
			Total_V_y_rms = sqrt(Total_V_y_rms);

			Total_V_r_rms = Total_V_r_rms / ensemble.GetNumberOfIons();
			Total_V_r_rms = sqrt(Total_V_r_rms);

			// Storing the postion, for next calculation
			for(int N=0; N < ensemble.GetNumberOfIons(); N++)
			{
				
				for(int dim=0; dim < 3; dim++)
				{
					PosOneTauAgo[N][dim] = PosTwoTauAgo[N][dim];
				}
			}

	

			// Rescale velocities.
			ensemble.RescaleVelocityXYZ(Total_V_x_rms,Total_V_y_rms,Total_V_z_rms);

		}	
		// update 3d histogram
		
		if (t >= StartRecordingHistogram)
		{
			if(t == StartRecordingHistogram)
			{
				cout << "Record started " << endl; 
			}
			ensemble.UpdateHistogram();
		}
		
		if ( t % 10000 == 0 ) {
			cout << "****** t = "<< t<< " *****"<< endl;
		}

	}

	*/
	//Periodefile.close();
	Temperaturefile.close();

	// cleaning up on GPU
	CudaCoulombFree(PosX_d, PosY_d, PosZ_d, ForceX_d, ForceY_d, ForceZ_d,Charge_d);
	// cleaning up in general
	delete PosX;
	delete PosY;
	delete PosZ;
	delete ForceX;
	delete ForceY;
	delete ForceZ;
	delete Charge;
	// delete SecVel; remember to fix this
	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempPos[i];
	delete TempPos;

	for(int i = 0; i < ensemble.GetNumberOfIons();i++)
		delete [] TempVel[i];
	delete TempVel;
}


