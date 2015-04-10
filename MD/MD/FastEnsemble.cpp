#include "stdafx.h" // Pre-compiled header.
#include "FastEnsemble.h"
#include "FastIon.h"
#include "constants.h"
#include <iostream>

#include <stdarg.h>


using namespace std;

// Constructor for one ion.
FastEnsemble::FastEnsemble(int m1, int n1,double ioncharge1, int m2, int n2,double ioncharge2, double PixelToD, int Histx, int Histy, int Histz,Trap Thetrap)
{
	trap = Thetrap;
	HistNx = Histx;
	HistNy = Histy;
	HistNz = Histz;
	PixelToDistance = PixelToD;
	SteadyStateTemperature = 0.0;//SteadyStateTzSec;
	NumberOfIons=n1+n2;
	ions = new FastIon [NumberOfIons];
	ReducedMass = ((m1*m2)/(m1+m2))*u2kg;
	IonOneN = n1;
	IonTwoN = n2;

	for(int n = 0; n < NumberOfIons; n++) {
		if (n<n1 ){
		
  		ions[n].initialize(m1,ioncharge1); // setting mass, secular velocity and allocating memory for the ions pos and vel arrays
		}
		else {
		
  		ions[n].initialize(m2,ioncharge2); // setting mass, secular velocity and allocating memory for the ions pos and vel arrays
		}
	}
}

FastEnsemble::FastEnsemble(int m1, int n1,double ioncharge, double PixelToD, int Histx, int Histy, int Histz,Trap Thetrap)
{
	trap = Thetrap;
	HistNx = Histx;
	HistNy = Histy;
	HistNz = Histz;
	PixelToDistance = PixelToD;
	SteadyStateTemperature = 0.0;//SteadyStateTzSec;
	NumberOfIons=n1;
	ReducedMass = m1*u2kg;
	ions = new FastIon [NumberOfIons];

	for(int n = 0; n < NumberOfIons; n++) {
		
		ions[n].initialize(m1,ioncharge); // setting mass, secular velocity and allocating memory for the ions pos and vel arrays
	}
}


	
// Member functions
void FastEnsemble::CleanUpEnsemble()
{
	for(int n = 0; n  < NumberOfIons; n++)
		ions[n].CleanUpIon();


	for(int i = 0; i < HistNx;i++)
	{
		for(int j = 0; j < HistNy;j++)
		{
			delete [] histogram[i][j];
			delete [] VelHistogram[i][j];
		}

		delete [] histogram[i];
		delete [] VelHistogram[i];
	}
	delete histogram;
	delete VelHistogram;

}

double FastEnsemble::GetCurrentTemperature()
{
	return (pow(VzSecRMS,2)+pow(VrSecRMS,2))*Mass(0)/Kb/3; // REMEBER TO CHANGE THIS FOR MASSES
}

void FastEnsemble::RescaleVelocityXYZ(double Total_V_x_rms,double Total_V_y_rms,double Total_V_z_rms)
{
		// calculate current temperature
	double CurrentTempZ = (pow(Total_V_z_rms,2)*Mass(0))/(Kb);
	double CurrentTempX = (pow(Total_V_x_rms,2)*Mass(0))/(Kb);
	double CurrentTempY = (pow(Total_V_y_rms,2)*Mass(0))/(Kb);

	// rescale velocity distribution (Maybe change which fraction of SST is given)
	 
	//cout << "T :" << (pow(TotalV,2) * Mass(0))/(3*Kb) <<		sqrt(pow(CurrentTempZ,2)+pow(CurrentTempX,2) + pow(CurrentTempY,2) )<< endl;
	double T_fraction = SteadyStateTemperature;//(sqrt(3)*SteadyStateTemperature)/3; // This was just SteadyState Before

	double az = sqrt((T_fraction)/CurrentTempZ);
	double ax = sqrt((T_fraction)/CurrentTempX);
	double ay = sqrt((T_fraction)/CurrentTempY);

	cout << "T's :" << CurrentTempX << " "<< CurrentTempY << " " << CurrentTempZ << endl;
	cout << ax <<" " << ay << " " << az << endl; 

	double lowerLimit = 0.98;//sqrt(100/102.5); // Hvornår der ikke skal køles mere.
	double upperLimit = 1.02;//sqrt(100/97.5); // Hvornår der ikke skal varmes op mere.

	//cout << lowerLimit << " " << upperLimit << endl;
	// Good value for raise temp is 1.005

	/*
	if ( ax < 0.995 ) 
	{
		ax = 0.9990;
	}
	else if (ax > 1.005 )
	{
		ax = 1.0005;
	}
	else 
	{
		ax = 1.0;
	}
	if ( ay < 0.995 ) 
	{
		ay = 0.9990;
	}
	else if (ay > 1.005 )
	{
		ay = 1.0005;
	}
	else 
	{
		ay = 1.0;
	}
	if ( az < 0.995 ) 
	{
		az = 0.9990;
	}
	else if (az > 1.005 )
	{
		az = 1.0005;
	}
	else 
	{
		az = 1.0;
	}
	*/


		if (ax > upperLimit)
		{
			//ax = 1.0000002;
			//ax = 1.000002; // We use this!
			ax = 1.000002;

		} 
		else if (ax < lowerLimit)
		{
			ax = 0.996;


		}
		else {
			ax = 1.0;
		}

		if (ay > upperLimit)
		{
			//ay = 1.0000002;
			//ay = 1.000002; // We use this
			ay = 1.000002;//1.00002; //Testet med lidt højere opkøling.
		}
		else if (ay < lowerLimit)
		{
			ay = 0.996;
		}
		else {

			ay = 1.0;
		}

		if (az > upperLimit)
		{
			//az = 1.0000002;
			//az = 1.000002; // We use this!
			az = 1.000002;
		}
		else if (az < lowerLimit)
		{
			az = 0.996;
		}
		else {
			az = 1.0;
		}
		

		/*
		if (ax < lowerLimit)
		{
			ax = 0.9990;
		}
		else {

			if ( ReachedTempArea == true)
			{
			ax = 0.9999;
			ReachedTempArea = false;
			}
			else
			{
			ax = 1.0001;
			ReachedTempArea = true;
			}
		}

		if (ay < lowerLimit)
		{
			ay = 0.9990;
		}
		else {

			if ( ReachedTempArea == true)
			{
			ay = 0.9999;
			ReachedTempArea = false;
			}
			else
			{
			ay = 1.0001;
			ReachedTempArea = true;
			}
		}

		if (az < lowerLimit)
		{
			az = 0.9990;
		}
		else {
			if ( ReachedTempArea == true)
			{
			az = 0.99999;
			ReachedTempArea = false;
			}
			else
			{
			az = 1.00001;
			ReachedTempArea = true;
			}
		}
		
		if (ax < lowerLimit)
		{
			ax = 0.9990;
		}
		else {

			ax = 1.0000;
		}

		if (ay < lowerLimit)
		{
			ay = 0.9990;
		}
		else {

	       ay = 1.0000;
		}

		if (az < lowerLimit)
		{
			az = 0.9990;
		}
		else {
			az = 1.0000;
		}
		*/

	for(int n = 0; n < NumberOfIons; n++)
	{

		
		//SquareRoot(xd*xd + yd*yd + zd*zd)
		double RadiusFromCenter = sqrt(pow(ions[n].Position(1),2) +pow(ions[n].Position(2),2) + pow(ions[n].Position(3),2));
		//cout << "Radius" << RadiusFromCenter << endl;
		//cout << "radiustimes1000  " <<(RadiusFromCenter)*1000 << endl;
		//cout <<" Old: " <<0.999 - (RadiusFromCenter)*500 << endl;
		//cout << " Current :" << 0.999-(RadiusFromCenter*100) << endl;
		//cout << 0.999-(RadiusFromCenter*75) << endl;
		if (ax == 0.996)
		{
			//ax = 0.991 - (RadiusFromCenter)*1000; // Used 0.991 before and 1000
			  //ax = 0.996 - (RadiusFromCenter/Radius)*0.006;
			  //ax = 0.996 - exp(-Radius/RadiusFromCenter)*0.006;
			  ax = 0.999-(RadiusFromCenter*100);
			  //ax = 0.999 - (RadiusFromCenter)*250;
		}

		if (ay == 0.996)
		{
			//ay = 0.991 - (RadiusFromCenter)*1000;
			//ay = 0.996 - (RadiusFromCenter/Radius)*0.006;
			//ay =  0.996 - exp(-Radius/RadiusFromCenter)*0.006;
			ay = 0.999 - (RadiusFromCenter)*100;
		}

		if (az == 0.996)
		{
			//az = 0.991 - (RadiusFromCenter)*1000;
			//az = 0.996 - (RadiusFromCenter/Radius)*0.006;
			//az = 0.996 - exp(-Radius/RadiusFromCenter)*0.006;
			az = 0.999 - (RadiusFromCenter)*100;
		}


		if (ax == 1.000002)
		{
			//ax = 1.0002 - (RadiusFromCenter); // Used 0.991 before and 1000
			//ax = 0.998 - (RadiusFromCenter/Radius)*0.04;
			//ax = 0.999 - exp(-Radius/RadiusFromCenter)*0.04;
			//ax = 1.000002 - (RadiusFromCenter/Radius)*0.000001;
			//ax = 1.000002 - exp(-Radius/RadiusFromCenter)*0.000001;
			ax = 1.000002 -(RadiusFromCenter/100);
		}

		if (ay == 1.000002)
		{
			//ay = 1.0002 - (RadiusFromCenter);
			//ay = 0.998 - (RadiusFromCenter/Radius)*0.04;
			//ay = 0.999 - exp(-Radius/RadiusFromCenter)*0.04;
			//ay = 1.000002 - (RadiusFromCenter/Radius)*0.000001;
			//ay = 1.000002 - exp(-Radius/RadiusFromCenter)*0.000001;
			ay = 1.000002 -(RadiusFromCenter/100);
		}

		if (az == 1.000002)
		{
			//az = 1.0002 - (RadiusFromCenter);
			//az = 0.998 - (RadiusFromCenter/Radius)*0.04;
			//az = 0.999 - exp(-Radius/RadiusFromCenter)*0.04;
			//az = 1.000002 - (RadiusFromCenter/Radius)*0.000001;
			//az = 1.000002 - exp(-Radius/RadiusFromCenter)*0.000001;
			az = 1.000002 -(RadiusFromCenter/100);
		}

		//cout << "exp af rad: " << exp(-Radius/RadiusFromCenter) << endl;
		//cout << 1.000002-(RadiusFromCenter/100) << endl;

		
		
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ax);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ay);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*az);
	}

	cout <<"ax = " << ax << "  ay = " << ay << " az= " << az  << endl;
}

void FastEnsemble::InitialiseVelocityHistogram()
{
	//allocating memory to histogram
	VelHistogram = new long double ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		VelHistogram[i] = new long double *[HistNy];
		for(int j = 0; j < HistNy;j++)
			VelHistogram[i][j] = new long double [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				VelHistogram[i][j][k] = 0;
}

void FastEnsemble::InitialiseCountHistogram()
{
	//allocating memory to histogram
	CountHistogram = new long int ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		CountHistogram[i] = new long int *[HistNy];
		for(int j = 0; j < HistNy;j++)
			CountHistogram[i][j] = new long int [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				CountHistogram[i][j][k] = 0;
}

void FastEnsemble::InitialiseHistogram()
{
	//allocating memory to histogram
	histogram = new long int ** [HistNx];
	for(int i = 0; i < HistNx;i++)
	{
		histogram[i] = new long int *[HistNy];
		for(int j = 0; j < HistNy;j++)
			histogram[i][j] = new long int [HistNz];
	}
	for(int i=0; i < HistNx; i++)
		for(int j=0; j < HistNy; j++)
			for(int k=0; k < HistNz; k++)
				histogram[i][j][k] = 0;

}

void FastEnsemble::UpdateHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			histogram[Nx][Ny][Nz]++;
		}
	}
}

void FastEnsemble::MyUpdateVelocityHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			VelHistogram[Nx][Ny][Nz]+= ((long double) pow(ions[i].Velocity(),2));
		}
	}
}

void FastEnsemble::UpdateCountHistogram()
{
	for (int i = 0; i < NumberOfIons; i++)
	{
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			CountHistogram[Nx][Ny][Nz]++;
		}
	}
}

double FastEnsemble::ReturnHist(int i, int j, int k)
{
	return histogram[i][j][k];
}

double FastEnsemble::ReturnVelHist(int i, int j, int k)
{
	return ((double) VelHistogram[i][j][k]);
}

double FastEnsemble::ReturnCountHist(int i, int j, int k)
{
	return ((double) CountHistogram[i][j][k]);
}

void FastEnsemble::SetSteadyStateTemperature(double Val)
{
	SteadyStateTemperature = Val;
}

int FastEnsemble::GetNumberOfIons()
{
	return NumberOfIons;
}

double FastEnsemble::Mass(int N)
{
	return ions[N].GetMass();
}

double FastEnsemble::getRho0(double Vrf){

	return eps0*Vrf*Vrf/(Mass(0)*pow(trap.r0,4)*trap.OmegaRF*trap.OmegaRF);
}

void FastEnsemble::CrystalGenerator(double Vrf, double Vend)
{
	// Calculating pseudo trap frequencies
	double wz2=2*trap.eta*e*Vend/pow(trap.z0,2)/Mass(0); // ERROR in formula from Magnus' thesis! he means potential and not electric potential.
	double wr2=pow(e*Vrf/Mass(0)/trap.OmegaRF,2)/2/pow(trap.r0,4)-trap.eta*e*Vend/Mass(0)/pow(trap.z0,2);

	// Solving plasma model aspect-ratio equation

	double x0 = 0.001;
	double xn = 100.0;
	double eps = 0.0001;
	BisectionRootFinder Roots;
	Roots.SetWsquare(wz2/wr2);

	double alpha = Roots.bisection(x0,xn,eps);//CalculateAspectRatio(wz2/wr2);
	//cout << "alpha " << alpha <<" and ratio is "<< wz2/wr2  << endl;
	//double alpha = CalculateAspectRatio(2);
	//cout << "aspect ratio is " << alpha << '\n';
	/*cout << "out putting ratio(alpha)\n";
	for (double ratio = 0.01; ratio < 2; ratio+=0.01)
		if (ratio != 1)
			cout << CalculateAspectRatio(ratio*ratio) << ';' << '\n';
	 */
	// Calculating number density
	double rho0 = (eps0*Vrf*Vrf)/(Mass(0)*pow(trap.r0,4)*trap.OmegaRF*trap.OmegaRF);

	// Calculating volume of crystal
	double V = ((double) NumberOfIons)/rho0;
	//cout << "volume of crystal is " << V << '\n';

	// Calculating length and radius of crystal
	double L = pow(6*V/(PI*alpha*alpha), ((double) 1) /  ((double) 3));
	bool CrystalNotDone = true;
	while(CrystalNotDone)
	{
		double R = alpha*L/2;

		//cout << "radius is " << R << '\n';
		//cout << "Length is " << L << '\n';
		//cout << "2R/L=" << 2*R/L << '\n';

		// unit cell cube length for bcc structure
		double a = pow( 2 / rho0, ((double) 1) / ((double) 3));

		//cout << "unit cell length is " << a << '\n';

		//cout << "Number of ions " << NumberOfIons << '\n';

		int IonNumber = 0;
		for(int i = -1*((int) ceil(R/a)); i <= ((int) ceil(R/a));i++)
			for(int j = -1*((int) ceil(R/a)); j <= ((int) ceil(R/a));j++)
				for(int k = -1*((int) ceil(L/a/2)); k <= ((int) ceil(L/a/2));k++)
				{
					// placing ion in unit cell "origo" corner and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*((double) i)/R,2) + pow(a*((double) j)/R,2) + pow(a*((double) k)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, ((double) i)*a);
						ions[IonNumber].SetPosition(1, ((double) j)*a);
						ions[IonNumber].SetPosition(2, ((double) k)*a);
						IonNumber++;
					}


					// placing ion in unit cell center and checking if all ions have a position
					if((IonNumber < NumberOfIons) && (pow(a*(((double) i) + 0.5)/R,2) + pow(a*(((double) j) + 0.5)/R,2) + pow(a*(((double) k) + 0.5)*2/L,2) <= 1)) //if ion is inside ellipsoid
					{
						ions[IonNumber].SetPosition(0, (((double) i) + 0.5)*a);
						ions[IonNumber].SetPosition(1, (((double) j) + 0.5)*a);
						ions[IonNumber].SetPosition(2, (((double) k) + 0.5)*a);
						IonNumber++;
					}

				}

		if(IonNumber < NumberOfIons)
		{
			//cout << NumberOfIons - IonNumber << " ion(s) have not been positioned\n";
			L=L*1.005; // increasing crystal length with 0.5 percent
		}
		else
		{
			Radius = R;
			Length =L;
			//cout << "All ions have been positioned\n";
			CrystalNotDone=false;
		}


		//cout << a << " a_WS" <<endl;
	}

	// NORMALT ER DEN NUL! LIGE NU SVARER DEN TIL 1K !
	// setting all ions to zero-velocity
	for (int i = 0; i < NumberOfIons; i++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			double vel = rand() % 10;
			if(vel < 5 ) {
				vel = -1.0;
			}
			else {
				vel = 1.0;
			}
			//ions[i].SetVelocity(dim, vel*5.5838);
			ions[i].SetVelocity(dim, vel* (sqrt( ((3*Kb*trap.StartVelOfIons) /ions[i].m) / 3)));
		}

	}
	cout << Radius << " = Radius of crystal " << Length << " = length of crystal"<< endl; //REMOVE THIS!!
}

double FastEnsemble::Ekin()
{
	double Ekin = 0;
	for(int n = 0; n < NumberOfIons; n++)
		Ekin += ions[n].Ekin();

	return Ekin;
}

double FastEnsemble::Ttot()
{
	return Ekin()/(1.5*NumberOfIons*Kb);
}

double FastEnsemble::Position(int dim, int N)
{
	return ions[N].Position(dim);
}

double FastEnsemble::Velocity(int dim, int N)
{
	return ions[N].Velocity(dim);
}

void FastEnsemble::SavePositionToFile()
{
	FILE *f;
	string name = trap.dirName+"\\MDpos.xyz";
	const char * c = name.c_str();
	f=fopen(c,"w");

	fprintf(f, "%d\n%s\n", NumberOfIons, "Ca");
	for(int N=1; N <= NumberOfIons; N++)
		fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, Position(0, N-1)*1e6, Position(1, N-1)*1e6, Position(2, N-1)*1e6);

	fclose(f);
}


void FastEnsemble::SaveTwoIonPositionToFile()
{
	FILE *f;
	string name = trap.dirName+"\\MDpos.xyz";
	const char * c = name.c_str();
	f=fopen(c,"w");

	fprintf(f, "%d\n%s\n",  NumberOfIons, "Ca");
	for(int N=1; N <= NumberOfIons; N++)
	{
		if (N <= IonOneN) 
		{
			fprintf(f, "%s%d\t%f\t%f\t%f\n ", "Ca", N, Position(0, N-1)*1e6, Position(1, N-1)*1e6, Position(2, N-1)*1e6);
		}
		if (N > IonOneN)
		{
			fprintf(f, "%s%d\t%f\t%f\t%f\n ", "C", IonOneN-N+1, Position(0, N-1)*1e6, Position(1, N-1)*1e6, Position(2, N-1)*1e6);
		}

	}
	fclose(f);
}





void FastEnsemble::SaveIonDataToFile()
{
	// Saves ion data to file
	ofstream Ionfile (trap.dirName+"\\IonData.txt");

	Ionfile << "N \t mass \t x \t y \t z \t Vx \t Vy \t Vz \t Vsec \t norm(V)" << endl; 
	Ionfile << "R: "<< Radius <<" L: " << Length << endl; 
	for(int N = 0; N < GetNumberOfIons(); N++)
	{
		if (Ionfile.is_open())
		{
			Ionfile <<			N	+ 1			<< ",  " <<
								ions[N].GetMass()		<< ",  " <<
								ions[N].Position(0)		<< ",  " <<
								ions[N].Position(1)		<< ",  " <<
								ions[N].Position(2)		<< ",  " <<
								ions[N].Velocity(0)		<< ",  " <<
								ions[N].Velocity(1)		<< ",  " <<
								ions[N].Velocity(2)		<< ",  " <<
								ions[N].Velocity()		<< ",  " <<
								endl;
		}   

	}
	Ionfile.close();
    
}

// For two ion types. 
void FastEnsemble::InitialiseHistogramsForTwoIonTypes()
{
	int SpeciesOfIons = 2;
	histograms = new long double *** [9];
	for (int x = 0; x < 9; x++)
	{
		histograms[x] = new long double ** [HistNx];
		for(int i = 0; i < HistNx;i++)
		{
			histograms[x][i] = new long double *[HistNy];
			for(int j = 0; j < HistNy;j++)
			{
				histograms[x][i][j] = new long double [HistNz];
			}
		}
	}

	for (int x = 0; x<9;x++)
		for(int i=0; i < HistNx; i++)
			for(int j=0; j < HistNy; j++)
				for(int k=0; k < HistNz; k++)
					histograms[x][i][j][k] = 0;
	
}

void FastEnsemble::UpdateVelandCountHistograms()
{ 
	for (int i = 0; i < NumberOfIons; i++)
	{

		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));

		// Fyld i det samlede histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			histograms[7][Nx][Ny][Nz]++;
			histograms[8][Nx][Ny][Nz]+= ((long double) pow(ions[i].Velocity(),2));
		}
		// Fyld i ion 1 histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0 && i < IonOneN)
		{
			histograms[1][Nx][Ny][Nz]++;
			histograms[2][Nx][Ny][Nz]+= ((long double) pow(ions[i].Velocity(),2));
		}
		// Fyld i ion 2 histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0 && i > IonOneN)
		{
			histograms[4][Nx][Ny][Nz]++;
			histograms[5][Nx][Ny][Nz]+= ((long double) pow(ions[i].Velocity(),2));
		}

	}
}

void FastEnsemble::UpdateHistograms()
{ 
	for (int i = 0; i < NumberOfIons; i++)
	{
		
		int Nx = ((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2));
		int Ny = ((int) ((ions[i].Position(1))/PixelToDistance+((double) HistNy)/2));
		int Nz = ((int) ((ions[i].Position(2))/PixelToDistance+((double) HistNz)/2));
		
		// Fyld i det samlede histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0)
		{
			
			histograms[6][Nx][Ny][Nz]++;
			
		}
		// Fyld i ion 1 histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0 && i < IonOneN)
		{
			histograms[0][Nx][Ny][Nz]++;
		}
		// Fyld i ion 2 histogram
		if (Nx < HistNx && Ny < HistNy && Nz < HistNz && Nx > 0 && Ny > 0 && Nz > 0 && i > IonOneN)
		{
			histograms[3][Nx][Ny][Nz]++;
		}
		
	}
}

double FastEnsemble::ReturnFromTwoTypeIonHist(int i, int j, int k, int dim)
{
	return histograms[dim][i][j][k];
}

void FastEnsemble::RescaleVelocityXYZTwoSpecies(double Total_V_x_rms1,double Total_V_y_rms1,double Total_V_z_rms1,double Total_V_x_rms2,double Total_V_y_rms2,double Total_V_z_rms2)
{
	// calculate current temperature for species one.
	double CurrentTempZ = (pow(Total_V_z_rms1,2)*Mass(0))/(Kb);
	double CurrentTempX = (pow(Total_V_x_rms1,2)*Mass(0))/(Kb);
	double CurrentTempY = (pow(Total_V_y_rms1,2)*Mass(0))/(Kb);

	// rescale velocity distribution (Maybe change which fraction of SST is given)
	 
	//cout << "T :" << (pow(TotalV,2) * Mass(0))/(3*Kb) <<		sqrt(pow(CurrentTempZ,2)+pow(CurrentTempX,2) + pow(CurrentTempY,2) )<< endl;
	double T_fraction = SteadyStateTemperature;//(sqrt(3)*SteadyStateTemperature)/3; // This was just SteadyState Before

	double az = sqrt((T_fraction)/CurrentTempZ);
	double ax = sqrt((T_fraction)/CurrentTempX);
	double ay = sqrt((T_fraction)/CurrentTempY);

	//cout << "T's :" << CurrentTempX << " "<< CurrentTempY << " " << CurrentTempZ << endl;
	//cout << ax <<" " << ay << " " << az << endl; 

	double lowerLimit = 0.98;//sqrt(100/102.5); // Hvornår der ikke skal køles mere.
	double upperLimit = 1.02;//sqrt(100/97.5); // Hvornår der ikke skal varmes op mere.



		if (ax > upperLimit)
		{

			ax = 1.000002;

		} 
		else if (ax < lowerLimit)
		{
			ax = 0.996;


		}
		else {
			ax = 1.0;
		}

		if (ay > upperLimit)
		{

			ay = 1.000002;
		}
		else if (ay < lowerLimit)
		{
			ay = 0.996;
		}
		else {

			ay = 1.0;
		}

		if (az > upperLimit)
		{
			az = 1.000002;
		}
		else if (az < lowerLimit)
		{
			az = 0.996;
		}
		else {
			az = 1.0;
		}		

	for(int n = 0; n < IonOneN; n++)
	{

		double RadiusFromCenter = sqrt(pow(ions[n].Position(1),2) +pow(ions[n].Position(2),2) + pow(ions[n].Position(3),2));

		if (ax == 0.996)
		{

			  ax = 0.999-(RadiusFromCenter*100);

		}

		if (ay == 0.996)
		{

			ay = 0.999 - (RadiusFromCenter)*100;
		}

		if (az == 0.996)
		{

			az = 0.999 - (RadiusFromCenter)*100;
		}


		if (ax == 1.000002)
		{
			ax = 1.000002 -(RadiusFromCenter/100);
		}

		if (ay == 1.000002)
		{

			ay = 1.000002 -(RadiusFromCenter/100);
		}

		if (az == 1.000002)
		{

			az = 1.000002 -(RadiusFromCenter/100);
		}
			
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ax);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ay);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*az);
	}

	//--------------------------------------
		// calculate current temperature for species one.
	CurrentTempZ = (pow(Total_V_z_rms2,2)*Mass(IonOneN))/(Kb);
	CurrentTempX = (pow(Total_V_x_rms2,2)*Mass(IonOneN))/(Kb);
	CurrentTempY = (pow(Total_V_y_rms2,2)*Mass(IonOneN))/(Kb);

	// rescale velocity distribution (Maybe change which fraction of SST is given)
	 
	//cout << "T :" << (pow(TotalV,2) * Mass(0))/(3*Kb) <<		sqrt(pow(CurrentTempZ,2)+pow(CurrentTempX,2) + pow(CurrentTempY,2) )<< endl;
	T_fraction = SteadyStateTemperature;//(sqrt(3)*SteadyStateTemperature)/3; // This was just SteadyState Before

	az = sqrt((T_fraction)/CurrentTempZ);
	ax = sqrt((T_fraction)/CurrentTempX);
	ay = sqrt((T_fraction)/CurrentTempY);

	//cout << "T's :" << CurrentTempX << " "<< CurrentTempY << " " << CurrentTempZ << endl;
	//cout << ax <<" " << ay << " " << az << endl; 

	lowerLimit = 0.98;//sqrt(100/102.5); // Hvornår der ikke skal køles mere.
	upperLimit = 1.02;//sqrt(100/97.5); // Hvornår der ikke skal varmes op mere.



		if (ax > upperLimit)
		{

			ax = 1.000002;

		} 
		else if (ax < lowerLimit)
		{
			ax = 0.996;


		}
		else {
			ax = 1.0;
		}

		if (ay > upperLimit)
		{

			ay = 1.000002;
		}
		else if (ay < lowerLimit)
		{
			ay = 0.996;
		}
		else {

			ay = 1.0;
		}

		if (az > upperLimit)
		{
			az = 1.000002;
		}
		else if (az < lowerLimit)
		{
			az = 0.996;
		}
		else {
			az = 1.0;
		}		

	for(int n = IonOneN; n < NumberOfIons; n++)
	{

		double RadiusFromCenter = sqrt(pow(ions[n].Position(1),2) +pow(ions[n].Position(2),2) + pow(ions[n].Position(3),2));

		if (ax == 0.996)
		{

			  ax = 0.999-(RadiusFromCenter*100);

		}

		if (ay == 0.996)
		{

			ay = 0.999 - (RadiusFromCenter)*100;
		}

		if (az == 0.996)
		{

			az = 0.999 - (RadiusFromCenter)*100;
		}


		if (ax == 1.000002)
		{
			ax = 1.000002 -(RadiusFromCenter/100);
		}

		if (ay == 1.000002)
		{

			ay = 1.000002 -(RadiusFromCenter/100);
		}

		if (az == 1.000002)
		{

			az = 1.000002 -(RadiusFromCenter/100);
		}
			
		ions[n].SetVelocity(0, ions[n].Velocity(0)*ax);
		ions[n].SetVelocity(1, ions[n].Velocity(1)*ay);
		ions[n].SetVelocity(2, ions[n].Velocity(2)*az);
	}

}