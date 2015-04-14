// MolecularDynamics.cpp : Defines the entry point for the console application.
//

// include from making project in VS-10
#include "stdafx.h"
#include <stdio.h>

// From the original MD.cpp
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "FastEnsemble.h"
#include "integrator.h"
#include "constants.h"
#include "Trap.h"

// To save files.
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <stdarg.h>
#include <windows.h>
#include <string>
using namespace std;


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M", &tstruct);

    return buf;
}





// unix main
int main(int argc, char* argv[])
{

	if (argc < 19) 
	{
		cout << "Not enough arguments: " <<  argc << endl;
		cout << " press a button to stop the program" << endl;
		cin.get();
		return 0;
	}

	if (argc == 19 )
	{
	
		 string dirName = "Simulation-"+currentDateTime();
 
		if (CreateDirectory(dirName.c_str() , NULL))
		{
			cout << "Folder for data files created" << endl; 
			// Directory created
		}
		else if (ERROR_ALREADY_EXISTS == GetLastError())
		{
			// Directory already exists
			cout << "Folder exists, please remove or delete to start simulation" << endl;
			cout << "Press any key to close program " << endl;
			cin.get();
			return 0;
		}
		else
		{
			cin.get();
			return 0;
			 // Failed for some other reason
		}




	clock_t tStart = clock();
	cout << "The program is trying to start with one type of ion!" << endl;

	double RFVoltage			= strtod(argv[1], 0);
	double ECVoltage			= strtod(argv[2], 0);
	double SimulatedTemperatur	= strtod(argv[3], 0);
	double StartVelOfIons		= strtod(argv[4], 0);  // Not used yet
	int NumberOfIons			=  atoi(argv[5]);
	int MassOfIons				=  atoi(argv[6]);

	
	int Timesteps					= atoi(argv[7]);
	int StartRecordingOfHistogram	= atoi(argv[8]);
	int StepsPrRFPeriode			= atoi(argv[9]);  
	int SizeOfHistogramsX			= atoi(argv[10]);
	int SizeOfHistogramsY			= atoi(argv[11]);
	int SizeOfHistogramsZ			= atoi(argv[12]); 

	double OmegaRF					= strtod(argv[13], 0)*1e6; 
	double r0						= strtod(argv[14], 0);    
	double z0						= strtod(argv[15], 0);	 
	double eta						= strtod(argv[16], 0);    
	double binSize					= strtod(argv[17], 0);
	double IonsCharge				= strtod(argv[18],0);  // Not used.    
		
	cout << "Allocating memory..\n";
	Trap trap(OmegaRF,r0,z0,eta,StartVelOfIons,dirName);
	//int m1, int n1, int m2, int n2
	
	
	FastEnsemble crystal(MassOfIons,NumberOfIons,IonsCharge,binSize,SizeOfHistogramsX,SizeOfHistogramsY,SizeOfHistogramsZ, trap);
	
	
	cout << "Generating crystal...\n";
	crystal.CrystalGenerator(RFVoltage,ECVoltage);
	crystal.SetSteadyStateTemperature(SimulatedTemperatur);

	cout << " Simulating with " << crystal.GetNumberOfIons() << " ions" << endl;
	// Starting up the simulation 
	cout << "using steps with length = " << 1/(OmegaRF/2/PI)/StepsPrRFPeriode  << "s\n";

	// Create the integrator.
	//TemperatureRescaleLeFrogintegrator(crystal, Timesteps, StartRecordingOfHistogram, RFVoltage, ECVoltage, StepsPrRFPeriode, trap);
	NewTemperatureRescaleLeFrogintegrator(crystal, Timesteps, StartRecordingOfHistogram, RFVoltage, ECVoltage, StepsPrRFPeriode, trap);
	printf("Time taken for simulation: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "Now saving to datafiles\n";
	
	// Saves histogram data to file
	ofstream Histogramfile (dirName+"\\HistogramData.txt");
	Histogramfile << "Bin  i j k" << endl ;

	ofstream VHistfile (dirName+"\\VelocityHistogramData.txt");
	VHistfile << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile (dirName+"\\CountHistogramData.txt");
	CHistfile << "bin  i j k" << endl ;

	for(int i=0; i < SizeOfHistogramsX; i++)
	{
		for(int j=0; j < SizeOfHistogramsY; j++)
		{
			for(int k=0; k < SizeOfHistogramsZ; k++)
			{
				if (Histogramfile.is_open())
				{
					Histogramfile	<< crystal.ReturnHist(i, j, k)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile.is_open())
				{
					VHistfile		<< crystal.ReturnVelHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile.is_open())
				{
					CHistfile		<< crystal.ReturnCountHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
			}
		}
	}

	Histogramfile.close();
	VHistfile.close();
	CHistfile.close();
    
	crystal.SaveIonDataToFile();


	ofstream Configfile (dirName+"\\Configuration.txt");
	Configfile << RFVoltage << " " << ECVoltage << " " << SimulatedTemperatur << " " 
			   << StartVelOfIons << " " << NumberOfIons << " " << MassOfIons << " "
			   << Timesteps << " " << StartRecordingOfHistogram << " " << StepsPrRFPeriode << " "
				<< SizeOfHistogramsX << " " << SizeOfHistogramsY << " " << SizeOfHistogramsZ << " "
				<< OmegaRF << " " << r0 << " " << z0 << " " << eta <<" " << binSize << " " << IonsCharge << endl;  
                                                       
	//crystal.SavePositionToFile();
	
	cout << "Quest completed! I mean data saved!\n";
	cout << "Press a key to end the program\n";



	cin.get();
	


	return 0;
	}
	else {

	cout << "The program is trying to start, with more than two types of ions!" << endl;	
	string dirName = "Simulation-"+currentDateTime();
 
		if (CreateDirectory(dirName.c_str() , NULL))
		{
			cout << "Folder for data files created" << endl; 
			// Directory created
		}
		else if (ERROR_ALREADY_EXISTS == GetLastError())
		{
			// Directory already exists
			cout << "Folder exists, please remove or delete to start simulation" << endl;
			cout << "Press any key to close program " << endl;
			cin.get();
			return 0;
		}
		else
		{
			cout << "Program could not create folder, Press any key to close program " << endl;
			cin.get();
			return 0;
			 // Failed for some other reason
		}

	double RFVoltage			= strtod(argv[1], 0);
	double ECVoltage			= strtod(argv[2], 0);
	double SimulatedTemperatur	= strtod(argv[3], 0);
	double StartVelOfIons		= strtod(argv[4], 0);
	int IonOneN					=  atoi(argv[5]);
	int IonOneMass				=  atoi(argv[6]);
	int IonTwoN					=  atoi(argv[7]);
	int IonTwoMass				=  atoi(argv[8]);
	
	int Timesteps					= atoi(argv[9]);
	int StartRecordingOfHistogram	= atoi(argv[10]);
	int StepsPrRFPeriode			= atoi(argv[11]);
	int SizeOfHistogramsX			= atoi(argv[12]);
	int SizeOfHistogramsY			= atoi(argv[13]);
	int SizeOfHistogramsZ			= atoi(argv[14]);

	double OmegaRF					= strtod(argv[15], 0)*1e6;
	double r0						= strtod(argv[16], 0);
	double z0						= strtod(argv[17], 0);
	double eta						= strtod(argv[18], 0);
	double binSize					= strtod(argv[19], 0);

	double IonsOneCharge			= strtod(argv[20],0);
	double IonsTwoCharge			= strtod(argv[21],0);

	cout << "Allocating memory..\n";
	clock_t tStart = clock();
	Trap trap(OmegaRF,r0,z0,eta,StartVelOfIons,dirName);
	
	
	FastEnsemble crystal(IonOneMass,IonOneN,IonsOneCharge,IonTwoMass,IonTwoN,IonsTwoCharge,binSize,SizeOfHistogramsX,SizeOfHistogramsY,SizeOfHistogramsZ, trap);

	
	cout << "Generating crystal...\n";
	crystal.CrystalGenerator(RFVoltage,ECVoltage);
	crystal.SetSteadyStateTemperature(SimulatedTemperatur);

	TwoIonTypesTemperatureRescaleLeFrogintegrator(crystal,Timesteps,StartRecordingOfHistogram, RFVoltage, ECVoltage, StepsPrRFPeriode, trap);
	
	printf("Time taken for simulation: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "Now saving to datafiles\n";

	// Saves histogram data to file
	ofstream Histogramfile1 (dirName+"\\HistogramDataIonTypeOne.txt");
	Histogramfile1 << "Bin  i j k" << endl ;

	ofstream VHistfile1 (dirName+"\\VelocityHistogramDataTypeOne.txt");
	VHistfile1 << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile1 (dirName+"\\CountHistogramDataTypeOne.txt");
	CHistfile1 << "bin  i j k" << endl ;

		// Saves histogram data to file
	ofstream Histogramfile2 (dirName+"\\HistogramDataIonTypeTwo.txt");
	Histogramfile2 << "Bin  i j k" << endl ;

	ofstream VHistfile2 (dirName+"\\VelocityHistogramDataTypeTwo.txt");
	VHistfile2 << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile2 (dirName+"\\CountHistogramDataTypeTwo.txt");
	CHistfile2 << "bin  i j k" << endl ;

		// Saves histogram data to file
	ofstream Histogramfile3 (dirName+"\\HistogramDataIonTypeBoth.txt");
	Histogramfile3 << "Bin  i j k" << endl ;

	ofstream VHistfile3 (dirName+"\\VelocityHistogramDataTypeBoth.txt");
	VHistfile3 << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile3 (dirName+"\\CountHistogramDataTypeBoth.txt");
	CHistfile3 << "bin  i j k" << endl ;


	cout << "Hist streams created...\n";


	for(int i=0; i < SizeOfHistogramsX; i++)
	{
		for(int j=0; j < SizeOfHistogramsY; j++)
		{
			for(int k=0; k < SizeOfHistogramsZ; k++)
			{
				if (Histogramfile1.is_open())
				{
					Histogramfile1	<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 0)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile1.is_open())
				{
					VHistfile1		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 2)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile1.is_open())
				{
					CHistfile1		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 1)		<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}

				if (Histogramfile2.is_open())
				{
					Histogramfile2	<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 0)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile2.is_open())
				{
					VHistfile2		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 2)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile2.is_open())
				{
					CHistfile2		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 1)		<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}

				if (Histogramfile3.is_open())
				{
					Histogramfile3	<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 0)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile3.is_open())
				{
					VHistfile3		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 2)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile3.is_open())
				{
					CHistfile3		<< crystal.ReturnFromTwoTypeIonHist(i, j, k, 1)		<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}


			}
		}
	}

	Histogramfile1.close();
	VHistfile1.close();
	CHistfile1.close();
		Histogramfile2.close();
	VHistfile2.close();
	CHistfile2.close();
		Histogramfile3.close();
	VHistfile3.close();
	CHistfile3.close();
    
	crystal.SaveIonDataToFile();

	
	ofstream Configfile (dirName+"\\Configuration.txt");
	Configfile << RFVoltage << " " << ECVoltage << " " << SimulatedTemperatur << " " 
			   << StartVelOfIons << " " <<  IonOneN <<  IonOneMass <<  IonTwoN << " " <<  IonTwoMass << " "
			   << Timesteps << " " << StartRecordingOfHistogram << " " << StepsPrRFPeriode << " "
				<< SizeOfHistogramsX << " " << SizeOfHistogramsY << " " << SizeOfHistogramsZ << " "
				<< OmegaRF << " " << r0 << " " << z0 << " " << eta <<" " << binSize << " " << IonsOneCharge << IonsTwoCharge << endl;  
     
	Configfile.close();

	crystal.SavePositionToFile();

	cout << "Quest completed! I mean data saved!\n";
	cout << "Press a key to end the program\n";


		cin.get();
		return 0;
	}
	/*
	clock_t tStart = clock();

	double Vrf = 220;
	double Vend = 18;//3.1; //Kugle = 18. Oval = 3.1
	// HUSK AT TJEKKE starthist når timesteps ændres! (LAV DET AUTO))
	int TimeSteps = 450000;//1850000;//2050000;//100000; //  Scale : 10000 = 0.1ms 
	double Temperature = 0.009335;//0.9335;	 // Scale is 0.01 = 10mK eller 0.001 = 1mK	    

	int StartHistograms = 360000;
	int StepsPrPeriode = 105;
	double PixelToDistance = 0.89e-6;
	int HistNx = 250;
	int HistNy = 250; // x and y should be the same...
	int HistNz = 250;



	cout << "Allocating memory..\n";
	//int m1, int n1, int m2, int n2
	FastEnsemble crystal(40,500,40,500,PixelToDistance,HistNx,HistNy,HistNz);
	
	cout << "Generating crystal...\n";
	crystal.CrystalGenerator(Vrf,Vend);
	crystal.SetSteadyStateTemperature(Temperature);

	cout << crystal.GetNumberOfIons();
	// Starting up the simulation 
	cout << " ions using steps with length = " << 1/(OmegaRF/2/PI)/StepsPrPeriode  << "s\n";

	// Create the integrator.
	TemperatureRescaleLeFrogintegrator(crystal, TimeSteps, StartHistograms, Vrf, Vend, StepsPrPeriode);
	//MADSDynamicTemperatureLeFrogintegrator(crystal, TimeSteps, Vrf, Vend);
	
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "DONE with simulation\n";
	cout << "Now saving to datafiles\n";
	
	// Saves histogram data to file
	ofstream Histogramfile ("HistogramData.txt");
	Histogramfile << "Bin  i j k" << endl ;

	ofstream VHistfile ("VelocityHistogramData.txt");
	VHistfile << "Sum of all velocity in this bin  i \t j \t k" << endl ;

	ofstream CHistfile ("CountHistogramData.txt");
	CHistfile << "bin  i j k" << endl ;

	for(int i=0; i < HistNx; i++)
	{
		for(int j=0; j < HistNy; j++)
		{
			for(int k=0; k < HistNz; k++)
			{
				if (Histogramfile.is_open())
				{
					Histogramfile	<< crystal.ReturnHist(i, j, k)				<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (VHistfile.is_open())
				{
					VHistfile		<< crystal.ReturnVelHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
				if (CHistfile.is_open())
				{
					CHistfile		<< crystal.ReturnCountHist(i, j, k)			<< ",  " 
									<< i										<< ",  " 
									<< j										<< ",  " 
									<< k				
									<< endl;

				}
			}
		}
	}

	Histogramfile.close();
	VHistfile.close();
	CHistfile.close();
    
	cout << "Histograms done" << endl;

	crystal.SaveIonDataToFile();


	cout << "Ion data done" << endl;
	// cleaning up
    //crystal.FreeTemperatureArrays(); // maybe having this in, helps.
	
	cout << "Arrays freed" << endl;

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout << "Quest completed! I mean data saved!\n";
	cout << "Press a key to end the program\n";


	
	cin.get();
	
	


	return 0;
	*/
}
