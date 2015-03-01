#include "Simulation.h"


Simulation::Simulation(void)
{
}


Simulation::~Simulation(void)
{
}

Simulation::Simulation(double RFV,double ECV,double SimT, double StartVel,
						int IonOneCount,int IonOneM,int IonTwoCount,int IonTwoMass,int Timesteps,
						int StartRecordingOfHistogram,int StepsPrRFPeriode,int SizeOfHistogramsX,
						int SizeOfHistogramsY,int SizeOfHistogramsZ,double OmegaRF,double r0,
						double z0,double eta,double binSize,double IonsOneCharge,double IonsTwoCharge)
{
	RFVoltage = RFV;
	ECVoltage = ECV;
	SimulatedTemperatur = SimT;
	StartVelOfIons = StartVel ;
	IonOneN = IonOneCount ;
	IonOneMass = IonOneM;
	IonTwoN = IonTwoCount ;
	int IonTwoMass;
	
	int Timesteps;
	int StartRecordingOfHistogram;
	int StepsPrRFPeriode;
	int SizeOfHistogramsX;
	int SizeOfHistogramsY;
	int SizeOfHistogramsZ;

	double OmegaRF;	//				= strtod(argv[15], 0)*1e6 ; //Convert to Hz from MHz
	double r0;
	double z0;
	double eta;
	double binSize;

	double IonsOneCharge;
	double IonsTwoCharge;
}

Simulation::Simulation(double RFVoltage,double ECVoltage,double SimulatedTemperatur, double StartVelOfIons,
						int IonN,int IonMass,int Timesteps,
						int StartRecordingOfHistogram,int StepsPrRFPeriode,int SizeOfHistogramsX,
						int SizeOfHistogramsY,int SizeOfHistogramsZ,double OmegaRF,double r0,
						double z0,double eta,double binSize,double IonsCharge)
{
}