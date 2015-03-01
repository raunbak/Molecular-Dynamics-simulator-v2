#ifndef FASTION_H_
#define FASTION_H_

#pragma once

class Simulation
{
public:
	double RFVoltage;
	double ECVoltage;
	double SimulatedTemperatur;
	double StartVelOfIons;
	int IonOneN;
	int IonOneMass;
	int IonTwoN;
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


	Simulation(void);
	~Simulation(void);
	Simulation(double RFVoltage,double ECVoltage,double SimulatedTemperatur, double StartVelOfIons,
						int IonOneN,int IonOneMass,int IonTwoN,int IonTwoMass,int Timesteps,
						int StartRecordingOfHistogram,int StepsPrRFPeriode,int SizeOfHistogramsX,
						int SizeOfHistogramsY,int SizeOfHistogramsZ,double OmegaRF,double r0,
						double z0,double eta,double binSize,double IonsOneCharge,double IonsTwoCharge);

	Simulation(double RFVoltage,double ECVoltage,double SimulatedTemperatur, double StartVelOfIons,
						int IonN,int IonMass,int Timesteps,
						int StartRecordingOfHistogram,int StepsPrRFPeriode,int SizeOfHistogramsX,
						int SizeOfHistogramsY,int SizeOfHistogramsZ,double OmegaRF,double r0,
						double z0,double eta,double binSize,double IonsCharge);

};

#endif