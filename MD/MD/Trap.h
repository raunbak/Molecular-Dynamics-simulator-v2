#pragma once
class Trap
{



public:
	Trap();
	~Trap();
	Trap(double Omega, double R0, double Z0, double ETA, double IonVel,std::string dirstr);
	double OmegaRF;					
	double r0;						
	double z0;						
	double eta;		
	double StartVelOfIons;
	std::string dirName;


	void setOmegaRF(double value) {
		OmegaRF = value;
	};
	void setRo(double value) {
		r0 = value;
	};

	void setZ0(double value) {
		z0 = value;
	};
	void setEta(double value) {
		eta = value;
	};
	void setStartVelOfIons(double value) {
		StartVelOfIons = value;

	};

	//double getOmegaRF();
	//double getRo();
	//double getZ0();
	//double getEta();
	//int getStepsPrRFPeriode();


};

