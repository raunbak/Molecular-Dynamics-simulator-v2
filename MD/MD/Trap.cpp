#include "stdafx.h"
#include "Trap.h"


Trap::Trap(void)
{
}


Trap::~Trap(void)
{
}

Trap::Trap(double Omega, double R0, double Z0, double ETA, double IonVel, std::string dirstr)
{
	OmegaRF = Omega;					
	r0 = R0;						
	z0 = Z0;						
	eta = ETA;		
	StartVelOfIons = IonVel;
	dirName =  dirstr; 


}
